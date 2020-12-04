#!/usr/bin/env python
infile_arch="ar122_taxonomy.tsv"
infile_bac="bac120_taxonomy.tsv"


"""
suggestion 1 for eukaryote part of taxtree:
virus: domain=Viruses
fungi: domain=Eukaryota; kingdom (between)= fungi;
vertebrates: domain=Eukaryota; kingdom (between)=Metazoa; phylum=Chordata
invertevrates: domain=Eukaryota; kingdom (between)=Metazoa
plants: domain=Eukaryota; kingdom (between)= viridiplantae ??, 
protozoa: domain=Eukaryota, kindom (netween) = "wt"f???

suggestion 2 for eukaryote part of taxtree:
virus: domain=Viruses
fungi: domain=Eukaryota; euk_category (eukcat)= fungi;
vertebrates: domain=Eukaryota; euk_category (eukcat)= vertebrates;
invertevrates: domain=Eukaryota; euk_category (eukcat)=invertebrates
plants: domain=Eukaryota; euk_category (eukcat) = plants
protozoa: domain=Eukaryota; euk_category (eukcat) = protozoa


basic structure of lca dicts
taxdict[taxid] = { "parent" : parent, \
					"rank" : rank, \
					"taxname" : taxname }
LCA_walktree[taxid] = { "level" : level, \
						"children" : [] } 
						
						
currently downloading the somewhat arbitrary subsets "fungi", "invertebrate", "plant" "protozoa", vertebrate_mammalian, vertebrate_other and viral from the ftp refseq RELEASE folder (btw: what about "red algae" where are they?)
the folder structure there may well change in the foreseeable future, so i should keep an eye on that
the more constant alternative would be to download all "Known" eukaryote entires from refseq usint entrez. BUT this is highly inefficient/slow (partially also because it downloads uncompressed fastas)
Another (but more work intensive) posibility could be to download all protein files from the "complete" refseq release subfolder and filter only the eukaryotic and viral sequences from there.
So rstri seemed the most efficient compromise for the time being.	

notes for downloading SILVA: 
The database releases are structured into two datasets for each gene: SILVA Parc and SILVA Ref. The Parc datasets comprise the entire SILVA databases for the respective gene, whereas the Ref datasets represent a subset of the Parc comprising only high-quality nearly full-length sequences.				
--> Download REF datasets, not Parc!
--> Download the NR (nonredundant) versions of each dataset!
 https://www.arb-silva.de/fileadmin/silva_databases/current/Exports/SILVA_138.1_LSURef_NR99_tax_silva.fasta.gz (65 Mb)
 https://www.arb-silva.de/fileadmin/silva_databases/current/Exports/SILVA_138.1_SSURef_NR99_tax_silva.fasta.gz (189 Mb)
--> Download taxonomy here:
	https://www.arb-silva.de/fileadmin/silva_databases/current/Exports/taxonomy/taxmap_slv_lsu_ref_nr_138.1.txt.gz
	https://www.arb-silva.de/fileadmin/silva_databases/current/Exports/taxonomy/taxmap_slv_ssu_ref_nr_138.1.txt.gz
	
"""

ncbi_ftp_server = "ftp://ftp.ncbi.nlm.nih.gov"
ftp_adress_refseqrelease = "{}/refseq/release".format(ncbi_ftp_server) #when using ftp lib, ncbi_ftp_server adress needs to be seperate from actual subfolder path. when using wget, both should be supplied as one path/url. using the wget way here 
_refseq_vireukcat_list = ["fungi", "invertebrate", "plant", "protozoa", "vertebrate_mammalian", "vertebrate_other", "viral"] 
refseq_dbsource_dict = { refseq_vireukcat : {"url":"{}/{}/".format(ftp_adress_refseqrelease, refseq_vireukcat), "pattern" : "*.protein.faa.gz"} for refseq_vireukcat in _refseq_vireukcat_list }
refseq_dbsource_dict["crc"] = { "url": "{}/release-catalog/".format(ftp_adress_refseqrelease), "pattern" : "release*.files.installed" } #TODO: add eukaryotic 18S/28S + ITS seqeunces to this! either silva or ncbi?


gtdb_server = "https://data.ace.uq.edu.au/public/gtdb/data/releases/latest"
gtdb_source_dict = { "gtdb_taxfiles" : { "url": "{}/".format(gtdb_server), "pattern" : "*_taxonomy.tsv,MD5SUM" }, \
					 "gtdb_fastas" : { "url": "{}/genomic_files_reps".format(gtdb_server), "pattern" : "gtdb_genomes_reps.tar.gz,gtdb_proteins_aa_reps.tar.gz" }, \
					 "gtdb_vs_ncbi_lookup" : { "url" : "{}/auxillary_files".format(gtdb_server), "pattern" : "*_vs_*.xlsx" } }

# TODO: individual fastas also downloadable at silva (https://www.arb-silva.de/no_cache/download/archive/current/Exports/), together with taxonomy_lokuptables (https://www.arb-silva.de/no_cache/download/archive/current/Exports/taxonomy/)
#    --> TODO: Download fastas (NOT the truncated alignments "*_trunc.fasta.gz", NOT the aligned Fastas, and NOT the *Parc_*.fastas, but just the *Ref_NR99_tax_silva.fasta.gz) !!!!!
#    --> Grep and filter only EUkaryote sequences from those --> merge with gtdb dataset OR merge them all (if not too large) and make sure taxonomy is updated!'accordingly!
import os
import sys
import traceback
import time
import misc
from misc import openfile


def _batchdownload_unixwget(sourcedirurl, pattern, targetdir=None): #adds wget >1.19 to dependencies, is probably not the best way to do this, but easier and more stable that ftp/urrlib2, automaticcaly resumes failed downloads etc ... working-solution for now, BUT CONSIDER SWITCHING!
	"""
	Downloads all files from a remote server (specified by 'sourcedirurl') with filenames that match with 'pattern'.
	"Pattern' should be a bash-wildcard pattern or comma-seperated list of such patterns. Pass 'pattern="*"' for downlading all files in the ftp-directory
	For downloading only a single file, just use the specific filename as 'pattern'.
	uses the unix tool wget (v.1.19+) for this, because this already has built-in retry on broken connection, basic verification and even the possibility to continue incomplete downloads
	an alternative, purely python-based function (usinf ftplib) may be available in future versions, but this sems good enough for now
	A target directory needs to be specified. This should be a dedicated, temporary folder used only for the downloads of this script. The target directory will be created if it doesn't already exist.  
	This function assumes that only one wget-download is occuring at a given time (based on the existance of a ".listing" file in the target-directory) and could interfere with other parallel downloads
	It will try to resume incomplete downloads and will not repeat downloads of already existing, complete files. If broken downloads are already present in the target folder, they should be manually deleted before calling this function
	if a ".listing" file still exists in the target from a previous aborted download, but downloads should be resumed anyway, pass a 'force=True' argument to indicate that it can be ignored  
	"""
	import subprocess
	assert targetdir, "\nERROR: A target directory needs to be specified\n"
	wgetcmd = ["wget", "-nd", "-np", "-q", "--tries=20", "--wait=1", "--reject-regex=\?C=", "--show-progress", "--progress=bar:force", "-c", "-r", "-l", "1", "-P", targetdir, "-A", pattern, sourcedirurl ]
	# "--reject-regex=?C=" is specifically for download from https (but should not interfere with ftp downloads)
	print(" ".join(wgetcmd))
	wget_proc = subprocess.Popen(wgetcmd, stdout = subprocess.PIPE, stderr = subprocess.PIPE, text = True, universal_newlines=True)
	while True:
		output = wget_proc.stderr.readline()
		if wget_proc.poll() != None:
			break
		if output:
			sys.stderr.write("\r" + output.rstrip("\n")) #find a way to have newlines everytime a new file is downloaded
	print("finished")
	return wget_proc.poll()


def calculate_md5hash(infile): # TODO: probably move to "misc.py"?
	"""
	calculate md5-hash of a binary file
	"""
	import hashlib #for comparing md5-checksums of downloaded databases
	blocksize = 2**20 #chunks to read file in
	with open(infile, "rb") as f:
		filehash = hashlib.md5()
		while True:
			data = f.read(blocksize)
			if not data:
				break
			filehash.update(data)
	return filehash.hexdigest()	

def calculate_crc32hash(infile): # TODO: probably move to "misc.py"?
	"""
	calculate CRC32-hash of a binary file
	"""
	import zlib
	blocksize = 2**20 #chunks to read file in
	with open(infile, "rb") as f:
		crcvalue = 0
		while True:
			data = f.read(blocksize)
			if not data:
				break
			crcvalue = zlib.crc32(data, crcvalue) & 0xffffffff
	return crcvalue

def get_cksum(infile): #for gods sake, the refseq-checksums are equivalent to cksum results rather than md5-checksums (or even crc32)! Cannot recreate, and not a single predefined module for calculating these in python! giving up and calling cksum for this!
	"""
	AAAAARGH!GODAMIT!ARGH!WHY?!
	"""
	import subprocess
	wget_proc = subprocess.run(["cksum", infile], stdout = subprocess.PIPE, stderr = subprocess.PIPE, text = True, universal_newlines=True)
	outline = wget_proc.stdout
	#print(outline)
	return outline.split()[0]
	

def check_crcfile(crc_file, targetdir, patternstring):  # TODO: merge with check_md5file() to ONE function (not that hard)!
	"""
	as with _download_unixwget(), 'pattern' can be a comma seperated list of patterns (concatenated to a single string)
	raises AssertionException if something is wrong with the downloaded files, returns a list of downloaded filenames in everything is fine
	"""
	import re
	import fnmatch
	patternlist = patternstring.split(",")
	infile = openfile(crc_file)
	allisfine = True
	ok_filelist, bad_filelist, missing_filelist = [], [], []
	for line in infile:
		#print(line)
		tokens = line.strip().split()
		#print(tokens)
		checksum = tokens[0]
		filename = tokens[1]
		for pattern in patternlist:
			if fnmatch.fnmatch(filename, pattern):
				expectedfile = os.path.join(targetdir, filename)
				if not os.path.isfile(expectedfile):
					allisfine = False
					missing_filelist.append(filename)
					break
				actualhash = get_cksum(expectedfile)
				if not str(actualhash) == checksum:
					allisfine = False
					print("BAD FILE: {}! {} != {}".format(expectedfile, actualhash, checksum))
					bad_filelist.append(filename)
					os.remove(expectedfile)
					break
				#print("yay they match! {}! {} != {}".format(expectedfile, actualhash, checksum))
				ok_filelist.append(expectedfile)
	infile.close()
	assert len(ok_filelist) != 0, "\nERROR: None of the expected files appear to have been downloaded. Do you have read/write permissions for the targetfolder?\n"
	assert allisfine, "\nERROR: something went wrong during download: {} expected files are missing, and {} files have mismatching MD5checksums!\n  --> Missing files: {}\n  --> Corrupted files: {}\n".format(len(missing_filelist), len(bad_filelist), ",".join(missing_filelist), ",".join(bad_filelist))
	return ok_filelist			

def check_gtdbmd5file(md5_file, targetdir, patternstring): # TODO: merge with check_crcfile() to ONE function (not that hard)!
	"""
	as with _download_unixwget(), 'pattern' can be a comma seperated list of patterns (concatenated to a single string)
	raises AssertionException if something is wrong with the downloaded files, returns a list of downloaded filenames in everything is fine
	"""
	import re
	import fnmatch
	patternlist = patternstring.split(",")
	infile = openfile(md5_file)
	allisfine = True
	ok_filelist, bad_filelist, missing_filelist = [], [], []
	for line in infile:
		tokens = line.strip().split()
		checksum = tokens[0]
		filename = os.path.basename(tokens[1])
		for pattern in patternlist:
			if fnmatch.fnmatch(filename, pattern):
				expectedfile = os.path.join(targetdir, filename)
				if not os.path.isfile(expectedfile):
					allisfine = False
					missing_filelist.append(filename)
					break
				actualhash = calculate_md5hash(expectedfile)
				if not str(actualhash) == checksum:
					allisfine = False
					bad_filelist.append(filename)
					os.remove(expectedfile)
					break
				#print("yay they match! {}! {} != {}".format(expectedfile, actualhash, checksum))
				ok_filelist.append(expectedfile)
	infile.close()
	assert len(ok_filelist) != 0, "\nERROR: None of the expected files appear to have been downloaded. Do you have read/write permissions for the targetfolder?\n"
	assert allisfine, "\nERROR: something went wrong during download: {} expected files are missing, and {} files have mismatching CRCchecksums!\n  --> Missing files: {}\n  --> Corrupted files: {}\n".format(len(missing_filelist), len(bad_filelist), ",".join(missing_filelist), ",".join(bad_filelist))
	return ok_filelist

	
def download_refseq_eukaryote_prots(sourcedict=refseq_dbsource_dict, targetfolder=None):
	"""
	downloads all protein-fastas from the refseq release subfolders specified in source_dict.
	sourcedict should have categories as keys and a dictionary of corresponding remote folder urls and searchpatterns as values
	sourcedict must also have an entry "crc" pointing to the "release-catalog" folder of the refseq-release (which, unfortunately contains CRC-checksums instead of MD5)
	returns a list of the names of all files downloaded if successful, None if not.
	"""
	import glob
	assert targetfolder and os.path.abspath(targetfolder) != os.getcwd(), "\nERROR: You must provide a temporary targetfolder name!\n Temporary Target folder will be created, if it doesn't already exist.\n Temporary Target folder can NOT be the current working directory, as it will be deleted later!\n"
	temp_download_dict = {}
	#### first get the the crc-checksum file:
	sys.stderr.write("\nFirst downloading CRCchecksum file\n")
		### checking for (and deleting) older crc files (shit happens. first download-attempt may have happened before update, and then resumed after update, leading to two different CRCchecksum files)
	preexisting_crc = glob.glob(targetfolder + "/release*.files.installed") 
	if len(preexisting_crc) != 0:
		sys.stderr.write("\nWARNING: found the following pre-existing refseq-CRCchecksum-files in the targetfolder '{}' : {}.\n These are probably remnants of earlier download-attempts and will be deleted now\n".format(targetfolder, ", ".join(preexisting_crc)))
		for f in preexisting_crc:
			os.remove(f)
		### now downloading CRC
	_batchdownload_unixwget(sourcedict["crc"]["url"], sourcedict["crc"]["pattern"], targetdir=targetfolder) 
	downloaded_crc = glob.glob(targetfolder + "/release*.files.installed")
	assert len(downloaded_crc) == 1, "\nERROR: Expected exactly 1 'release*.files.installed' file in {} after download, but found {}! Be honest, are you messing with me?\n".format(targetfolder, len(downloaded_crc))
	crcfile = downloaded_crc[0]
	#### Now that We have the crc file, download the rest:
	for vireukcat in [ vec for vec in sourcedict if vec != "crc" ]:
		sys.stderr.write("\nNow downloading refseq release category \"{}\"...\n".format(vireukcat))
		returncode = _batchdownload_unixwget(sourcedict[vireukcat]["url"], sourcedict[vireukcat]["pattern"], targetdir=targetfolder)
		if returncode != 0:
			sys.stderr.write("\nWARNING: wget returned non-zero returncode '{}' after downloading {} \n".format(returncode, vireukcat))
	download_list = check_crcfile(crcfile, targetfolder, ",".join(["{}*.protein.faa.gz".format(x) for x in sourcedict if x != "crc"])) # todo: correct md5 to crc, to avoid confusion
	return download_list

def download_gtdb_stuff(sourcedict = gtdb_source_dict, targetfolder=None):
	"""
	downloads all protein-fastas from the refseq release subfolders specified in source_dict.
	sourcedict should have categories as keys and a dictionary of corresponding remote folder urls and searchpatterns as values
	returns a dictionary listing the names of all files downloaded for each category/subfolder if successful, None if not.
	"""
	for gtdbcat in sourcedict:
		sys.stderr.write("\nNow downloading from gtdb: \"{}\"...\n".format(gtdbcat))
		returncode = _batchdownload_unixwget(sourcedict[gtdbcat]["url"], sourcedict[gtdbcat]["pattern"], targetdir=targetfolder)
		if returncode != 0:
			sys.stderr.write("\nWARNING: wget returned non-zero returncode '{}' after downloading {} \n".format(returncode, gtdbcat))
	download_dict = { gtdbcat : check_gtdbmd5file(os.path.join(targetfolder, "MD5SUM"), targetfolder, sourcedict[x]["pattern"]) for x in sourcedict } #todo: need to do something about MD5SUM file ??
	return download_dict
		
def _empty_taxdicts(): #creating basic "pro-forma" enries for Eukaryotes
	"""
	simply returns prefilled "starterdictionaries" for padding GTDB taxonomy with Eukaryote entries
	ranks and levels have been chosen to be wuivalent with the ncbi taxonomy.
	the bogus rank "eukcat" was defined here to simpify and represent the datasets available thorugh the rfseq/release page of the ncbi ftp-server (because downloading eukaryotes through entrez would take AGES AND would be uncompressed. --> better to find datasets on the ftp server
	"""
	taxdict = {"root": { "parent": "root", "rank":0, "taxname" : "root"},\
			   "r__Cellular_organisms": {"parent" : "root", "rank" : 0, "taxname" : "Cellular organisms"}, \
			   "d__Viruses": {"rank": 10, "taxname" : "Virus"}, \
			   "d__Eukaryota": {"parent": "r__Cellular_organisms", "rank" : 10, "taxname" : "Eukaryota"}, \
			   "eukcat__Fungi": {"parent": "d__Eukaryota", "rank" : -1, "taxname" : "Fungi"}, \
			   "eukcat__Vertebrates": {"parent": "d__Eukaryota", "rank" : -1, "taxname" : "Vertebrates"}, \
			   "eukcat__Invertebrates": {"parent": "d__Eukaryota", "rank" : -1, "taxname" : "Invertebrates"}, \
			   "eukcat__Plants": {"parent": "d__Eukaryota", "rank" : -1, "taxname" : "Plants"}, \
			   "eukcat__Protozoa": {"parent": "d__Eukaryota", "rank" : -1, "taxname" : "Protozoa"} }


	LCA_walktree = {"root": { "level" : 1, "children" : ["r__Cellular_organisms", "r__Viruses"]}, \
					"r__Cellular_organisms": {"level" : 2, "children" : ["d__Bacteria", "d__Archaea", "d__Eukaryota"]}, \
					"d__Viruses": {"level" : 2, "children" : []}, \
					"d__Eukaryota": {"level" : 3, "children" : ["eukcat_Fungi", "eukcat_Vertebrates", "eukcat_Invertebrates", "eukcat_Plants", "eukcat_Protozoa"]}, \
					"eukcat_Fungi": {"level" : 4, "children" : []}, \
					"eukcat_Vertebrates": {"level" : 4, "children" : []}, \
					"eukcat_Invertebrates": {"level" : 4, "children" : []}, \
					"eukcat_Plants": {"level" : 4, "children" : []}, \
					"eukcat_Protozoa": {"level" : 4, "children" : []} } #using only refseq-release unofficial categories here because ncbi refseq selection is not very detailed yet: "Plants" (do they inclue red algae??), "vertebrates" "invertebrates", "fungi", "virus", "protozoa" (what does that even mean in ncbi taxonomy?)
	return taxdict, LCA_walktree

def read_gtdb_taxonomy_from_tsv(infilename, taxdict=None, LCA_walktree=None):
	"""
	assumes that taxdict, LCA_walktree either BOTH pre-exist and are passed together, or that NONE of them preexist!
	returns a taxonomy dictionary, a LCA_walktree dictionary and a preliminary unsorted acc2taxid-filename
	"""
	if taxdict == None or LCA_walktree == None:
		assert taxdict == None and LCA_walktree == None, "\nError: either taxdict AND LCA_walktree must be passed, or NONE of them!\n"
		taxdict, LCA_walktree = _empty_taxdicts()
	acc2taxidfilename = os.path.join(os.path.dirname(infilename), "temp_acc2taxid_" + os.path.basename(infilename) + ".acc2taxid.gz") #writing preliminary sorted accession2taxid lookup file. TODO: concatenate this with other sorted acc2taxid-lookupfiles using "_create_sorted_acc2taxid_lookup()" in getdb.py. Using "dummy" fields for columns 1 & 4 (actual info in columns 2&3) to make it compatible with "_create_sorted_acc2taxid_lookup()" in getdb.py
	acc2taxidfile = openfile(acc2taxidfilename, "wt")
	infile=misc.openfile(infilename)
	
	linecount = 0
	for line in infile:
		tokens =line.strip().split("\t")
		if len(tokens) != 2: #TODO: not safe to assume this validation check (expect 2 tab-seperated columns per line) will always be valid. this may change
			continue
		acc=tokens[0]
		taxlist = tokens[1].split(";")
		#tdomain = taxlist[0] #tphylum = taxlist[1] #tclass = taxlist[2] #torder = taxlist[3] #tfamily = taxlist[4] #tgenus = taxlist[5] #tspecies = taxlist[6]
		for i in range(len(taxlist)):
			taxid = taxlist[i]
			#print(taxid)
			level = i + 2
			rank = (i + 1) * 10
			if taxid not in taxdict:
				if i == 0:
					parent = taxid
				else:
					parent = taxlist[i-1]
				taxdict[taxid] = { "parent" : parent, "rank" : rank, "taxname" : taxid[3:] }
			elif i != 0:
				assert taxdict[taxid]["parent"] == taxlist[i-1], "ERROR, found contradicting 'parent' entries for taxon {}! \n full taxon line: '{}'\n".format(taxid, tokens[3])
			if taxid not in LCA_walktree:
				if i < len(taxlist) -1:
					children = [taxlist[i+1]]
				else:
					children = []
				LCA_walktree[taxid] = {"level" : level, "children" : children}
			elif i < len(taxlist) -1 and taxlist[i+1] not in LCA_walktree[taxid]["children"]:
				LCA_walktree[taxid]["children"].append(taxid)
		acc2taxidfile.write("0\t{}\t{}\t0\n".format(acc, taxlist[-1])) #"dummy" fields with value 0 in columns 1 and 4, to make it look like the original acc2taxid files expected by "_create_sorted_acc2taxid_lookup()" in getdb.py (interesting infos in columns 2&3)
		linecount += 1
		if linecount % 500 == 0:
			sys.stderr.write("\rread {} lines and added {} taxa so far..".format(linecount, len(taxdict)))
			sys.stderr.flush()
	sys.stderr.write("\rread {} lines and added {} taxa so far..\n".format(linecount, len(taxdict)))
	infile.close()
	acc2taxidfile.close()
	return taxdict, LCA_walktree, acc2taxidfilename #this file should be used to crossreference protein-ids and contig ids with taxids in the gtdb reference-fastas later on! #the filename should also be used as input for "_create_sorted_acc2taxid_lookup()" in getdb.py later on




def read_silva_taxonomy_from_tsv(infilename, taxdict, LCA_walktree): #IMPORTANT! ALWAYS PARSE GTSB FIRST, THEN SILVA!
	"""
	assumes that taxdict and LCA_walktree pre-exist from parsing the gtdb files first!
	TODO: add a check for this!
	Originally planned to use the gtdb-to-silva mapping file https://data.ace.uq.edu.au/public/gtdb/data/silva/latest/silva_mapping.tsv
	BUT that file doesn't help with taxa that are exclusive to gtdb OR exclusive to silva!
	also there currently is no mapping file for LSU only for SSU.
	Therefore parsing taxonomy myself the hard way
	
	returns a taxonomy dictionary, a LCA_walktree dictionary and a preliminary unsorted acc2taxid-filename
	"""
	"""
	Additional NOTE
	there are many, many conflicts between sivla and gtdb and worse many inconsistent species designations between silva ssu and silva LSU!
	The Strategy for for here now:
		- for each silva entry:
			--> for each taxonomic level in the silva-PATH:
				--> check if the taxon designation already exists in the gtdb-taxid_dict --> if Yes, assign acc to that taxon
				Result: in the worst case, some taxa can only be assigned to Bacterium level due to unclear taxon mapping between gtdb and silva
	"""
		
	rank_prefix_dict = { 0 : "d__",\
						 1 : "p__",\
						 2 : "c__",\
						 3 : "o__",\
						 4 : "f__",\
						 5 : "g__",\
						 6 : "s__" }
	
	acc2taxidfilename = os.path.join(os.path.dirname(infilename), "temp_acc2taxid_" + os.path.basename(infilename) + ".acc2taxid.gz") #writing preliminary sorted accession2taxid lookup file. TODO: concatenate this with other sorted acc2taxid-lookupfiles using "_create_sorted_acc2taxid_lookup()" in getdb.py. Using "dummy" fields for columns 1 & 4 (actual info in columns 2&3) to make it compatible with "_create_sorted_acc2taxid_lookup()" in getdb.py
	acc2taxidfile = openfile(acc2taxidfilename, "wt")		 
	infile=misc.openfile(infilename)
	linecount = 0
	#print("here")
	for line in infile:
		#print(line)
		tokens =line.strip().split("\t")
		#print(tokens)
		if len(tokens) != 6 or tokens[0] == "primaryAccession": #TODO: not safe to assume this validation check (expect 6 tab-seperated columns per line) will always be valid. this may change
			continue
		acc=".".join(tokens[:3]) # silva taxmap format: column1 = contigname, column2 = start, column3 = stop. --> accession = contigname.start.stop
		if tokens[3].startswith("Eukaryota"):
			taxtokens = [tokens[3].rstrip(";").split(";")[0]] #The eukaryotic taxonomy system is totally different than the bacterial/archaeal one. Since it is already hard enough to consolidate silva ad gtdb taxonomies, for now (to make things easier) ignoring all specific classifications above domain for eukaryotes! (This tool is meant for bacteria/archaea, eukaryotes schould just be recognized as contamination)
		else:
			taxtokens = [t.replace(" ", "_") for t in tokens[3].rstrip(";").split(";") ] #
		slv_species = "s__" + " ".join(tokens[4].split()[:2]) #4th token is silvas species designation. This often has more than 2 space seperated words, while gtdb only has 2. SO only using first two words of each species designation. this species designation is ignored UNLESS it already exists in the gtdb taxonomy
		taxlist = [ rank_prefix_dict[x] + taxtokens[x] for x in range(len(taxtokens)) ] + [slv_species] #sticking the taxonlevel-prefixes in front of the taxon-names (gtdb-style)
		#tdomain = taxlist[0] #tphylum = taxlist[1] #tclass = taxlist[2] #torder = taxlist[3] #tfamily = taxlist[4] #tgenus = taxlist[5] #tspecies = taxlist[6]


		# ~ for i in reversed(range(len(taxlist))):
		## The following section was only needed, when planning to add silva-exclusive taxa to the tree.
		## since i now decided to only allow taxa that are already in the gtdb taxonomy, i can ignore this
		# ~ badwords = ["uncultured", "endosymbiont", "Incertae sedis", "problematica", "metagenom"] #slv_taxonomy is full of unhelpful taxon designations such as "metagenome" or "uncultured". Could use silva-Taxids but the taxids make it hard to cross-reference with gtdb (which has some bins without an silva entry). until there is a uniform taxid system shared by gtdb and silva, it seems safer to just skip/ignore such taxons and look for are more specific higher ranking taxon designation		
			# ~ for fword in badwords:
				# ~ if fword in taxlist[i]: #slv_taxonomy is full of unhelpful taxon designations such as "metagenome" or "uncultured". Could use silva-Taxids but the taxids make it hard to cross-reference with gtdb (which has some bins without an silva entry). until there is a uniform taxid system shared by gtdb and silva, it seems safer to just skip/ignore such taxons and look for are more specific higher ranking taxon designation
					# ~ continue
			# ~ taxid = taxlist[i].split()[0] # silva has some special taxon infos such as e.g. "NIVA-CYA 15" found in the genus designation "Planktothrix NIVA-CYA 15" not recognized by gtdb taxonomy. removing the extra parts here
		for i in reversed(range(len(taxlist))):			
			taxid = taxlist[i]
			level = i + 2
			rank = (i + 1) * 10
			#taxid = rank_prefix_dict[rank] + taxid
			if taxid not in taxdict:#decided to avoid conflicts between silva and gtdb by basically ignoring silva taxons not found in gtdb (only using highest common marker). This could be much better,if tehre was a current, upt-to-date gtdb to silva-ssu AND silva-lsu lookuptable! will try to see if such a thing can be generated in the future...
				assert i>0, "\nERROR: unrecognized tacon on domain level in silva??? --> {}\n".format(taxid)
				continue
			assert taxid in LCA_walktree, "\nError: the following taxid is present in taxdict but not in LCA_walktree: {}\n".format(taxid)
			break
		acc2taxidfile.write("0\t{}\t{}\t0\n".format(acc, taxid)) #"dummy" fields with value 0 in columns 1 and 4, to make it look like the original acc2taxid files expected by "_create_sorted_acc2taxid_lookup()" in getdb.py (interesting infos in columns 2&3)
		linecount += 1
		if linecount % 500 == 0:
			sys.stderr.write("\rread {} lines and added {} taxa so far..".format(linecount, len(taxdict)))
			sys.stderr.flush()
	sys.stderr.write("\rread {} lines and added {} taxa so far..\n".format(linecount, len(taxdict)))
	infile.close()
	acc2taxidfile.close()
	return taxdict, LCA_walktree, acc2taxidfilename #the filename should be used as input for "_create_sorted_acc2taxid_lookup()" in getdb.py later on

def make_gtdb_blastdb():
	"""
	make blastdbs out of gtdb_reps + refseq_release + SILVA datasets
	for gtdb_reps: 
		- either assign all proteins of each genome to their respective genome/assembly-id (could be slightly faster)
		- OR: assign the CONTIGS to an taxonomy id and cut of the consecutive prodigal CDS numbering of during parsing
		- OR: assign each protein-ID to an taxonomy ID (seems wasteful but more compatible with potential future changes of the GTDB-data structure)
	the last option seems wasteful, but is probably the best way to go (in case gtdb decides to change their reference protein naming scheme)
	However, choosing the second option for now (compromise bewteen wastefulness and flexibility) and will have to keep an eye on potential changes in the accession-naming at gtdb...
	"""
	pass		
				
## test functions start here
def test_ftpwget():
	sys.stderr.write("\nTESTING WGET WITH VIRAL DATA\n")
	sourcedirurl = "ftp://ftp.ncbi.nlm.nih.gov/refseq/release/viral/"
	pattern = "*.protein.faa.gz"
	targetdir = "deltest2"
	retcode = _ftpdownload_unixwget(sourcedirurl, pattern, targetdir)	
	print("\n")
	print(retcode)
	print("\n")

def test_httpwget():
	sys.stderr.write("\nTESTING WGET WITH gtdb-taxonomy-tables\n")
	sourcedirurl = "https://data.ace.uq.edu.au/public/gtdb/data/releases/latest/"
	pattern = "*.tsv"
	targetdir = "deltesthttps"
	retcode = _httpsdownload_unixwget(sourcedirurl, pattern, targetdir)	
	print("\n")
	print(retcode)
	print("\n")


def test_refseq_download():
	print("testing the download")
	testdict = {}
	testdict["viral"] = refseq_dbsource_dict["viral"]
	testdict["crc"] = refseq_dbsource_dict["crc"]
	print(download_refseq_eukaryote_prots(sourcedict=testdict, targetfolder="deltestwtf"))
	print("\nfinished\n")

def test_gtdb_download():
	print("testing the download")
	testdict = {}
	testdict["gtdb_taxfiles"] = gtdb_source_dict["gtdb_taxfiles"]
	testdict["gtdb_vs_ncbi_lookup"] = gtdb_source_dict["gtdb_vs_ncbi_lookup"]
	print(download_gtdb_stuff(sourcedict=testdict, targetfolder="deltestwtf"))
	print("\nfinished\n")

def test_get_lookup_dicts():
	print("testing the creation of lookup_dicts")
	silva_lsu_lookup = "taxmap_slv_lsu_ref_nr_138.1.txt.gz"
	silva_ssu_lookup = "taxmap_slv_ssu_ref_nr_138.1.txt.gz"
	gtdb_lookup_bact = "bac120_taxonomy.tsv"
	gtdb_lookup_arch = "ar122_taxonomy.tsv"
	taxdict, LCA_walktree = None, None
	outfilelist = []
	for f in [ gtdb_lookup_bact, gtdb_lookup_arch]:
		print("reading {}".format(f))
		taxdict, LCA_walktree, filename = read_gtdb_taxonomy_from_tsv(f, taxdict, LCA_walktree)
		outfilelist.append(filename)
	for f in [silva_ssu_lookup, silva_lsu_lookup]:
		print("reading {}".format(f))
		taxdict, LCA_walktree, filename = read_silva_taxonomy_from_tsv(f, taxdict, LCA_walktree)
		outfilelist.append(filename)
	print("outfiles:")
	print("\n".join(outfilelist))
	
def main():
	#download_refseq_eukaryote_prots(".")
	#test_wget()
	#test_refseq_download()
	#test_gtdb_download()
	test_get_lookup_dicts()
	#test_httpwget()
	
main()	
