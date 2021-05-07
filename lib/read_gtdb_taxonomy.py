#!/usr/bin/env python
from _version import __version__
#todo: rename progress files to progess_getgtdb_..... (create own type of progressfile for different sections of pipeline, because download of dbs is only done once (and differently for gtdb and ncbi) but lca etc is done again and again for each analysis
#todo: add a check for the correct wget version (1.19+) by calling "wget --version"
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

silva_server = "https://www.arb-silva.de/fileadmin/silva_databases/current"
silva_source_dict = { "silva_version" : { "url" : "{}/".format(silva_server), "wishlist" : [ "VERSION.txt" ]}, \
					  "silva_taxfiles" : { "url" : "{}/Exports/taxonomy/".format(silva_server), "wishlist" : ["taxmap_slv_lsu_ref_nr_{}.txt.gz", "taxmap_slv_lsu_ref_nr_{}.txt.gz.md5", "taxmap_slv_ssu_ref_nr_{}.txt.gz", "taxmap_slv_ssu_ref_nr_{}.txt.gz.md5"] }, \
					  "silva_fastas" : { "url" : "{}/Exports/".format(silva_server), "wishlist" : ["SILVA_{}_LSURef_NR99_tax_silva.fasta.gz", "SILVA_{}_LSURef_NR99_tax_silva.fasta.gz.md5", "SILVA_{}_SSURef_NR99_tax_silva.fasta.gz", "SILVA_{}_SSURef_NR99_tax_silva.fasta.gz.md5"] } } #currently, silva does not seem to allow recursive downloads based on filename-patterns --> Using this workaround instead. format function will need to replace '{}' with the database version later
#    --> Consider Grepping and filtering only EUkaryote sequences from these --> merge with gtdb dataset OR merge them all (if not too large) and make sure taxonomy is updated!'accordingly!

taxdb_outfilebasename = "gtdb_taxonomy_br.json.gz" #todo: "_br" still stands for "binrefiner". change that!
acc2taxid_outfilebasename = "gtdb_all.accession2taxid.sorted"
lcawalkdb_outfilebasename = "gtdb_lcawalkdb_br.db" #todo: "_br" still stands for "binrefiner". change that!

_progress_steps = {"download": [None, "01a", "02a", "03a"], \
				  "prepare": ["04a", "04b", "05a", "06a", "06b", "06c", "07a"]} #todo: add more steps for individual dbs (instead of concatenating all to one db



import os
import sys
import traceback
import time
import misc
from misc import openfile
import getdb #TODO: obsolete when i integrate this all into getdb.py

def _download_unixwget(sourceurl, pattern=None, targetdir=None): #adds wget >1.19 to dependencies, is probably not the best way to do this, but easier and more stable that ftp/urrlib2, automaticcaly resumes failed downloads etc ... working-solution for now, BUT CONSIDER SWITCHING!
	"""
	Downloads files from a remote server (specified by 'sourceurl') either specified directly in the url or via a filename-pattern specified with 'pattern'.
	'sourceurl' should either pint to a specific remote file, or to a remote folder. In the latter case, 'pattern'needs to be specified, to indicate which files in that folder to download
	
	'Pattern' should be a bash-wildcard pattern or comma-seperated list of such patterns. Pass 'pattern="*"' for downlading all files in the ftp-directory. 
	Pass 'pattern=None' to download only a specific file indicated in sourceurl.
	For downloading only a single file, either specify file in sourceurl and set pattern=None OR just use the specific filename as 'pattern' (IF the server allows recursive downloads/directory listing).
	
	A target directory needs to be specified. This should be a dedicated, temporary folder used only for the downloads of this script. The target directory will be created if it doesn't already exist.  
		
	This function uses the unix tool wget (v.1.19+) for this, because this already has built-in 'retry on broken connection', basic file verifications and even the possibility to continue incomplete downloads.
	An alternative, purely python-based function (using ftplib) may be available in future versions, but this seems good enough (or even more efficient) for now
	
	This function will try to resume incomplete downloads and will not repeat downloads of already existing, complete files. If broken downloads are already present in the target folder, they should be manually deleted before calling this function 
	"""
	#TODO: maybe move this to misc.py?
	import subprocess
	assert targetdir, "\nERROR: A target directory needs to be specified\n" #may allow a default value in the future (possibly "."), but not yet
	wget_basecmd = ["wget", "-nd", "-q", "--tries=20", "--wait=1", "--show-progress", "--progress=bar:force", "-c", "-N", "-P", targetdir]
	wget_batchargs = ["-r", "-np", "-l", "1", "--reject-regex=\?C=", "-A", pattern]
	# "--reject-regex=?C=" is specifically for download from https (but should not interfere with ftp downloads)
	if pattern: #if pattern exists, assume recursive download and resume that sourceurl points to a FOLDER 
		wgetcmd = wget_basecmd + wget_batchargs + [sourceurl]
	else: #otherwise, if pattern==None,  assume that sourceurl points to a file and disable recursive download
		wgetcmd = wget_basecmd + [sourceurl]
	sys.stderr.write("\n" + " ".join(wgetcmd) + "\n")
	sys.stderr.flush()
	wget_proc = subprocess.Popen(wgetcmd, stdout = subprocess.PIPE, stderr = subprocess.PIPE, text = True, universal_newlines=True)
	while True:
		output = wget_proc.stderr.readline()
		if wget_proc.poll() != None:
			break
		if output:
			sys.stderr.write("\r" + output.rstrip("\n")) #find a way to have newlines everytime a new file is downloaded
	#print("finished")
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

def get_cksum(infile): #for gods sake, contrary to the docs on the ftp server, the refseq-checksums are equivalent to cksum results rather than md5-checksums (or even crc32)! Cannot recreate, and not a single predefined module for calculating these in python! giving up and calling cksum for this!
	"""
	AAAAARGH!GODAMIT!ARGH!WHY?!
	"""
	import subprocess
	wget_proc = subprocess.run(["cksum", infile], stdout = subprocess.PIPE, stderr = subprocess.PIPE, text = True, universal_newlines=True)
	outline = wget_proc.stdout
	#print(outline)
	return outline.split()[0]

	
def download_refseq_eukaryote_prots(sourcedict=refseq_dbsource_dict, targetfolder=None):
	"""
	downloads all protein-fastas from the refseq release subfolders specified in source_dict.
	sourcedict should have categories as keys and a dictionary of corresponding remote folder urls and searchpatterns as values
	sourcedict must also have an entry "crc" pointing to the "release-catalog" folder of the refseq-release (which, unfortunately contains CRC-checksums instead of MD5)
	returns a list of the names of all files downloaded if successful, None if not.
	"""
	#TODO: consider the following: Refseq only provides downloads for about one representative per species. probably that is a good thing, because it prevents overinterpreting results from overrepresented species. but possibly I should also consider downloading the sequence data for the other examples too?
	# start ot nested subfunctions
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
		# end of nested subfunctions
		
	#TODO: check if a flag file exists, that indicates that this already ran successfully
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
	_download_unixwget(sourcedict["crc"]["url"], sourcedict["crc"]["pattern"], targetdir=targetfolder) 
	downloaded_crc = glob.glob(targetfolder + "/release*.files.installed")
	assert len(downloaded_crc) == 1, "\nERROR: Expected exactly 1 'release*.files.installed' file in {} after download, but found {}! Be honest, are you messing with me?\n".format(targetfolder, len(downloaded_crc))
	crcfile = downloaded_crc[0]
	#### Now that We have the crc file, download the rest:
	for vireukcat in [ vec for vec in sourcedict if vec != "crc" ]:
		sys.stderr.write("\nNow downloading refseq release category \"{}\"...\n".format(vireukcat))
		returncode = _download_unixwget(sourcedict[vireukcat]["url"], sourcedict[vireukcat]["pattern"], targetdir=targetfolder)
		if returncode != 0:
			sys.stderr.write("\nWARNING: wget returned non-zero returncode '{}' after downloading {} \n".format(returncode, vireukcat))
	download_list = check_crcfile(crcfile, targetfolder, ",".join(["{}*.protein.faa.gz".format(x) for x in sourcedict if x != "crc"])) # todo: correct md5 to crc, to avoid confusion
	#TODO: create a flag file to indicate that this has already run
	return download_list

def download_gtdb_stuff(sourcedict = gtdb_source_dict, targetfolder=None):
	"""
	downloads all protein-fastas from the refseq release subfolders specified in source_dict.
	sourcedict should have categories as keys and a dictionary of corresponding remote folder urls and searchpatterns as values
	returns a dictionary listing the names of all files downloaded for each category/subfolder if successful, None if not.
	"""
	#TODO: check for a flag file, indicating that this already ran successfully
	# start ot nested subfunctions
	def check_gtdbmd5file(md5_file, targetdir, patternstring): # TODO: merge with check_crcfile() to ONE function (not that hard)!
		"""
		as with _download_unixwget(), 'pattern' can be a comma seperated list of patterns (concatenated to a single string)
		raises AssertionException if something is wrong with the downloaded files, returns a list of downloaded filenames in everything is fine
		"""
		#import re
		import fnmatch
		#print("checking for '{}'".format(patternstring))
		#sys.stdout.flush()
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
		#print("returning {}".format(ok_filelist)) 
		return ok_filelist
		# end of nested subfunctions
		
	for gtdbcat in sourcedict:
		sys.stderr.write("\nNow downloading from gtdb: \"{}\"...\n".format(gtdbcat))
		returncode = _download_unixwget(sourcedict[gtdbcat]["url"], sourcedict[gtdbcat]["pattern"], targetdir=targetfolder)
		if returncode != 0:
			sys.stderr.write("\nWARNING: wget returned non-zero returncode '{}' after downloading {} \n".format(returncode, gtdbcat))
	download_dict = { x : check_gtdbmd5file(os.path.join(targetfolder, "MD5SUM"), targetfolder, sourcedict[x]["pattern"]) for x in sourcedict } #todo: need to do something about MD5SUM file ??
	#TODO: create a flag file, indicating that this already ran successfully
	return download_dict

def download_silva_stuff(sourcedict = silva_source_dict, targetfolder=None):
	# start of nested subfunctions
	def check_silvamd5file(filelist, wishlist):
		allisfine = True
		successlist = []
		#print(filelist)
		for f in filelist:
			if os.path.exists(f) and os.path.basename(f) in wishlist:
				if os.path.basename(f) != "VERSION.txt":
					md5file = f + ".md5"
					#calculate md5hash for downloaded file:
					md5ist = calculate_md5hash(f)
					#read hash from md5-checkfile:
					with open(md5file, "r") as m:
						md5soll = m.readline().split()[0]	
					#compare both with eachother:
					if md5soll != md5ist:
						os.remove(f)
						os.remove(md5file)
						allisfine = False
						continue
				successlist.append(f)
			elif os.path.basename(f) in wishlist:
				sys.stderr.write("\nExpcted downloaded file '{}', but could not find it ...\n".format(f))
				allisfine = False
		assert allisfine, "\nError: Download of silva data failed!\nPlease check connection and retry again later\n"
		return successlist
	
	def getsilvaversion(urldict, targetfolder): #yes I know. This could go easier with urllib2 yadayadayada. But for the other stuff wget is better, so ticking with that for now!
		sys.stderr.write("\n Now trying to get current silva release version...\n")
		versionfilename = os.path.join(targetfolder, urldict["wishlist"][0])
		#if os.path.exists(versionfilename):
		#	sys.stderr.write("\nDeleting pre-existing {}\n".format(versionfilename))
		#	os.remove(versionfilename)
		_download_unixwget(urldict["url"] + urldict["wishlist"][0], targetdir=targetfolder)
		#print("donethat")
		
		with open(versionfilename) as versionfile:
			version = versionfile.read().strip()
			#print("wtf")
		return version
	# end of nested subfunctions
	
	print("silvastuff")
	version = getsilvaversion(sourcedict["silva_version"], targetfolder)
	prelim_downloadlist = [os.path.join(targetfolder, sourcedict["silva_version"]["wishlist"][0])]
	wishdict = {silvacat : [ w.format(version) for w in sourcedict[silvacat]["wishlist"] ] for silvacat in sourcedict }
	#print(prelim_downloadlist)
	for silvacat in ["silva_taxfiles", "silva_fastas"]:
		#print(silvacat)
		for wish in wishdict[silvacat]:
			#print(w)
			#print("---- {} ---".format(wish))
			sys.stderr.write("\nNow downloading from silva: \"{}\"...\n".format(wish))
			returncode = _download_unixwget(sourcedict[silvacat]["url"] + wish, pattern = None, targetdir=targetfolder)
			if returncode != 0:
				sys.stderr.write("\nWARNING: wget returned non-zero returncode '{}' after downloading {} \n".format(returncode, wish))
			prelim_downloadlist.append(os.path.join(targetfolder, wish))
			#print(prelim_downloadlist)
	download_dict = { x : check_silvamd5file([df for df in prelim_downloadlist if not df.endswith(".md5")], wishdict[x]) for x in wishdict} # listing all the files that are not md5 files #TODO: something wrong here. lists all the files again in every key
	#TODO: create a flag file, indicating that this already ran successfully
	return download_dict
		
		
def _empty_taxdicts(): #creating basic "pro-forma" enries for Eukaryotes
	"""
	simply returns prefilled "starterdictionaries" for padding GTDB taxonomy with Eukaryote entries
	ranks and levels have been chosen to be wuivalent with the ncbi taxonomy.
	the bogus rank "eukcat" was defined here to simpify and represent the datasets available thorugh the rfseq/release page of the ncbi ftp-server (because downloading eukaryotes through entrez would take AGES AND would be uncompressed. --> better to find datasets on the ftp server
	"""
	taxdict = {"root": { "parent": "root", "rank":0, "taxname" : "root"},\
			   "r__Cellular_organisms": {"parent" : "root", "rank" : 0, "taxname" : "Cellular organisms"}, \
			   "r__Viruses": {"parent" : "root", "rank": 10, "taxname" : "Virus"}, \
			   "d__Eukaryota": {"parent": "r__Cellular_organisms", "rank" : 10, "taxname" : "Eukaryota"}, \
			   "eukcat__fungi": {"parent": "d__Eukaryota", "rank" : -1, "taxname" : "Fungi"}, \
			   "eukcat__vertebrate": {"parent": "d__Eukaryota", "rank" : -1, "taxname" : "Vertebrates"}, \
			   "eukcat__invertebrate": {"parent": "d__Eukaryota", "rank" : -1, "taxname" : "Invertebrates"}, \
			   "eukcat__plant": {"parent": "d__Eukaryota", "rank" : -1, "taxname" : "Plants"}, \
			   "eukcat__protozoa": {"parent": "d__Eukaryota", "rank" : -1, "taxname" : "Protozoa"}, \
			   "eukcat__viral": {"parent": "d__Viruses", "rank" : -1, "taxname" : "Viruses"}, \
			   "eukcat__vertebrate_mammalian": {"parent": "eukcat__vertebrate", "rank" : -1, "taxname" : "Mammalian"}, \
			   "eukcat__vertebrate_other": {"parent": "eukcat__vertebrate", "rank" : -1, "taxname" : "Non-mammalian_vertebrates"}}

	LCA_walktree = {"root": { "level" : 1, "children" : ["r__Cellular_organisms", "r__Viruses"]}, \
					"r__Cellular_organisms": {"level" : 2, "children" : ["d__Bacteria", "d__Archaea", "d__Eukaryota"]}, \
					"r__Viruses": {"level" : 2, "children" : ["eukcat__viral"]}, \
					"d__Eukaryota": {"level" : 3, "children" : ["eukcat__fungi", "eukcat__vertebrate", "eukcat__invertebrate", "eukcat__plant", "eukcat__protozoa"]}, \
					"eukcat__fungi": {"level" : 4, "children" : []}, \
					"eukcat__vertebrate": {"level" : 4, "children" : ["eukcat__vertebrate_mammalian", "eukcat__vertebrate_other"]}, \
					"eukcat__invertebrate": {"level" : 4, "children" : []}, \
					"eukcat__plant": {"level" : 4, "children" : []}, \
					"eukcat__protozoa": {"level" : 4, "children" : []}, \
					"eukcat__viral": {"level" : 4, "children" : []}, \
					"eukcat__vertebrate_mammalian": {"level" : 5, "children" : []}, \
					"eukcat__vertebrate_other": {"level" : 5, "children" : []} } #using only refseq-release unofficial categories here because ncbi refseq selection is not very detailed yet: "Plants" (do they inclue red algae??), "vertebrates" "invertebrates", "fungi", "virus", "protozoa" (what does that even mean in ncbi taxonomy?)
	return taxdict, LCA_walktree

def read_gtdb_taxonomy_from_tsv(infilename, taxdict=None, LCA_walktree=None):#todo let this return a acc2taxid dictionary after all
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
			level = i + 3 # TODO: !!!!! this was originnaly wrongly set at "i+2" leading to problems with the "level" values: "d__Bacteria" and "d__Archaea" were being assigned same level as r_Cellular_organisms == 2, and all subsequent levels also being shifted one too low!!! should be 3! test if correction now works! (lca of eukaryotes and bacteria should yield "r__cellular_organisms" instead of "d__Bacteria"!)

			rank = (i + 1) * 10
			if taxid not in taxdict:
				if i == 0:
					parent = "r__Cellular_organisms" #based on the (safe) assumption that only cellular organisms appear in the gtdb-taxonomy overview files
				else:
					parent = taxlist[i-1]
				taxdict[taxid] = { "parent" : parent, "rank" : rank, "taxname" : taxid[3:] }
			elif i != 0:
				#print(tokens)
				#print("=========")
				assert taxdict[taxid]["parent"] == taxlist[i-1], "ERROR, found contradicting 'parent' entries for taxon {}! \n full taxon line: '{}'\n".format(taxid, tokens[1])
			if taxid not in LCA_walktree:
				if i < len(taxlist) -1:
					children = [taxlist[i+1]]
				else:
					children = []
				LCA_walktree[taxid] = {"level" : level, "children" : children}
			elif i < len(taxlist) -1 and taxlist[i+1] not in LCA_walktree[taxid]["children"]:
				LCA_walktree[taxid]["children"].append(taxlist[i+1])
				#LCA_walktree[taxid]["children"] = list(set(LCA_walktree[taxid]["children"].append(taxlist[i+1]))) #todo: possible alternative solution if still problem with duplicate "children" entries, delete if not necessary anymore
		# ~ if "UBA1103" in line:
			# ~ print("-------")
			# ~ print(line)
			# ~ print(tokens)
			# ~ print(taxlist)
			# ~ print(taxid)
			# ~ print("------")
		acc2taxidfile.write("0\t{}\t{}\t0\n".format(acc, taxlist[-1])) #"dummy" fields with value 0 in columns 1 and 4, to make it look like the original acc2taxid files expected by "_create_sorted_acc2taxid_lookup()" in getdb.py (interesting infos in columns 2&3)
		linecount += 1
		if linecount % 500 == 0:
			sys.stderr.write("\rread {} lines and added {} taxa so far..".format(linecount, len(taxdict)))
			sys.stderr.flush()
	sys.stderr.write("\rread {} lines and added {} taxa so far..\n".format(linecount, len(taxdict)))
	infile.close()
	acc2taxidfile.close()
	return taxdict, LCA_walktree, acc2taxidfilename #this file should be used to crossreference protein-ids and contig ids with taxids in the gtdb reference-fastas later on! #the filename should also be used as input for "_create_sorted_acc2taxid_lookup()" in getdb.py later on
	#todo: should return acc2taxid_dict after all
	
def read_silva_taxonomy_from_tsv(infilename, taxdict, LCA_walktree): #IMPORTANT! ALWAYS PARSE GTSB FIRST, THEN SILVA!
	"""
	assumes that taxdict and LCA_walktree pre-exist from parsing the gtdb files first!

	Originally planned to use the gtdb-to-silva mapping file https://data.ace.uq.edu.au/public/gtdb/data/silva/latest/silva_mapping.tsv
	BUT that file doesn't help with taxa that are exclusive to gtdb OR exclusive to silva!
	also there currently is no mapping file for LSU only for SSU.
	Therefore parsing taxonomy myself the hard way

	there are many, many conflicts between sivla and gtdb and worse many inconsistent species designations between silva ssu and silva LSU!
	The Strategy for for here now:
		- for each silva entry:
			--> for each taxonomic level in the silva-PATH:
				--> check if the taxon designation already exists in the gtdb-taxid_dict --> if Yes, assign acc to that taxon
				Result: in the worst case, some taxa can only be assigned to Bacterium level due to unclear taxon mapping between gtdb and silva
					
	returns a taxonomy dictionary, a LCA_walktree dictionary and a preliminary unsorted acc2taxid-filename


	"""
	# TODO: there are sme unfortunate cases in silva, where the "organism name" field, that basically gives the species and strain designation, contradicts the taxon path! Example: "Lactobacillus rossiae DSM 15814" which, in one case, is a Eukaryote!! --> add a check that verifies the organism name with it's own path, or ignore "organism name" 
		
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
	sys.stderr.write("\nreading silva data from {}\n".format(infilename))
	sys.stderr.flush()
	
	for line in infile:
		tokens =line.strip().split("\t")
		if len(tokens) != 6 or tokens[0] == "primaryAccession": #TODO: not safe to assume this validation check (expect 6 tab-seperated columns per line) will always be valid. this may change
			continue
		acc=".".join(tokens[:3]) # silva taxmap format: column1 = contigname, column2 = start, column3 = stop. --> accession = contigname.start.stop
		if tokens[3].startswith("Eukaryota"):
			taxtokens = [tokens[3].rstrip(";").split(";")[0]] #The eukaryotic taxonomy system is totally different than the bacterial/archaeal one. Since it is already hard enough to consolidate silva ad gtdb taxonomies, for now (to make things easier) ignoring all specific classifications above domain for eukaryotes! (This tool is meant for bacteria/archaea, eukaryotes schould just be recognized as contamination)
		else:
			taxtokens = [t.replace(" ", "_") for t in tokens[3].rstrip(";").split(";") ] #
		slv_species = "s__" + " ".join(tokens[4].split()[:2]) #4th token is silvas species designation. This often has more than 2 space seperated words, while gtdb only has 2. SO only using first two words of each species designation. this species designation is ignored UNLESS it already exists in the gtdb taxonomy
		taxlist = [ rank_prefix_dict[x] + taxtokens[x] for x in range(len(taxtokens)) ] + [slv_species] #sticking the taxonlevel-prefixes in front of the taxon-names (gtdb-style)
		for i in reversed(range(len(taxlist))):			
			taxid = taxlist[i]
			level = i + 2
			rank = (i + 1) * 10
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
	The first option is probaly the fastest/least wasteful because it keeps the lookup-dict as small as possible, but is also very unflexible in case anything changes in the way gtdb stores it's reference files and taxonomy
	the last option seems wasteful, but is probably the most flexible way to go (in case gtdb decides to change their reference protein naming scheme)
	-->choosing the SECOND option for now (compromise bewteen wastefulness and flexibility) and will have to keep an eye on potential changes in the accession-naming at gtdb...
	"""
	pass		

def _concat_fastas(contig_fastalist, outfastahandle, return_headerdict = False, remove_prodigalIDs = False, remove_descriptions = False): #todo: make solution for nucleotide fastas (maybe just create individual blast-dbs because ncbi-blast can handle that)?
	"""
	concatenates the fastafiles given in contig_fastalist
	if return_headers == True, thisreturns dictionary with assemblyIDs as keys and lists of corresponding fasta-sequenceIDs as values. otherwise returns None
	if remove_prodigalIDs == True, it will remove the protein id from the fastaids (based on the regex "_\d+ ")
	(--> naming the proteins simply after the contigs they originated from. additional protein info will remain in the sequence description)
	if remove_descriptions == True, it will remove all description info from the headers (keep only the sequenceID)
	"""
	import re
	assemblyIDpattern = re.compile("GC._\d+\.\d")
	prodigalinpattern = re.compile("(^>\w+\.\d+)(_\d+)(.*)")
	prodigalreplacement = r"\1 \3"
	headerdict = {}
	filecounter = seqcounter = 0
	for cf in contig_fastalist:
		if os.path.isdir(cf):
			continue	
		filebasename = os.path.basename(cf)
		assemblyid = re.search(assemblyIDpattern, filebasename).group(0)
		infile = openfile(cf)
		for line in infile:
			if line.startswith(">"):
				seqcounter += 1
				if remove_descriptions:
					line = line.split()[0]
				if remove_prodigalIDs:
					line = re.sub(prodigalinpattern, prodigalreplacement, line)
				line = ">{}_{}\n".format(assemblyid, line[1:].rstrip()) #prepend assembly id before seqid, because some have nonunique contignames such as "contig_295" (sigh)
				if return_headerdict:
					seqid = line.split()[0][1:]
					if assemblyid in headerdict:
						headerdict[assemblyid].append(seqid)
					else:
						headerdict[assemblyid] = [seqid]
			outfastahandle.write(line)
		infile.close()
		os.remove(cf) #delete file because it is no longer needed #todo: uncomment this! only commented out for test-runs
		filecounter += 1
		if seqcounter % 1000 == 0:
			sys.stderr.write("\rread {} of {} files so far, wrote {} sequences to {}".format(filecounter, len(contig_fastalist), seqcounter, outfastahandle.name))
			sys.stderr.flush()
	sys.stderr.write("\rread {} of {} files so far, wrote {} sequences to {}\n".format(filecounter, len(contig_fastalist), seqcounter, outfastahandle.name))
	sys.stderr.flush()
	return headerdict if return_headerdict else None
				
def gtdb_contignames2taxids(contig_fastalist, acc2taxidinfilelist, acc2taxidoutfilename, outfastahandle):
	"""
	Todo: read acc2taxidfile line by line (taxids for assemblie-ids)
	Todo: open_contigfastas (named after assembly-ids, read contig names), read seqeunces and write to concatenated file
	Todo: assign each contig to taxid based on assembly-id
	write to new acc2taxidfile
	return new acc2taxid-filename
	Note:	sadly the refrence file at gtdb do not follow the same naming scheme for protein and genomic files (or their taxid-lookuptables). 
			the seqids in their taxid-lookuptables as well as the protein files have a "RS_" or "GB_" prefix, but the genomic files do not.
			have to remove that prefix here
	"""
	import re
	assemblyIDpattern = re.compile("GC._\d+\.\d")
	headerdict = _concat_fastas(contig_fastalist, outfastahandle, return_headerdict = True, remove_prodigalIDs = False, remove_descriptions = True)
	import getdb #todo: delete this. is only for debuggung
	getdb.dict2jsonfile(headerdict, "delme_headerdict.json.gz") #todo: delete this. is only for debuggung (or keep it as additional progress saving point
	outaccfile = openfile(acc2taxidoutfilename, "wt")
	linecounter = contigcounter = 0
	reptuple_list = [] #todo: this is only for checking sanity of gtdb downloads (which genomes are considered "representative" and have seqdate provided? which are missing? are all species/genera/phyla covered? which are missing)? comment this out later
	for acc2taxidinfilename in acc2taxidinfilelist:
		inaccfile = openfile(acc2taxidinfilename)
		for line in inaccfile:
			linecounter += 1
			tokens = line.split("\t") #note: split on tab if complete taxname should be used. split on any whitespace, if only "first part" should be used
			seqid = re.search(assemblyIDpattern,tokens[1]).group(0)
			taxid = tokens[2]
			if seqid not in headerdict:
				#reptuple_list.append((seqid, "NOT_representative")) #todo: this is only for checking sanity of gtdb downloads (which genomes are considered "representative" and have seqdate provided? which are missing? are all species/genera/phyla covered? which are missing)? comment this out later
				continue
			#cd temreptuple_list.append((seqid, "YES_REPRESENTATIVE")) #todo: this is only for checking sanity of gtdb downloads (which genomes are considered "representative" and have seqdate provided? which are missing? are all species/genera/phyla covered? which are missing)? comment this out later
			headerlist = headerdict[seqid] 
			outaccfile.write("\n".join(["0\t{}\t{}\t0\n".format(contigname, taxid) for contigname in headerlist]))
			contigcounter += len(headerlist)
			if linecounter % 100 == 0:
				sys.stderr.write("\rassigned {} contigs to {} taxids".format(contigcounter, linecounter))
				sys.stderr.flush()
		inaccfile.close()
	outaccfile.close()
	with open("delmetemp_repcheckfile.tsv", "wt") as repcheckfile: #todo: this is only for checking sanity of gtdb downloads (which genomes are considered "representative" and have seqdate provided? which are missing? are all species/genera/phyla covered? which are missing)? comment this out later
		repcheckfile.write("{}\n".format("\n".join(["{}\t{}".format(x[0], x[1]) for x in reptuple_list]))) #todo: this is only for checking sanity of gtdb downloads (which genomes are considered "representative" and have seqdate provided? which are missing? are all species/genera/phyla covered? which are missing)? comment this out later
	sys.stderr.write("\rassigned {} contigs to {} taxids\n".format(contigcounter, linecounter))
	sys.stderr.flush()	
	return acc2taxidoutfilename

def refseq_contignames2taxids(fastalist, outfasta_filehandle, acc2taxidoutfilename): #todo: rename in refseq_protids2taxids!
	"""
	creates simplified accesison2taxid file for refseq viral and eukaryotic data
	also concatenates the reference sequence files and delets them afterwards
	returns an acc2taxidoutfilename, and keeps the outfasta_filehandle (does not close it)
	"""
	import re
	suffix_pattern = re.compile("\.\d+\.protein.faa.gz")
	acc2taxidfile = openfile(acc2taxidoutfilename, "wt")
	for filename in fastalist:
		euvircat = "eukcat__" + re.sub(suffix_pattern, "", os.path.basename(filename))
		infile = openfile(filename)
		for line in infile:
			if line.startswith(">"):
				acc=line.split()[0][1:]
				acc2taxidfile.write("0\t{}\t{}\t0\n".format(acc, euvircat))
			outfasta_filehandle.write(line)
		infile.close()
		#os.remove(filename) #todo: uncomment this
	acc2taxidfile.close()
	return acc2taxidoutfilename
		
		
		
		
def _download_dbdata_nonncbi(targetdir, progressdump):
	"""
	Downloads GTDB, SILVA and Refseq Eukaryote data, together with taxonomic info, to create a smaller and more condensed reference dataset compared to NCBI-nr
	The resulting taxonomy will be based almost exclusively on the GTDB taxonomy (except for Eukaryotes)
creates files for storing taxdict, LCA_walktree and acc2taxid files and returns the corresponding filenames [TODO: NOPE IT DOESNT]
	if 'continueflag' is set to True, it will assume an existing targetdir to be from an aborted previous run and will try to continue (or restart) from there
	if 'continueflag' is set to False and targetdir already exists, it will return an AssertionError
	returns gtdb_download_dict, refseq_files, silva_download_dict (in that order)
	"""
	sys.stderr.write("\nWARNING! This download may take a LONG time! However, you may be able to resume if it is aborted prematurely...\n")
	assert os.path.abspath(targetdir) != os.getcwd(), "\nERROR: targetdir may not be the current working directory\n"
	assert progressdump["step"] in _progress_steps["download"], "\nError: loaded progress step '{}' is not involved in Download steps {}!\n".format(progressdump["step"], str(_progress_steps["download"]))
	steporder = _progress_steps["download"]

	#step 1: download gtdb
	step = steporder[1]
	laststep = None 
	sys.stderr.write("\nsaved step from previous run {}--\n".format(laststep))
	sys.stderr.write("\n--{}--\n".format(step))
	currentprogressmarker = "progress_step{}.json".format(step)
	if progressdump["step"] == laststep:
		progressdump["step"] = step
		progressdump["gtdb_download_dict"] = download_gtdb_stuff(gtdb_source_dict, targetdir)
		getdb.dict2jsonfile(progressdump, os.path.join(targetdir, currentprogressmarker))
	else:
		sys.stderr.write("-->already did this earlier. skipping this step\n")
		sys.stderr.flush()
		
	#step 2: download refseq-release eukaryotes
	laststep = step
	step = steporder[steporder.index(step) + 1]
	sys.stderr.write("\n--{}--\n".format(step))
	lastprogressmarker = currentprogressmarker 
	currentprogressmarker = "progress_step{}.json".format(step)
	if progressdump["step"] == laststep:
		progressdump["step"] = step
		progressdump["refseq_files"] = download_refseq_eukaryote_prots(refseq_dbsource_dict, targetdir)
		getdb.dict2jsonfile(progressdump, os.path.join(targetdir, currentprogressmarker))
		#os.remove(os.path.join(targetdir, lastprogressmarker)) #todo: uncomment this
	else:
		sys.stderr.write("-->already done this earlier. skipping this step\n")
		sys.stderr.flush()
		
	#step 3: download silva stuff
	laststep = step
	step = steporder[steporder.index(step) + 1]
	sys.stderr.write("\n--{}--\n".format(step))
	lastprogressmarker = currentprogressmarker 
	currentprogressmarker = "progress_step{}.json".format(step)
	if progressdump["step"] == laststep:
		progressdump["step"] = step
		progressdump["silva_download_dict"] = download_silva_stuff(silva_source_dict, targetdir)
		getdb.dict2jsonfile(progressdump, os.path.join(targetdir, currentprogressmarker))
		#os.remove(os.path.join(targetdir, lastprogressmarker)) #todo: uncomment this
	else:
		sys.stderr.write("-->already done this earlier. skipping this step\n")
		sys.stderr.flush()

	return progressdump

def _check_progressmarker(targetdir):
	"""
	checks if a progressmarkerfile (aka progressdumpfile) exists in the targetdir.
	If such a progressmarkerfile exists, the current progress will be read from it and returned as dictionary.
	otherwise it returns an empty progressdump_dict
	"""
	import re
	import fnmatch
	if not os.path.isdir(targetdir):
		print("targetdir '{}' does not exist yet --> this is a new run!".format(targetdir)) 
		return {"step" : None}
	markerfile = markertag = None
	markerpattern="progress_step*.json*"
	#markertagsubpattern = re.compile("step(\d+[a-z]*)\.")
	markerfilelist = sorted([ m for m in os.listdir(targetdir) if fnmatch.fnmatch(m, os.path.join(markerpattern)) ])
	print("found the following markerfiles: {}".format(markerfilelist))
	if len(markerfilelist) > 0:
		if len(markerfilelist) > 1:
			sys.stderr.write("\nwarning: found multiple progress-markers: {}. Expected only one, so will be using only the latest one: '{}'. if this is incorrect, please delete progressmarkers and run again\n".format(str(markerfilelist), markerfilelist[-1]))
		markerfile = markerfilelist[-1] #if for unexpected reasons multiple progressmarkers exist, accept only the latest one
		#markertag = re.search(markertagsubpattern, os.path.basename(markerfile)).group(1)
		print("highest progressmarker is: '{}'".format(markerfile))
	sys.stdout.flush()
	return getdb.jsonfile2dict(os.path.join(targetdir, markerfile)) if markerfile else {"step" : None}
	
	
	
def _prepare_dbdata_nonncbi(targetdir, progressdump): #Todo:add continueflag argument
	import blasthandler
	import getdb #todo: for using dict2jsonthis will be obsolete, when this is moved there #edit: no it won't. keeping download for gttdb and ncbi data seperate. common stuff goes to misc or getdb
	import time
	def step4a():
		sstart = time.time()
		progressdump["acc2taxidfilelist"] = []
		progressdump["taxdict"], progressdump["lca_walktree"] = None, None
		for gtdbtax in gtdb_download_dict["gtdb_taxfiles"]:
			progressdump["taxdict"], progressdump["lca_walktree"], acc2taxidfile = read_gtdb_taxonomy_from_tsv(gtdbtax, progressdump["taxdict"], progressdump["lca_walktree"]) #todo: change read_gtdb_taxonomy_from_tsv() so that it accepts a filehandle to write acc2taxid-lookup to
			progressdump["acc2taxidfilelist"].append(acc2taxidfile)
		send = time.time()
		print("this took {} seconds".format(send -sstart))
		
	def step4b():#todo take progressdump as input, return progressdump as output
		#TODO: create seperate concat_fastas (and seperate diamond-dbs) for gtdb and refseq (maybe even for each refseq-vireukat category) 
		sstart = time.time()
		assert len(gtdb_download_dict["gtdb_fastas"]) == 2, "\nERROR: number of downloaded gtdb tars different than expected: should be 2 but is {}\n".format(len(gtdb_download_dict)) #Notifying myself, so that i don't forget to adapt this if i ever choose to download more or less reference categories from gtdb
		gtdb_genomefiles = misc.untar(gtdb_download_dict["gtdb_fastas"][0], targetdir, removetar=True) 
		gtdb_proteinfiles = misc.untar(gtdb_download_dict["gtdb_fastas"][1], targetdir, removetar=True) 
		progressdump["concatprotfasta"] = os.path.join(targetdir, "concat_refprot.faa")
		progressdump["concatgenomefasta"] = os.path.join(targetdir, "concat_refgenomes.fasta")
		concatgenomefastahandle = openfile(progressdump["concatgenomefasta"], "wt")
		concatprotfastahandle = openfile(progressdump["concatprotfasta"], "wt")
		progressdump["acc2taxidfilelist"].append(gtdb_contignames2taxids(gtdb_genomefiles, progressdump["acc2taxidfilelist"], os.path.join(targetdir, "temp_acc2taxid_gtdb.acc2taxid.gz"), concatgenomefastahandle))
		concatgenomefastahandle.close()
		_concat_fastas(gtdb_proteinfiles, concatprotfastahandle, return_headerdict = False, remove_prodigalIDs = True, remove_descriptions = False) #does not return an acc2taxidfile ,because that is already covered by the contig-accessions
		send = time.time()
		print("this took {} seconds".format(send -sstart))	
		return concatgenomefastahandle, concatprotfastahandle 
	
	def step5a(): #todo take progressdump as input, return progressdump as output
		sstart = time.time()
		for f in silva_download_dict["silva_taxfiles"]:
			sys.stderr.write("reading {}\n".format(f)) #todo: debugging message
			sys.stderr.flush()
			taxdict, LCA_walktree, filename = read_silva_taxonomy_from_tsv(f, progressdump["taxdict"], progressdump["lca_walktree"])
			progressdump["acc2taxidfilelist"].append(filename)
		send = time.time()
		print("this took {} seconds".format(send -sstart))
		return taxdict, LCA_walktree, progressdump["acc2taxidfilelist"] #todo: instead of returng these variables, just put them into pgrogressdump and return THAT


	def step6a(concatprotfastahandle):
		sstart = time.time()
		progressdump["acc2taxidfilelist"].append(refseq_contignames2taxids(refseq_files, concatprotfastahandle, os.path.join(targetdir, "temp_acc2taxid_refseq.acc2taxid.gz")))
		send = time.time()
		print("this took {} seconds".format(send -sstart))

	def step7a():
		sstart = time.time() #todo: all these time.time() calls are just for debugging. delete later
		progressdump["finalacc2taxidfile"] = getdb._create_sorted_acc2taxid_lookup(progressdump["acc2taxidfilelist"], os.path.join(targetdir, acc2taxid_outfilebasename))
		#todo: when successful, delete files in progressdump["acc2taxidfilelist"] (this already done in getdb._create_sorted_acc2taxid_lookup() but need to uncomment that line
		progressdump["acc2taxidfilelist"] = []
		send = time.time()
		print("this took {} seconds".format(send -sstart))
				
	#####end of nested subfunctions 
	#Todo: Make this more elegant if there is any time for that, later
	#step 4: get gtdb_taxonomy
	#TODO: track all progress n form of progressdump. That will make saving and restoring progress much more straightforward!
	#TODO: make seperate diamond databases for GTDB-bacteria, GTDB-Archaea and Refseq-eukaryotes. Maybe even a sperate one for each eukaryote-category! Test if parallel blasts against these are faster than blasting against a single large db 
	steporder = _progress_steps["prepare"]
	laststep = progressdump["step"]
	assert laststep <= max(steporder), "\nError: last loaded progress step '{}' comes AFTER the preparation steps {}!\n".format(laststep, str(steporder))
	#print(progressdump.keys())
	gtdb_download_dict, refseq_files, silva_download_dict = progressdump["gtdb_download_dict"], progressdump["refseq_files"], progressdump["silva_download_dict"]
	lastprogressmarker = "progress_step{}.json".format(laststep)
	step = steporder[0]
	sys.stderr.write("\n--{}--\n".format(step))
	currentprogressmarker = "progress_step{}.json".format(step)
	#print("{} == {} ? --> {}".format(laststep, _progress_steps["download"][-1], laststep == _progress_steps["download"][-1]))
	sys.stdout.flush() # todo: only for debugging
	#step 04a
	if laststep == _progress_steps["download"][-1]:
		print("doing step '{}'".format(step))
		progressdump["step"] = step
		step4a() ##step 4a: read gtdb_taxonomy files (which use the assembly/genome_accessions as identifiers). NOTE: preliminary refseq taxonomy (knows only domain-level for Eujaryota and Viruses) autmatically added to taxodnomy_dict at this point
		getdb.dict2jsonfile(progressdump, os.path.join(targetdir, currentprogressmarker))
		#os.remove(os.path.join(targetdir, lastprogressmarker)) #todo: uncomment this
	else:
		sys.stderr.write("-->already done this earlier. skipping this step\n")
		sys.stderr.flush()
		
	lastprogressmarker, laststep = currentprogressmarker, step
	step = steporder[steporder.index(step) + 1]
	sys.stderr.write("\n--{}--\n".format(step))
	currentprogressmarker = "progress_step{}.json".format(step)
 	#step 04b
	if progressdump["step"] == laststep:
		print("doing step '{}'".format(step))
		sys.stdout.flush() # todo: only for debugging
		progressdump["step"] = step
		concatgenomefastahandle, concatprotfastahandle = step4b() ##step 4b: uncompress gtdb-tar-files and read contig names for each genome-file --> assign contigs to taxids and write to acc2taxidfile	
		getdb.dict2jsonfile(progressdump, os.path.join(targetdir, currentprogressmarker))
		#os.remove(os.path.join(targetdir, lastprogressmarker)) #todo: uncomment this
	else:
		#concatgenomefastahandle =  openfile(progressdump["concatgenomefasta"], "at")
		concatprotfastahandle = openfile(progressdump["concatprotfasta"], "at")
		sys.stderr.write("-->already done this earlier. skipping this step\n")
		sys.stderr.flush()
		
	#step5: get silva taxonomy
	lastprogressmarker, laststep = currentprogressmarker, step 
	step = steporder[steporder.index(step) + 1]
	sys.stderr.write("\n--{}--\n".format(step))
	#step 05a
	if progressdump["step"] == laststep:
		print("doing step '{}'".format(step))
		sys.stdout.flush()
		sys.stdout.flush() # todo: only for debugging
		currentprogressmarker = "progress_step{}.json".format(step)
		progressdump["step"] = step
		print("adding silva...")
		sys.stdout.flush()
		progressdump["taxdict"], progressdump["lca_walktree"], progressdump["acc2taxidfilelist"] = step5a() ##step 5a: read silva_taxonomy files and add this info to taxdict and LCA_walktree. taxdict and lca_walktree should be complete at this point 
		print("writing taxdict file...")
		sys.stdout.flush()		
		taxdict_file = getdb.dict2jsonfile(progressdump["taxdict"], os.path.join(targetdir, taxdb_outfilebasename)) #this is the final db db file for this
		progressdump["taxdict_file"] = taxdict_file
		print("writing walktree file...")
		sys.stdout.flush()	
		progressdump["lca_walktree_file"] = getdb.dict2jsonfile(progressdump["lca_walktree"], os.path.join(targetdir, "lca_walktreedict.json.gz")) #for saving progress. still needs to be turned to actual lcawalkdb
		print("building walktree...")
		sys.stdout.flush()
		progressdump["lcawalkdb_file"] = getdb.build_lca_db(progressdump["lca_walktree"], os.path.join(targetdir, lcawalkdb_outfilebasename), "root")
		progressdump["lca_walktree"] = progressdump["taxdict"] = "" #saving memory
		getdb.dict2jsonfile(progressdump, os.path.join(targetdir, currentprogressmarker))
		#os.remove(os.path.join(targetdir, lastprogressmarker)) #todo: uncomment this
	else:
		#concatgenomefastahandle =  openfile(progressdump["concatgenomefasta"], "at")
		concatprotfastahandle = openfile(progressdump["concatprotfasta"], "at")
		sys.stderr.write("-->already done this earlier. skipping this step\n")
		sys.stderr.flush()
			
	#step 6: cocatenate all reference datasets, create blastdbs, and delete associated fastafiles
	lastprogressmarker, laststep = currentprogressmarker, step 
	step = steporder[steporder.index(step) + 1]
	sys.stderr.write("\n--{}--\n".format(step))
	currentprogressmarker = "progress_step{}.json".format(step)
	#step 06a
	if progressdump["step"] == laststep:
		progressdump["step"] = step
		step6a(concatprotfastahandle) ##step6a: concatenate all refseq protein reference-fasta files (to outhandle concatfasta), keep outhandle (concatfasta) and return acc2taxidfile(name)
		concatprotfastahandle.close()
		getdb.dict2jsonfile(progressdump, os.path.join(targetdir, currentprogressmarker))
		#os.remove(os.path.join(targetdir, lastprogressmarker)) #todo: uncomment this
	else:
		#concatgenomefastahandle =  openfile(progressdump["concatgenomefasta"], "at")
		sys.stderr.write("-->already done this earlier. skipping this step\n")
		sys.stderr.flush()
	
	lastprogressmarker, laststep = currentprogressmarker, step
	step = steporder[steporder.index(step) + 1]
	sys.stderr.write("\n--{}--\n".format(step))
	currentprogressmarker = "progress_step{}.json".format(step)
	#step 06b
	print("HERE ShOULD COME STEP 6b!")
	print(step)
	print(currentprogressmarker)
	if progressdump["step"] == laststep:
		progressdump["step"] = step
		progressdump["outprotdbname"] = os.path.join(targetdir, "gtdbplus_protdb")
		blasthandler.make_diamond_db(progressdump["concatprotfasta"], progressdump["outprotdbname"])  ##step6b: create protein diamond-db
		getdb.dict2jsonfile(progressdump, os.path.join(targetdir, currentprogressmarker))	
		os.remove(os.path.join(targetdir, lastprogressmarker)) #todo: uncomment this
	else:
		#concatgenomefastahandle =  openfile(progressdump["concatgenomefasta"], "at") #todo: make concatenated gtdb genomes-fasta optional (just using 1S/18S data takes up much less space)
		sys.stderr.write("-->already done this earlier. skipping this step\n")
		sys.stderr.flush()	
			
	lastprogressmarker, laststep = currentprogressmarker, step
	step = steporder[steporder.index(step) + 1]
	sys.stderr.write("\n--{}--\n".format(step))
	currentprogressmarker = "progress_step{}.json".format(step)
	#step 06c #todo: make a step06c() subfunction for this
	if progressdump["step"] == laststep:
		progressdump["step"] = step
		dbfasta_list = [progressdump["concatgenomefasta"]] + silva_download_dict["silva_fastas"] #todo: create a new checkpoint for each created database... until then: UNCOMMENT THIS LINE!
		#dbfasta_list = silva_download_dict["silva_fastas"] #todo: for debugging only. choose line above for full run or implement seperate checkpoints for each...
		progressdump["outnucldblist"] = []
		for dbfasta in dbfasta_list:
			#todo : make the following less dumb and ugly:
			outnucldbname = dbfasta
			for suff in [".gz", ".fasta", ".fa", ".faa", ".fna"]:
				if outnucldbname.endswith(suff):
					outnucldbname = outnucldbname[:-len(suff)]
			progressdump["outnucldblist"].append(blasthandler.make_blast_db_from_gz(dbfasta, outnucldbname, db_type = "nucl")) ##step6d: create nucleotide diamond or blastdb (test which one is faster, blast may be faster here, because there are not so many queries as in the protein-blasts)
		getdb.dict2jsonfile(progressdump, os.path.join(targetdir, currentprogressmarker))
		os.remove(os.path.join(targetdir, lastprogressmarker)) #todo: uncomment this
	else:
		sys.stderr.write("-->already done this earlier. skipping this step\n")
		sys.stderr.flush()
	#TODO: step 7: create combined sorted taxid_lookupfile
	#step 7a
	print("HERE SHOULD COME STEP 7a!")
	lastprogressmarker, laststep = currentprogressmarker, step
	step = steporder[steporder.index(step) + 1]
	sys.stderr.write("\n--{}--\n".format(step))
	currentprogressmarker = "progress_step{}.json".format(step)
	print(step)
	print(currentprogressmarker)
	if progressdump["step"] == laststep:
		progressdump["step"] = step
		step7a()
		getdb.dict2jsonfile(progressdump, os.path.join(targetdir, currentprogressmarker))
	else:
		sys.stderr.write("-->already done this earlier. skipping this step\n")
		sys.stderr.flush()		
	print("finished for now!")
	sys.stdout.flush()
	return progressdump

def getNprepare_dbdata_nonncbi(targetdir, continueflag=False):
	#TODO: pack all data into pgrogessdump
	#  Todo: pass only progressdump to subsequent functions
	#  Todo: have a dictionary of progresssteps and their respecitve functions.
	print(targetdir)
	progressdump = _check_progressmarker(targetdir)
	#steps1-3:
	print("HI THERE")
	#print(progressdump)
	print(progressdump["step"])
	if progressdump["step"] in _progress_steps["download"][:-1]:
		print("download")
		progressdump = _download_dbdata_nonncbi(targetdir, progressdump)
	#steps 4-6
	if progressdump["step"] in [_progress_steps["download"][-1]] + _progress_steps["prepare"][:-1]:
		print("prepare")
		progressdump = _prepare_dbdata_nonncbi(targetdir, progressdump)
		
## test functions start here
def test3():
	getNprepare_dbdata_nonncbi("tempdir", continueflag=True)
	
def maintest():
	test3()

def main():
	import argparse
	myparser = argparse.ArgumentParser(prog=os.path.basename(sys.argv[0]), description = "Download and prepare GTDB- and SILVA-data for MDMcleaner")
	myparser.add_argument("-o", "--outdir", action = "store", dest = "outdir", default = "mdmcleaner_gtdb-data", help = "target directory for gtdb-/silva-data. may not be the current working directory. Default = 'mdmcleaner_gtdb-data'")
	myparser.add_argument("-c", "--continue", action = "store_true", dest = "continue_from_checkpoint", default = False, help = "continue from last progress (if targetdir already exists). Default = False") 
	args = myparser.parse_args()
	getNprepare_dbdata_nonncbi(args.outdir, args.continue_from_checkpoint)

if __name__ == '__main__':
#maintest()	#todo: remove this test-funtion
	main()
