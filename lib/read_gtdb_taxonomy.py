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
	

def check_crcfile(crc_file, targetdir, patternstring): #as with _download_unixwget(), 'pattern' can be a comma seperated list of patterns (concatenated to a single string)
	#this will replace the returning of filelists from each download-call 
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
					#os.remove(expectedfile)
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
	sourcedict should have categories as keys and corresponding remote folder urls as values
	sourcedict must also have an entry "crc" pointing to the "release-catalog" folder of the refseq-release (which, unfortunately contains CRC-checksums instead of MD5)
	returns a dictionary listing the names of all files downloaded for each category/subfolder if successful, None if not.
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
	_batchdownload_unixwget(sourcedict["crc"]["url"], sourcedict["crc"]["pattern"], targetdir=targetfolder) #TODO: consider simply removing the .listfile as security check. Not sure what exactly I am protecting with that anyway... 
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
	for gtdbcat in sourcedict:
		sys.stderr.write("\nNow downloading from gtdb: \"{}\"...\n".format(gtdbkcat))
		returncode = _batchdownload_unixwget(sourcedict[gtdb]["url"], sourcedict[gtdbcat]["pattern"], targetdir=targetfolder)
		if returncode != 0:
			sys.stderr.write("\nWARNING: wget returned non-zero returncode '{}' after downloading {} \n".format(returncode, gtdbcat))
		download_dict = { gtdbcat : check_md5file(md5file, targetfolder, ",".join(["{}{}".format(x, sourcedict[x]["pattern"]) for x in sourcedict ])) } #todo: need to do something about MD5SUM file ??
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
	if taxdict == None or LCA_walktree == None:
		taxdict, LCA_walktree = _empty_taxdicts()
	infile=misc.openfile(infilename)
	for line in infile:
		tokens =line.strip().split("\t")
		if len(tokens) != 2:
			continue
		acc=tokens[0]
		taxlist = tokens[1].split(";")
		#tdomain = taxlist[0] #tphylum = taxlist[1] #tclass = taxlist[2] #torder = taxlist[3] #tfamily = taxlist[4] #tgenus = taxlist[5] #tspecies = taxlist[6]
		for i in range(len(taxlist)):
			taxid = taxlist[i]
			level = i + 2
			rank = (i + 1) * 10
			if taxid not in taxdict:
				if i == 0:
					parent = taxid
				else:
					parent = taxlist[i-1]
				taxlist[taxid] = { "parent" : parent, "rank" : rank, "taxname" : taxid[3:] }
			if taxid not in LCA_walktree:
				if i < len(taxlist) -1:
					children = [taxlist[i+1]]
				else:
					children = []
				LCA_walktree[taxid] = {"level" : level, "children" : children}
			elif i < len(taxlist) -1 and taxlist[i+1] not in LCA_walktree[taxid]["children"]:
				LCA_walktree[taxid][children].append(taxid)
				
				
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

def main():
	#download_refseq_eukaryote_prots(".")
	#test_wget()
	test_refseq_download()
	#test_httpwget()
	
main()	
