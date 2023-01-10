#!/usr/bin/env python
from mdmcleaner._version import __version__
import re
#todo: rename progress files to progess_getgtdb_..... (create own type of progressfile for different sections of pipeline, because download of dbs is only done once (and differently for gtdb and ncbi) but lca etc is done again and again for each analysis
#todo: add a check for the correct wget version (1.19+) by calling "wget --version"
#todo: in order to reduce protein-db size: cluster each protein based on 99% identity within each genome, to reduce paralogs. keep only longest copy per genome

zenodo_publication_db = "https://zenodo.org/record/5698995/files/MDMcleanerDB.tar.bz2"
zenodo_publication_db_md5 = "b3862577d25000805cc5830f0b1d50e3"

ncbi_ftp_server = "ftp://ftp.ncbi.nlm.nih.gov"
ftp_adress_refseqrelease = "{}/refseq/release".format(ncbi_ftp_server) #when using ftp lib, ncbi_ftp_server adress needs to be seperate from actual subfolder path. when using wget, both should be supplied as one path/url. using the wget way here 
_refseq_vireukcat_list = ["fungi", "invertebrate", "plant", "protozoa", "vertebrate_mammalian", "vertebrate_other", "viral"] 
refseq_dbsource_dict = { refseq_vireukcat : {"url":"{}/{}/".format(ftp_adress_refseqrelease, refseq_vireukcat), "pattern" : "*.protein.faa.gz"} for refseq_vireukcat in _refseq_vireukcat_list }
refseq_dbsource_dict["crc"] = { "url": "{}/release-catalog/".format(ftp_adress_refseqrelease), "pattern" : "release*.files.installed" } #TODO: add eukaryotic 18S/28S + ITS seqeunces to this! either silva or ncbi?

MD5FILEPATTERN_GTDB = "MD5SUM*"
gtdb_server = "https://data.ace.uq.edu.au/public/gtdb/data/releases/latest"
gtdb_source_dict = { "gtdb_taxfiles" : { "url": "{}/".format(gtdb_server), "pattern" : "{}".format(",".join(["*_taxonomy.tsv", MD5FILEPATTERN_GTDB]))}, \
					 "gtdb_fastas" : { "url": "{}/genomic_files_reps".format(gtdb_server), "pattern" : "gtdb_genomes_reps.tar.gz,gtdb_proteins_aa_reps.tar.gz" }, \
					 "gtdb_vs_ncbi_lookup" : { "url" : "{}/auxillary_files".format(gtdb_server), "pattern" : "*_vs_*.xlsx" } } #todo: remove gtdb_vs_ncbi_lookuptables

silva_server = "https://www.arb-silva.de/fileadmin/silva_databases/current"
ALT_silva_server = "ftp://arb-silva.de/current" # apparently sometimes one or the other of the silva servers is not reachable. therefore always trying both alternately

silva_source_dict = { "silva_version" : { "url" : "{}/".format(silva_server), "alturl" : "{}/".format(ALT_silva_server), "wishlist" : [ "VERSION.txt" ]}, \
					  "silva_taxfiles" : { "url" : "{}/Exports/taxonomy/".format(silva_server), "alturl" : "{}/Exports/taxonomy/".format(ALT_silva_server), "wishlist" : ["taxmap_slv_lsu_ref_nr_{}.txt.gz", "taxmap_slv_lsu_ref_nr_{}.txt.gz.md5", "taxmap_slv_ssu_ref_nr_{}.txt.gz", "taxmap_slv_ssu_ref_nr_{}.txt.gz.md5"] }, \
					  "silva_fastas" : { "url" : "{}/Exports/".format(silva_server), "alturl" : "{}/Exports/".format(ALT_silva_server), "wishlist" : ["SILVA_{}_LSURef_NR99_tax_silva.fasta.gz", "SILVA_{}_LSURef_NR99_tax_silva.fasta.gz.md5", "SILVA_{}_SSURef_NR99_tax_silva.fasta.gz", "SILVA_{}_SSURef_NR99_tax_silva.fasta.gz.md5"] } } #currently, silva does not seem to allow recursive downloads based on filename-patterns --> Using this workaround instead. format function will need to replace '{}' with the database version later


#    --> Consider Grepping and filtering only EUkaryote sequences from these --> merge with gtdb dataset OR merge them all (if not too large) and make sure taxonomy is updated!'accordingly!

taxdb_outfilebasename = "gtdb_taxonomy_br.json.gz" #todo: "_br" still stands for "binrefiner". change that!
acc2taxid_outfilebasename = "gtdb_all.accession2taxid.sorted"
lcawalkdb_outfilebasename = "gtdb_lcawalkdb_br.db" #todo: "_br" still stands for "binrefiner". change that!

_progress_steps = {"download": [None, "01a", "02a", "03a"], \
				  "prepare": ["04a", "04b", "04c", "05a", "06a", "06b", "06c", "07a"], \
				  "finished": ["08a"]} #todo: add more steps for individual dbs (instead of concatenating all to one db

default_settings = { x:x for x in ["makeblastdb", "diamond"]}

import os
import sys
import traceback
import time
from mdmcleaner import misc
from mdmcleaner.misc import openfile
from mdmcleaner import getdb #TODO: obsolete when i integrate this all into getdb.py



def _download_unixwget(sourceurl, pattern=None, targetdir=None, verbose=False): #adds wget >1.19 to dependencies, is probably not the best way to do this, but easier and more stable that ftp/urrlib2, automaticcaly resumes failed downloads etc ... working-solution for now, BUT CONSIDER SWITCHING!
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
	wget_basecmd = ["wget", "-nd", "-q", "--tries=20", "--wait=1", "-c", "-N", "-P", targetdir]
	if verbose:
		wget_basecmd += ["--show-progress", "--progress=bar:force"]
	wget_batchargs = ["-r", "-np", "-l", "1", "--reject-regex=\?C=", "-A", pattern]
	# "--reject-regex=?C=" is specifically for download from https (but should not interfere with ftp downloads)
	if pattern: #if pattern exists, assume recursive download and resume that sourceurl points to a FOLDER 
		wgetcmd = wget_basecmd + wget_batchargs + [sourceurl]
	else: #otherwise, if pattern==None,  assume that sourceurl points to a file and disable recursive download
		wgetcmd = wget_basecmd + [sourceurl]
	if verbose:
		sys.stderr.write("\n" + " ".join(wgetcmd) + "\n")
		sys.stderr.flush()
	wget_proc = subprocess.Popen(wgetcmd, stdout = subprocess.PIPE, stderr = subprocess.PIPE, text = True, universal_newlines=True)
	while True:
		output = wget_proc.stderr.readline()
		if wget_proc.poll() != None:
			break
		if output:
			sys.stderr.write("\r" + output.rstrip("\n")) #todo: find a way to have newlines everytime a new file is downloaded
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

def is_md5(inhash): #ncbi keeps switching between md5 and crc (or actually rather cksum) file-hashes (and does not always explicitely document this) --> need simple check to verify if a hash might be md5.
	pattern = r"([a-fA-F\d]{32})"
	if re.match(pattern, inhash):
		return True
	return False
	
def get_cksum(infile): #Contrary to the docs on the ftp server, the refseq-checksums are equivalent to cksum results rather than md5-checksums (or even crc32)! Need to call cksum for this!
	"""
	calculate cksum-hash of a binary file
	"""
	import subprocess
	wget_proc = subprocess.run(["cksum", infile], stdout = subprocess.PIPE, stderr = subprocess.PIPE, text = True, universal_newlines=True)
	outline = wget_proc.stdout
	return outline.split()[0]

	
def download_refseq_eukaryote_prots(sourcedict=refseq_dbsource_dict, targetfolder=None, verbose=False, maxtries = 3):
	"""
	downloads all protein-fastas from the refseq release subfolders specified in source_dict.
	sourcedict should have categories as keys and a dictionary of corresponding remote folder urls and searchpatterns as values
	sourcedict must also have an entry "crc" pointing to the "release-catalog" folder of the refseq-release (which, unfortunately contains CRC-checksums instead of MD5)
	returns a list of the names of all files downloaded if successful, None if not.
	if crc-checksums do not match it will reattempt downloads as often as speciefied by 'maxtries' (default = 3 times)
	"""
	#TODO: consider the following: Refseq only provides downloads for about one representative per species. probably that is a good thing, because it prevents overinterpreting results from overrepresented species. but possibly I should also consider downloading the sequence data for the other examples too?
	# start ot nested subfunctions
	def check_hashfile(hash_file, targetdir, patternstring):  # TODO: merge with check_md5file() to ONE function (not that hard)!
		"""
		as with _download_unixwget(), 'pattern' can be a comma seperated list of patterns (concatenated to a single string)
		raises AssertionException if something is wrong with the downloaded files, returns a list of downloaded filenames in everything is fine
		"""
		import re
		import fnmatch
		hashismd5=False # ncbi seems to keep switching beween file hash types. this now verifies is hash is actually md5...
		patternlist = patternstring.split(",")
		infile = openfile(hash_file)
		allisfine = True
		ok_filelist, bad_filelist, missing_filelist = [], [], []
		for line in infile:
			tokens = line.strip().split()
			checksum = tokens[0]
			if is_md5(checksum):
				hashismd5 = True
			filename = tokens[1]
			for pattern in patternlist:
				if fnmatch.fnmatch(filename, pattern):
					expectedfile = os.path.join(targetdir, filename)
					if not os.path.isfile(expectedfile):
						allisfine = False
						missing_filelist.append(filename)
						break
					if hashismd5:
						actualhash = calculate_md5hash(expectedfile)
					else:
						actualhash = get_cksum(expectedfile)
					if not str(actualhash) == checksum:
						allisfine = False
						sys.stderr.write("\n\tERROR FILE NOT MATCHING CHECKSUM: {}! {} != {}\n".format(expectedfile, actualhash, checksum))
						bad_filelist.append(filename)
						os.remove(expectedfile)
						break
					ok_filelist.append(expectedfile)
		infile.close()
		assert len(ok_filelist) != 0, "\nERROR: None of the expected files appear to have been downloaded. Do you have read/write permissions for the targetfolder?\n"
		if not allisfine:
			sys.stderr.write("\nERROR: something went wrong during download: {} expected files are missing, and {} files have mismatching MD5checksums!\n  --> Missing files: {}\n  --> Corrupted files: {}\nPlease try again WITHOUT deleting anything in the target folder (only mismatching files will be redownloaded and run will continue from there)".format(len(missing_filelist), len(bad_filelist), ",".join(missing_filelist), ",".join(bad_filelist)))
		return ok_filelist, allisfine
		# end of nested subfunctions
		
	#TODO: check if a flag file exists, that indicates that this already ran successfully
	import glob
	assert targetfolder and os.path.abspath(targetfolder) != os.getcwd(), "\nERROR: You must provide a temporary targetfolder name!\n Temporary Target folder will be created, if it doesn't already exist.\n Temporary Target folder can NOT be the current working directory, as it will be deleted later!\n"
	temp_download_dict = {}
	#### first get the the crc-checksum file:
	sys.stderr.write("\n\tFirst downloading CRCchecksum file\n")
		### checking for (and deleting) older crc files (shit happens. first download-attempt may have happened before update, and then resumed after update, leading to two different CRCchecksum files)
	preexisting_crc = glob.glob(targetfolder + "/release*.files.installed") 
	if len(preexisting_crc) != 0:
		sys.stderr.write("\nWARNING: found the following pre-existing refseq-CRCchecksum-files in the targetfolder '{}' : {}.\n These are probably remnants of earlier download-attempts and will be deleted now\n".format(targetfolder, ", ".join(preexisting_crc)))
		for f in preexisting_crc:
			os.remove(f)
		### now downloading CRC
	_download_unixwget(sourcedict["crc"]["url"], sourcedict["crc"]["pattern"], targetdir=targetfolder, verbose=verbose) 
	downloaded_crc = glob.glob(targetfolder + "/release*.files.installed")
	assert len(downloaded_crc) == 1, "\nERROR: Expected exactly 1 'release*.files.installed' file in {} after download, but found {}! Be honest, are you messing with me?\n".format(targetfolder, len(downloaded_crc))
	crcfile = downloaded_crc[0]
	refseq_release_number = re.search("(release\d+)\.files\.installed", crcfile).group(1)
	#### Now that We have the crc file, download the rest:
	trycounter = 0
	download_list = None
	allisfine = False
	while not allisfine and trycounter < maxtries:
		if trycounter > 0:
			sys.stderr.write("\t-->Reattempting download\n")
		for vireukcat in [ vec for vec in sourcedict if vec != "crc" ]:
			sys.stderr.write("\n\tNow downloading refseq release category \"{}\" (attempt {})...\n".format(vireukcat, trycounter + 1))
			returncode = _download_unixwget(sourcedict[vireukcat]["url"], sourcedict[vireukcat]["pattern"], targetdir=targetfolder)
			if returncode != 0:
				sys.stderr.write("\nWARNING: wget returned non-zero returncode '{}' after downloading {} \n".format(returncode, vireukcat))
		download_list, allisfine = check_hashfile(crcfile, targetfolder, ",".join(["{}*.protein.faa.gz".format(x) for x in sourcedict if x != "crc"]))
		trycounter += 1
	#TODO: create a flag file to indicate that this has already run
	assert allisfine, "\nERROR: Still incomplete download or mismatching MD5sums after {} download-attempts. Please check connection and try again later...\n".format(trycounter+1)
	return download_list, refseq_release_number

def download_gtdb_stuff(sourcedict = gtdb_source_dict, targetfolder=None, verbose=False, maxtries=3):
	"""
	downloads all protein-fastas from the refseq release subfolders specified in source_dict.
	sourcedict should have categories as keys and a dictionary of corresponding remote folder urls and searchpatterns as values
	returns a dictionary listing the names of all files downloaded for each category/subfolder if successful, None if not.
	if md5sums do not match it will reattempt downloads as often as speciefied by 'maxtries' (default = 3 times)
	"""
	#TODO: check for a flag file, indicating that this already ran successfully
	# start of nested subfunctions
	def which_md5filename(targetdir):
		"""
		on the gtdb downloadserver md5files sometimes have '.txt' suffix and sometimes not. check which one it is this time...
		"""
		import glob
		return glob.glob(targetdir + "/" + MD5FILEPATTERN_GTDB)[0] # --> assumes there is only one hit, therefore takes only the first of the list returned by glob.glob(); todo: make sure md5sum file is always deleted after db-setup! otherwise there may be problems if preexisting dbs are updated
	
	def check_gtdbmd5file(md5_file, targetdir, patternstring): # TODO: merge with check_hashfile() to ONE function (not that hard)!
		"""
		as with _download_unixwget(), 'pattern' can be a comma seperated list of patterns (concatenated to a single string)
		raises AssertionException if something is wrong with the downloaded files, returns a list of downloaded filenames in everything is fine
		"""
		import re
		subpattern="_(r\d+)"
		import fnmatch
		patternlist = patternstring.split(",")
		infile = openfile(md5_file)
		allisfine = True
		ok_filelist, bad_filelist, missing_filelist = [], [], []
		for line in infile:
			tokens = line.strip().split()
			checksum = tokens[0]
			filename = os.path.basename(tokens[1])
			filename = re.sub(subpattern, "", filename) #gtdb renames all files in the latest folder (removing version from filename), but not in the corresponding md5-file... 
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
					ok_filelist.append(expectedfile)
		infile.close()
		ok_filelist.sort() #ensure alphabetial ordering of files, even if order is different in MD5-file
		assert len(ok_filelist) != 0, "\nERROR: None of the expected files appear to have been downloaded. Do you have read/write permissions for the targetfolder?\n"
		if not allisfine:
			sys.stderr.write("\nWARNING: something went wrong during download: {} expected files are missing, and {} files have mismatching MD5-checksums!\n  --> Missing files: {}\n  --> Corrupted files: {}\n".format(len(missing_filelist), len(bad_filelist), ",".join(missing_filelist), ",".join(bad_filelist)))
		return ok_filelist, allisfine

	def get_gtdbversion_from_MD5SUM(md5sumfile):
		import re
		infile = openfile(md5sumfile)
		versionpattern = "_(r\d+)\.tsv.gz" #just aiming on the compressed taxonomy tables which should always be there
		for line in infile:
			tokens=line.strip().split()
			if len(tokens) < 2:
				continue
			patternmatch = re.search(versionpattern, tokens[1])
			if patternmatch:
				version = patternmatch.group(1)
				return version
	# end of nested subfunctions
	
	def get_download_dict(sourcedict, targetfolder):
		download_dict = {}
		for x in sourcedict:
			okdownloadfilelist, allisfine = check_gtdbmd5file(which_md5filename(targetfolder), targetfolder, sourcedict[x]["pattern"])
			if not allisfine:
				return
			download_dict[x] = okdownloadfilelist
		return download_dict
			
	download_dict = None
	trycounter = 0
	while download_dict == None and trycounter < maxtries:
		if trycounter > 0:
			sys.stderr.write("\t-->Reattempting download\n")
		for gtdbcat in sourcedict:
			sys.stderr.write("\n\tNow downloading from gtdb: \"{}\" (attempt {})...\n".format(gtdbcat, trycounter+1))
			returncode = _download_unixwget(sourcedict[gtdbcat]["url"], sourcedict[gtdbcat]["pattern"], targetdir=targetfolder, verbose=verbose)
			if returncode != 0:
				sys.stderr.write("\nWARNING: wget returned non-zero returncode '{}' after downloading {} \n".format(returncode, gtdbcat))
		download_dict = get_download_dict(sourcedict, targetfolder)
		trycounter += 1
	assert download_dict, "\nERROR: Still incomplete download or mismatching MD5sums after {} download-attempts. Please check connection and try again later...\n".format(trycounter +1)
	version = get_gtdbversion_from_MD5SUM(which_md5filename(targetfolder))
	return download_dict, version

def download_silva_stuff(sourcedict = silva_source_dict, targetfolder=None, verbose=False, maxtries = 3):
	# start of nested subfunctions
	def check_silvamd5file(filelist, wishlist):
		allisfine = True
		successlist = []
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
		return successlist, allisfine
	
	def getsilvaversion(urldict, targetfolder): #yes I know. This could go easier with urllib2 yadayadayada. But for the other stuff wget is better, so sticking with that for now!
		sys.stderr.write("\n\tNow trying to get current silva release version...\n")
		versionfilename = os.path.join(targetfolder, urldict["wishlist"][0])
		#if os.path.exists(versionfilename):
		#	sys.stderr.write("\nDeleting pre-existing {}\n".format(versionfilename))
		#	os.remove(versionfilename)
		url = "url"
		returncode = _download_unixwget(urldict[url] + urldict["wishlist"][0], targetdir=targetfolder, verbose=verbose)
		if returncode != 0:
			sys.stderr.write("\thttps url not reachable. trying ftp-url\n")
			url = "alturl"
			returncode = _download_unixwget(urldict[url] + urldict["wishlist"][0], targetdir=targetfolder, verbose=verbose)
			if returncode != 0:
				sys.exit("\n\tERROR: can't reach silva database. please try again later\n")	
		with open(versionfilename) as versionfile:
			version = versionfile.read().strip()
		return version, url
	# end of nested subfunctions
	
	def get_download_dict(prelim_downloadlist, wishdict):
		download_dict = {}
		for x in wishdict:
			okdownloadlist, allisfine = check_silvamd5file([df for df in prelim_downloadlist if not df.endswith(".md5")], wishdict[x])
			if not allisfine:
				return
			download_dict[x] = okdownloadlist
		return download_dict
	
	version , url = getsilvaversion(sourcedict["silva_version"], targetfolder)
	prelim_downloadlist = [os.path.join(targetfolder, sourcedict["silva_version"]["wishlist"][0])]
	wishdict = {silvacat : [ w.format(version) for w in sourcedict[silvacat]["wishlist"] ] for silvacat in sourcedict }
	download_dict = None
	trycounter = 0
	while download_dict == None and trycounter < maxtries:
		if trycounter > 0:
			sys.stderr.write("\t-->Reattempting download\n")
		for silvacat in ["silva_taxfiles", "silva_fastas"]:
			for wish in wishdict[silvacat]:
				sys.stderr.write("\n\tNow downloading from silva: \"{}\" (attempt {})...\n".format(wish, trycounter +1))
				# ~ sys.stderr.write(sourcedict[silvacat][url] + wish + "\n")
				returncode = _download_unixwget(sourcedict[silvacat][url] + wish, pattern = None, targetdir=targetfolder, verbose=verbose)
				if returncode != 0:
					sys.stderr.write("\nWARNING: wget returned non-zero returncode '{}' after downloading {} \n".format(returncode, wish))
				prelim_downloadlist.append(os.path.join(targetfolder, wish))
		download_dict = get_download_dict(prelim_downloadlist, wishdict)
		trycounter += 1
	assert download_dict, "\nERROR: Still incomplete download or mismatching MD5sums after {} download-attempts. Please check connection and try again later...\n".format(trycounter +1)
	#TODO: create a flag file, indicating that this already ran successfully
	return download_dict, version
		
		
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
			   "eukcat__viral": {"parent": "r__Viruses", "rank" : -1, "taxname" : "Viruses"}, \
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

		acc2taxidfile.write("0\t{}\t{}\t0\n".format(acc, taxlist[-1])) #"dummy" fields with value 0 in columns 1 and 4, to make it look like the original acc2taxid files expected by "_create_sorted_acc2taxid_lookup()" in getdb.py (interesting infos in columns 2&3)
		linecount += 1
		if linecount % 500 == 0:
			sys.stderr.write("\r\tread {} lines and added {} taxa so far..".format(linecount, len(taxdict)))
			sys.stderr.flush()
	sys.stderr.write("\r\tread {} lines and added {} taxa so far..\n".format(linecount, len(taxdict)))
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
	sys.stderr.write("\n\treading silva data from {}\n".format(infilename))
	sys.stderr.flush()
	
	for line in infile:
		tokens =line.strip().split("\t")
		if len(tokens) != 6 or tokens[0] == "primaryAccession": #TODO: not safe to assume this validation check (expect 6 tab-seperated columns per line) will always be valid. this may change
			continue
		acc=".".join(tokens[:3]) # silva taxmap format: column1 = contigname, column2 = start, column3 = stop. --> accession = contigname.start.stop
		if tokens[3].startswith("Eukaryota"):
			taxtokens = [tokens[3].rstrip(";").split(";")[0]] #The eukaryotic taxonomy system is totally different than the bacterial/archaeal one. Since it is already hard enough to consolidate silva and gtdb taxonomies, for now (to make things easier) ignoring all specific classifications above domain for eukaryotes! (This tool is meant for bacteria/archaea, eukaryotes schould just be recognized as contamination)
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

def _concat_fastas(contig_fastalist, outfastahandle, return_headerdict = False, remove_prodigalIDs = False, remove_descriptions = False): 
	"""
	concatenates the fastafiles given in contig_fastalist
	if return_headers == True, thisreturns dictionary with assemblyIDs as keys and lists of corresponding fasta-sequenceIDs as values. otherwise returns None
	if remove_prodigalIDs == True, it will remove the protein id from the fastaids (based on the regex "_\d+ ")
	(--> naming the proteins simply after the contigs they originated from. additional protein info will remain in the sequence description)
	if remove_descriptions == True, it will remove all description info from the headers (keep only the sequenceID)
	"""
	import re
	# ~ print("\n_concat_fastas\n")
	# ~ print(outfastahandle.name + "\n")
	assemblyIDpattern = re.compile("GC._\d+\.\d")
	prodigalinpattern = re.compile("(^>\w+(\.\d+){0,1})(_\d+)(.*)")
	ignore_patterns = [ re.compile("\.tsv$") ] #stored as a list, to make it easier to add additional patterns if things change on the gtdb server
	prodigalreplacement = r"\1 \4"
	headerdict = {}
	filecounter = seqcounter = 0
	for cf in contig_fastalist:
		ignore = False
		if os.path.isdir(cf):
			continue	
		filebasename = os.path.basename(cf)
		for ip in ignore_patterns:
			if re.search(ip, filebasename):
				ignore = True
				break
		if ignore == False and re.search(assemblyIDpattern, filebasename) == None:
			sys.stderr.write("\nWARNING: encountered unexpected file '{}'! Probably something changed in the organization of the gtdb-download-server. Ignoring for now, BUT MAKE SURE THE WORKFLOW STILL WORKS!\n".format(filebasename))
			ignore = True
		if ignore:
			continue
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
			sys.stderr.write("\r\tread {} of {} files so far, wrote {} sequences to {}".format(filecounter, len(contig_fastalist), seqcounter, outfastahandle.name))
			sys.stderr.flush()
	sys.stderr.write("\r\tread {} of {} files so far, wrote {} sequences to {}\n".format(filecounter, len(contig_fastalist), seqcounter, outfastahandle.name))
	sys.stderr.flush()
	return headerdict if return_headerdict else None
				
def gtdb_contignames2taxids(contig_fastalist, acc2taxidinfilelist, acc2taxidoutfilename, outfastahandle):
	"""
	return new acc2taxid-filename
	Note:	sadly the refrence file at gtdb do not follow the same naming scheme for protein and genomic files (or their taxid-lookuptables). 
			the seqids in their taxid-lookuptables as well as the protein files have a "RS_" or "GB_" prefix, but the genomic files do not.
			have to remove that prefix here
	"""
	import re
	assemblyIDpattern = re.compile("GC._\d+\.\d")
	headerdict = _concat_fastas(contig_fastalist, outfastahandle, return_headerdict = True, remove_prodigalIDs = False, remove_descriptions = True)
	from mdmcleaner import getdb #todo: delete this. is only for debuggung
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
				sys.stderr.write("\r\tassigned {} contigs to {} taxids".format(contigcounter, linecounter))
				sys.stderr.flush()
		inaccfile.close()
	outaccfile.close()
	with open("delmetemp_repcheckfile.tsv", "wt") as repcheckfile: #todo: this is only for checking sanity of gtdb downloads (which genomes are considered "representative" and have seqdate provided? which are missing? are all species/genera/phyla covered? which are missing)? comment this out later
		repcheckfile.write("{}\n".format("\n".join(["{}\t{}".format(x[0], x[1]) for x in reptuple_list]))) #todo: this is only for checking sanity of gtdb downloads (which genomes are considered "representative" and have seqdate provided? which are missing? are all species/genera/phyla covered? which are missing)? comment this out later
	sys.stderr.write("\r\tassigned {} contigs to {} taxids\n".format(contigcounter, linecounter))
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
		os.remove(filename)
	acc2taxidfile.close()
	return acc2taxidoutfilename
		
		
		
		
def _download_dbdata_nonncbi(targetdir, progressdump, verbose=False): #todo: pass verbosity from calling function
	"""
	Downloads GTDB, SILVA and Refseq Eukaryote data, together with taxonomic info, to create a smaller and more condensed reference dataset compared to NCBI-nr
	The resulting taxonomy will be based almost exclusively on the GTDB taxonomy (except for Eukaryotes)
creates files for storing taxdict, LCA_walktree and acc2taxid files and returns the corresponding filenames [TODO: NOPE IT DOESNT]
	if 'continueflag' is set to True, it will assume an existing targetdir to be from an aborted previous run and will try to continue (or restart) from there
	if 'continueflag' is set to False and targetdir already exists, it will return an AssertionError
	returns gtdb_download_dict, refseq_files, silva_download_dict (in that order)
	"""
	sys.stderr.write("\nWARNING! This download may take a LONG time! However, you may be able to resume if it is aborted prematurely...\n")
	sys.stderr.flush()
	assert os.path.abspath(targetdir) != os.getcwd(), "\nERROR: targetdir may not be the current working directory\n"
	assert progressdump["step"] in _progress_steps["download"], "\nError: loaded progress step '{}' is not involved in Download steps {}!\n".format(progressdump["step"], str(_progress_steps["download"]))
	steporder = _progress_steps["download"]

	#step 1: download gtdb
	step = steporder[1]
	stepdesc = "download GTDB data"
	laststep = None 
	sys.stderr.write("\n--{}: {}--\n".format(step, stepdesc))
	currentprogressmarker = "progress_step{}.json".format(step)
	if progressdump["step"] == laststep:
		progressdump["step"] = step
		progressdump["gtdb_download_dict"], progressdump["gtdb_version"] = download_gtdb_stuff(gtdb_source_dict, targetdir, verbose=verbose)
		sys.stderr.write("\tGTDB_version = {}\n".format(progressdump["gtdb_version"]))
		getdb.dict2jsonfile(progressdump, os.path.join(targetdir, currentprogressmarker))
	else:
		sys.stderr.write("-->already did this step earlier. skipping it now! (delete progress marker '{}' and higher if you want to redo it)\n".format(currentprogressmarker))
		sys.stderr.flush()
		
	#step 2: download refseq-release eukaryotes
	laststep = step
	step = steporder[steporder.index(step) + 1]
	stepdesc = "downlad RefSeq-release Eukaryotes and Viruses"
	sys.stderr.write("\n--{}: {}--\n".format(step, stepdesc))
	lastprogressmarker = currentprogressmarker 
	currentprogressmarker = "progress_step{}.json".format(step)
	if progressdump["step"] == laststep:
		progressdump["step"] = step
		progressdump["refseq_files"], progressdump["refseq_release_number"] = download_refseq_eukaryote_prots(refseq_dbsource_dict, targetdir, verbose=verbose)
		sys.stderr.write("\tRefSeq release = {}\n".format(progressdump["refseq_release_number"]))
		getdb.dict2jsonfile(progressdump, os.path.join(targetdir, currentprogressmarker))
		os.remove(os.path.join(targetdir, lastprogressmarker))
	else:
		sys.stderr.write("-->already did this step earlier. skipping it now! (delete progress marker '{}' and higher if you want to redo it)\n".format(currentprogressmarker))
		sys.stderr.flush()
		
	#step 3: download silva stuff
	laststep = step
	step = steporder[steporder.index(step) + 1]
	stepdesc = "downlad SILVA data"
	sys.stderr.write("\n--{}: {}--\n".format(step, stepdesc))
	lastprogressmarker = currentprogressmarker 
	currentprogressmarker = "progress_step{}.json".format(step)
	if progressdump["step"] == laststep:
		progressdump["step"] = step
		progressdump["silva_download_dict"], progressdump["silva_version"] = download_silva_stuff(silva_source_dict, targetdir, verbose=verbose)
		sys.stderr.write("\tSilva_version = {}\n".format(progressdump["silva_version"]))
		getdb.dict2jsonfile(progressdump, os.path.join(targetdir, currentprogressmarker))
		os.remove(os.path.join(targetdir, lastprogressmarker))
	else:
		sys.stderr.write("-->already did this step earlier. skipping it now! (delete progress marker '{}' and higher if you want to redo it)\n".format(currentprogressmarker))
		sys.stderr.flush()

	with openfile(os.path.join(targetdir, "DB_versions.txt"), "wt") as versionsfile:
		versionsfile.write("GTDB version = {}\n".format(progressdump["gtdb_version"]))
		versionsfile.write("RefSeq release = {}\n".format(progressdump["refseq_release_number"]))
		versionsfile.write("silva_download_dict = {}\n".format(progressdump["silva_version"]))
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
		sys.stderr.write("\ntargetdir '{}' does not exist yet --> creating it!\n".format(targetdir)) #todo: add exception for insufficient write permissions
		return {"step" : None}
	markerfile = markertag = None
	markerpattern="progress_step*.json*"
	#markertagsubpattern = re.compile("step(\d+[a-z]*)\.")
	markerfilelist = sorted([ m for m in os.listdir(targetdir) if fnmatch.fnmatch(m, os.path.join(markerpattern)) ])
	if len(markerfilelist) > 0:
		if len(markerfilelist) > 1:
			sys.stderr.write("\nwarning: found multiple progress-markers: {}. Expected only one, so will be using only the latest one: '{}'. if this is incorrect, please delete progressmarkers and run again\n".format(str(markerfilelist), markerfilelist[-1]))
		markerfile = markerfilelist[-1] #if for unexpected reasons multiple progressmarkers exist, accept only the latest one
		#markertag = re.search(markertagsubpattern, os.path.basename(markerfile)).group(1)
	sys.stdout.flush()
	return getdb.jsonfile2dict(os.path.join(targetdir, markerfile)) if markerfile else {"step" : None}
	
	
	
def _prepare_dbdata_nonncbi(targetdir, progressdump, verbose=False, settings=None): #Todo:add continueflag argument #todo: make sure verbosity is implemented
	from mdmcleaner import blasthandler
	from mdmcleaner import getdb #todo: for using dict2jsonthis will be obsolete, when this is moved there #edit: no it won't. keeping download for gttdb and ncbi data seperate. common stuff goes to misc or getdb
	import time #todo probably not needed anymore

	if settings == None:
		settings = default_settings #by default assume all tools in pATH

	def step4a():
		progressdump["acc2taxidfilelist"] = []
		progressdump["taxdict"], progressdump["lca_walktree"] = None, None
		for gtdbtax in gtdb_download_dict["gtdb_taxfiles"]:
			progressdump["taxdict"], progressdump["lca_walktree"], acc2taxidfile = read_gtdb_taxonomy_from_tsv(gtdbtax, progressdump["taxdict"], progressdump["lca_walktree"]) #todo: change read_gtdb_taxonomy_from_tsv() so that it accepts a filehandle to write acc2taxid-lookup to
			progressdump["acc2taxidfilelist"].append(acc2taxidfile)
		
	def step4b():#todo take progressdump as input, return progressdump as output
		#TODO: create seperate concat_fastas (and seperate diamond-dbs) for gtdb and refseq (maybe even for each refseq-vireukat category) 
		assert len(gtdb_download_dict["gtdb_fastas"]) == 2, "\nERROR: number of downloaded gtdb tars different than expected: should be 2 but is {}\n".format(len(gtdb_download_dict)) #Notifying myself, so that i don't forget to adapt this if i ever choose to download more or less reference categories from gtdb
		sys.stderr.write("\t-->unpacking genome files : {}\n".format(gtdb_download_dict["gtdb_fastas"][0]))
		gtdb_genomefiles = misc.untar(gtdb_download_dict["gtdb_fastas"][0], targetdir, removetar=True) 
		sys.stderr.write("\t-->unpacking protein files : {}\n".format(gtdb_download_dict["gtdb_fastas"][1]))
		gtdb_proteinfiles = misc.untar(gtdb_download_dict["gtdb_fastas"][1], targetdir, removetar=True)
		progressdump["gtdb_genomefiles"] = gtdb_genomefiles
		progressdump["gtdb_proteinfiles"] = gtdb_proteinfiles
		return gtdb_genomefiles, gtdb_proteinfiles

	def step4c():
		# ~ import pdb; pdb.set_trace()	
		progressdump["concatprotfasta"] = os.path.join(targetdir, "concat_refprot.faa")
		progressdump["concatgenomefasta"] = os.path.join(targetdir, "concat_refgenomes.fasta")
		concatgenomefastahandle = openfile(progressdump["concatgenomefasta"], "wt")
		concatprotfastahandle = openfile(progressdump["concatprotfasta"], "wt")
		sys.stderr.write("\t-->concatenating genomes\n")
		progressdump["acc2taxidfilelist"].append(gtdb_contignames2taxids(gtdb_genomefiles, progressdump["acc2taxidfilelist"], os.path.join(targetdir, "temp_acc2taxid_gtdb.acc2taxid.gz"), concatgenomefastahandle))
		concatgenomefastahandle.close()
		sys.stderr.write("\t-->concatenating proteins\n")
		_concat_fastas(gtdb_proteinfiles, concatprotfastahandle, return_headerdict = False, remove_prodigalIDs = True, remove_descriptions = False) #does not return an acc2taxidfile ,because that is already covered by the contig-accessions
		return concatgenomefastahandle, concatprotfastahandle 
	
	def step5a(): #todo take progressdump as input, return progressdump as output
		for f in silva_download_dict["silva_taxfiles"]:
			sys.stderr.write("\treading {}\n".format(f)) #todo: debugging message
			sys.stderr.flush()
			taxdict, LCA_walktree, filename = read_silva_taxonomy_from_tsv(f, progressdump["taxdict"], progressdump["lca_walktree"])
			progressdump["acc2taxidfilelist"].append(filename)
		return taxdict, LCA_walktree, progressdump["acc2taxidfilelist"] #todo: instead of returng these variables, just put them into pgrogressdump and return THAT


	def step6a(concatprotfastahandle):
		progressdump["acc2taxidfilelist"].append(refseq_contignames2taxids(refseq_files, concatprotfastahandle, os.path.join(targetdir, "temp_acc2taxid_refseq.acc2taxid.gz")))


	def step7a():
		progressdump["finalacc2taxidfile"] = getdb._create_sorted_acc2taxid_lookup(progressdump["acc2taxidfilelist"], os.path.join(targetdir, acc2taxid_outfilebasename))
		progressdump["acc2taxidfilelist"] = []

				
	#####end of nested subfunctions 
	#Todo: Make this more elegant if there is any time for that, later
	#step 4: get gtdb_taxonomy
	#TODO: make seperate diamond databases for GTDB-bacteria, GTDB-Archaea and Refseq-eukaryotes. Maybe even a sperate one for each eukaryote-category! Test if parallel blasts against these are faster than blasting against a single large db 
	steporder = _progress_steps["prepare"]
	laststep = progressdump["step"]
	assert laststep <= max(steporder), "\nError: last loaded progress step '{}' comes AFTER the preparation steps {}!\n".format(laststep, str(steporder))
	gtdb_download_dict, refseq_files, silva_download_dict = progressdump["gtdb_download_dict"], progressdump["refseq_files"], progressdump["silva_download_dict"]
	lastprogressmarker = "progress_step{}.json".format(laststep)
	step = steporder[0]
	stepdesc = "reading gtdb taxonomy"
	sys.stderr.write("\n--{}: {}--\n".format(step, stepdesc))
	currentprogressmarker = "progress_step{}.json".format(step)
	#step 04a
	if laststep == _progress_steps["download"][-1]:
		progressdump["step"] = step
		step4a() ##step 4a: read gtdb_taxonomy files (which use the assembly/genome_accessions as identifiers). NOTE: preliminary refseq taxonomy (knows only domain-level for Eujaryota and Viruses) autmatically added to taxodnomy_dict at this point
		getdb.dict2jsonfile(progressdump, os.path.join(targetdir, currentprogressmarker))
		os.remove(os.path.join(targetdir, lastprogressmarker))
	else:
		sys.stderr.write("-->already did this step earlier. skipping it now! (delete progress marker '{}' and higher if you want to redo it)\n".format(currentprogressmarker))
		sys.stderr.flush()
		
	lastprogressmarker, laststep = currentprogressmarker, step
	step = steporder[steporder.index(step) + 1]
	stepdesc = "unpacking gtdb sequence files"
	sys.stderr.write("\n--{}: {}--\n".format(step, stepdesc))
	currentprogressmarker = "progress_step{}.json".format(step)
 	#step 04b (unpack gtdb tars und store extracted filenames
	if progressdump["step"] == laststep:
		progressdump["step"] = step
		gtdb_genomefiles, gtdb_proteinfiles = step4b() ##step 4b: uncompress gtdb-tar-files
		getdb.dict2jsonfile(progressdump, os.path.join(targetdir, currentprogressmarker))
		os.remove(os.path.join(targetdir, lastprogressmarker))
	else:
		#concatgenomefastahandle =  openfile(progressdump["concatgenomefasta"], "at")
		gtdb_genomefiles = progressdump["gtdb_genomefiles"]
		gtdb_proteinfiles = progressdump["gtdb_proteinfiles"]
		sys.stderr.write("-->already did this step earlier. skipping it now! (delete progress marker '{}' and higher if you want to redo it)\n".format(currentprogressmarker))
		sys.stderr.flush()
		
	lastprogressmarker, laststep = currentprogressmarker, step
	step = steporder[steporder.index(step) + 1]
	stepdesc = "read sequence identifiers and add to taxonomy lookupfile"
	sys.stderr.write("\n--{}: {}--\n".format(step, stepdesc))
	currentprogressmarker = "progress_step{}.json".format(step)
 	##tep 4c:read contig names for each uncompressed genome-file --> assign contigs to taxids and write to acc2taxidfile
	if progressdump["step"] == laststep:
		progressdump["step"] = step
		concatgenomefastahandle, concatprotfastahandle = step4c() 	
		getdb.dict2jsonfile(progressdump, os.path.join(targetdir, currentprogressmarker))
		os.remove(os.path.join(targetdir, lastprogressmarker))
	else:
		#concatgenomefastahandle =  openfile(progressdump["concatgenomefasta"], "at")
		concatprotfastahandle = openfile(progressdump["concatprotfasta"], "at")
		sys.stderr.write("-->already did this step earlier. skipping it now! (delete progress marker '{}' and higher if you want to redo it)\n".format(currentprogressmarker))
		sys.stderr.flush()

		
	#step5: get silva taxonomy
	lastprogressmarker, laststep = currentprogressmarker, step 
	step = steporder[steporder.index(step) + 1]
	stepdesc = "integrating Silva taxonomy and preparing taxonomy lookup files"
	sys.stderr.write("\n--{} : {}--\n".format(step, stepdesc))
	#step 05a
	if progressdump["step"] == laststep:
		currentprogressmarker = "progress_step{}.json".format(step)
		progressdump["step"] = step
		sys.stderr.write("\t-->adding silva...\n")
		progressdump["taxdict"], progressdump["lca_walktree"], progressdump["acc2taxidfilelist"] = step5a() ##step 5a: read silva_taxonomy files and add this info to taxdict and LCA_walktree. taxdict and lca_walktree should be complete at this point 
		sys.stderr.write("\t-->writing taxdict file...\n")
		sys.stdout.flush()		
		taxdict_file = getdb.dict2jsonfile(progressdump["taxdict"], os.path.join(targetdir, taxdb_outfilebasename)) #this is the final db db file for this
		progressdump["taxdict_file"] = taxdict_file
		sys.stderr.write("\t-->writing walktree file...\n")
		sys.stdout.flush()	
		progressdump["lca_walktree_file"] = getdb.dict2jsonfile(progressdump["lca_walktree"], os.path.join(targetdir, "lca_walktreedict.json.gz")) #for saving progress. still needs to be turned to actual lcawalkdb
		sys.stderr.write("\t-->building walktree...\n")
		sys.stdout.flush()
		progressdump["lcawalkdb_file"] = getdb.build_lca_db(progressdump["lca_walktree"], os.path.join(targetdir, lcawalkdb_outfilebasename), "root")
		progressdump["lca_walktree"] = progressdump["taxdict"] = "" #saving memory
		getdb.dict2jsonfile(progressdump, os.path.join(targetdir, currentprogressmarker))
		os.remove(os.path.join(targetdir, lastprogressmarker))
	else:
		#concatgenomefastahandle =  openfile(progressdump["concatgenomefasta"], "at")
		concatprotfastahandle = openfile(progressdump["concatprotfasta"], "at")
		sys.stderr.write("-->already did this step earlier. skipping it now! (delete progress marker '{}' and higher if you want to redo it)\n".format(currentprogressmarker))
		sys.stderr.flush()
			
	#step 6: cocatenate all reference datasets, create blastdbs, and delete associated fastafiles
	lastprogressmarker, laststep = currentprogressmarker, step 
	step = steporder[steporder.index(step) + 1]
	stepdesc = "concatenate all refseq references and add to taxonomy lookup files"
	sys.stderr.write("\n--{}: {}--\n".format(step, stepdesc))
	currentprogressmarker = "progress_step{}.json".format(step)
	#step 06a
	if progressdump["step"] == laststep:
		progressdump["step"] = step
		step6a(concatprotfastahandle) ##step6a: concatenate all refseq protein reference-fasta files (to outhandle concatfasta), keep outhandle (concatfasta) and return acc2taxidfile(name)
		concatprotfastahandle.close()
		getdb.dict2jsonfile(progressdump, os.path.join(targetdir, currentprogressmarker))
		os.remove(os.path.join(targetdir, lastprogressmarker))
	else:
		#concatgenomefastahandle =  openfile(progressdump["concatgenomefasta"], "at")
		sys.stderr.write("-->already did this step earlier. skipping it now! (delete progress marker '{}' and higher if you want to redo it)\n".format(currentprogressmarker))
		sys.stderr.flush()
	
	lastprogressmarker, laststep = currentprogressmarker, step
	step = steporder[steporder.index(step) + 1]
	stepdesc = "create protein diamond database"
	sys.stderr.write("\n--{}: {}--\n".format(step, stepdesc))
	currentprogressmarker = "progress_step{}.json".format(step)
	#step 06b
	if progressdump["step"] == laststep:
		progressdump["step"] = step
		progressdump["outprotdbname"] = os.path.join(targetdir, "gtdbplus_protdb")
		blasthandler.make_diamond_db(progressdump["concatprotfasta"], progressdump["outprotdbname"], diamond = settings["diamond"])  ##step6b: create protein diamond-db
		getdb.dict2jsonfile(progressdump, os.path.join(targetdir, currentprogressmarker))	
		os.remove(os.path.join(targetdir, lastprogressmarker))
	else:
		#concatgenomefastahandle =  openfile(progressdump["concatgenomefasta"], "at") #todo: make concatenated gtdb genomes-fasta optional (just using 16S/18S data takes up much less space)
		sys.stderr.write("-->already did this step earlier. skipping it now! (delete progress marker '{}' and higher if you want to redo it)\n".format(currentprogressmarker))
		sys.stderr.flush()	
			
	lastprogressmarker, laststep = currentprogressmarker, step
	step = steporder[steporder.index(step) + 1]
	stepdesc = "create nucleotide BLAST databases"
	sys.stderr.write("\n--{}: {}--\n".format(step, stepdesc))
	currentprogressmarker = "progress_step{}.json".format(step)
	#step 06c #todo: make a step06c() subfunction for this
	if progressdump["step"] == laststep:
		progressdump["step"] = step
		dbfasta_list = [progressdump["concatgenomefasta"]] + silva_download_dict["silva_fastas"] #todo: create a new checkpoint for each created database... until then: UNCOMMENT THIS LINE!
		#dbfasta_list = silva_download_dict["silva_fastas"] #todo: for debugging only. choose line above for full run or implement seperate checkpoints for each...
		progressdump["outnucldblist"] = []
		for dbfasta in dbfasta_list:
			#todo : make the following code less dumb and ugly:
			outnucldbname = dbfasta
			for suff in [".gz", ".fasta", ".fa", ".faa", ".fna"]:
				if outnucldbname.endswith(suff):
					outnucldbname = outnucldbname[:-len(suff)]
			progressdump["outnucldblist"].append(blasthandler.make_blast_db_from_gz(dbfasta, outnucldbname, db_type = "nucl",makeblastdb = settings["makeblastdb"])) ##step6d: create nucleotide diamond or blastdb (test which one is faster, blast may be faster here, because there are not so many queries as in the protein-blasts)
		getdb.dict2jsonfile(progressdump, os.path.join(targetdir, currentprogressmarker))
		os.remove(os.path.join(targetdir, lastprogressmarker)) 
		#TODO: make sure tehre are actually 3 databases created here! one for 16S , one for 23S and one for gtdbgenomes!
	else:
		sys.stderr.write("-->already did this step earlier. skipping it now! (delete progress marker '{}' and higher if you want to redo it)\n".format(currentprogressmarker))
		sys.stderr.flush()
	#TODO: step 7: create combined sorted taxid_lookupfile
	#step 7a
	# ~ print("HERE SHOULD COME STEP 7a!")
	lastprogressmarker, laststep = currentprogressmarker, step
	step = steporder[steporder.index(step) + 1]
	stepdesc = "combine and sort taxid2lookup-files"
	sys.stderr.write("\n--{}: {}--\n".format(step, stepdesc))
	currentprogressmarker = "progress_step{}.json".format(step)
	if progressdump["step"] == laststep:
		progressdump["step"] = step
		step7a()
		getdb.dict2jsonfile(progressdump, os.path.join(targetdir, currentprogressmarker))
		os.remove(os.path.join(targetdir, lastprogressmarker)) 
	else:
		sys.stderr.write("-->already did this step earlier. skipping it now! (delete progress marker '{}' and higher if you want to redo it)\n".format(currentprogressmarker))
		sys.stderr.flush()		
	sys.stderr.write("\n{}\nfinished downloading and processing reference database!\n".format("="*50))
	sys.stderr.flush()
	return progressdump

def test_or_create_targetdir(targetdir):
	try:
		os.makedirs(targetdir, exist_ok=True)
		if not os.access(targetdir, os.W_OK):
			raise PermissionError
	except FileExistsError as e:
		raise FileExistsError("\n{}\n\nERROR: File '{}' should be a directory but is a file! Please choose a target directory and try again\n".format(e, targetdir))
	except PermissionError as e:
		raise PermissionError("\n{}\n\nERROR: insufficient write permissions for '{}'! Please choose a different target directory and try again\n".format(e, targetdir))


def getNprepare_dbdata_nonncbi(targetdir, verbose=False, settings=None):
	if settings == None:
		settings = default_settings #by default assume all tools in pATH
	#TODO: pack all data into pgrogessdump
	#  Todo: pass only progressdump to subsequent functions
	#  Todo: have a dictionary of progresssteps and their respecitve functions.	
	sys.stderr.write("Specified target directory '{}'\n".format(targetdir))
	test_or_create_targetdir(targetdir)
	progressdump = _check_progressmarker(targetdir)
	#steps1-3 (downloading):
	sys.stderr.write("highest pre-existing progressmarker found: '{}'\n".format(progressdump["step"]))
	sys.stderr.flush()
	if progressdump["step"] in _progress_steps["download"][:-1]:
		sys.stderr.write("beginning/continuing at download stage\n")
		progressdump = _download_dbdata_nonncbi(targetdir, progressdump, verbose=verbose)
	#steps 4-6 (processing)
	if progressdump["step"] in [_progress_steps["download"][-1]] + _progress_steps["prepare"][:-1]:
		sys.stderr.write("beginning/continuing at processing stage\n")
		progressdump = _prepare_dbdata_nonncbi(targetdir, progressdump, verbose=verbose, settings = settings)
	if progressdump["step"] == _progress_steps["prepare"][-1]:
		#cleanup
		cleanupwhenfinished(progressdump, targetdir, verbose)
	elif progressdump["step"] == _progress_steps["finished"][-1]:
		sys.stderr.write("\nfully processed database already exists at '{}'. If this should be replaced, please delete and try again, or specify a different target directory\n".format(targetdir))

def cleanupwhenfinished(progressdump, targetdir, verbose=False):
	'''
	deletes remaining intermediary files and logs component database versions
	'''
	import re
	from mdmcleaner import getdb
	sys.stderr.write("\n--CLEANING UP--\n")
	def deletefiles(stufflist):
		restlist = []
		for stuff in stufflist:
			# ~ print(stuff)
			# ~ sys.stdout.flush()		
			if os.path.exists(stuff) and os.path.isfile(stuff):
				os.remove(stuff)
				if verbose:
					sys.stderr.write("\tdeleted {}\n".format(stuff))
			else:
				restlist.append(stuff)
				

	def deletedirs(stufflist):
		for stuff in stufflist:
			# ~ print(stuff)
			# ~ sys.stdout.flush()
			if os.path.exists(stuff) and os.path.isdir(stuff):
				try:
					os.rmdir(stuff)
					sys.stderr.write("\tdeleted {}\n".format(stuff)) #todo: only when verbose
				except OSError:
					sys.stderr.write("\nWARNING: could not delete intermediary directory {}!. Either not empty or permissions not sufficient? Please check and delete manually, if necessary!\n".format(stuff))

	def NestedDictValues(stuffdict):
		for stuff in stuffdict.values():
			if isinstance(stuff, dict):
				yield from NestedDictValues(stuff)
			else:
				yield stuff
	# ~ import pdb; pdb.set_trace()
	laststep = progressdump["step"]
	lastprogressmarker = "progress_step{}.json".format(laststep)
	progressdump["step"] = _progress_steps["finished"][-1]
	currentprogressmarker = "progress_step{}.json".format(progressdump["step"])
	delcategories = ["silva_download_dict", "gtdb_download_dict", "gtdb_genomefiles", "gtdb_proteinfiles", "refseq_files", "concatprotfasta", "concatgenomefasta"]
	for dc in delcategories:
		# ~ print("---{}".format(dc))
		if isinstance(progressdump[dc], dict):
			dellists = [ x for x in NestedDictValues(progressdump[dc])]
			flattened_dellist = [y for x in dellists for y in x] + [y + ".md5" for x in dellists for y in x]
		elif isinstance(progressdump[dc], list):
			flattened_dellist = progressdump[dc]
		else:
			flattened_dellist = [progressdump[dc]]
		deletefiles(flattened_dellist)
		deletedirs(sorted(flattened_dellist, key=lambda x:len(x), reverse=True))
		progressdump[dc] = ""
		
			 
	delversionfiles = ["MD5SUM", "release209.files.installed", "VERSION.txt"]
	delfastafiles = ["concat_refgenomes.fasta", "concat_refprot.faa"]
	deletefiles([os.path.join(targetdir, x) for x in delversionfiles + delfastafiles])

	getdb.dict2jsonfile(progressdump, os.path.join(targetdir, currentprogressmarker))
	os.remove(os.path.join(targetdir,lastprogressmarker))
#todo: DROP download of ncbi2gtdb and vice versa mapping files will not use them anyway!

def get_publication_set(args, configs):
	from mdmcleaner import misc
	sys.stderr.write("\nDownloading publicaton reference-dataset from Zenodo (Warning: this is definitively NOT the most recent reference dataset!)\n")
	sys.stderr.flush()
	if args.outdir == None:
		args.outdir = "./db"
	if os.path.exists(os.path.join(args.outdir,"gtdb")):
		sys.exit("\ntargetfolder already exists: '{}'. will not overwrite! Please delete it if you want to use this destination, or choose another target-folder!\n")
	tarfile = os.path.join(args.outdir, os.path.basename(zenodo_publication_db))
	_download_unixwget(zenodo_publication_db, pattern=None, targetdir=args.outdir, verbose=False)
	tar_hash = misc.calculate_md5hash(tarfile)
	if zenodo_publication_db_md5 != tar_hash:
		sys.exit("\nERROR during download: md5hash of downloaded tar ({}) does not match hash of zenodo link ({})\n".format(tar_hash, zenodo_publication_db_md5))
	sys.stderr.write("unacking tar file")
	sys.stderr.flush()
	tar_contents = misc.untar(tarfile, targetdir=args.outdir, removetar = True)
	checksumdict = misc.check_md5file(os.path.join(args.outdir, "gtdb", "md5sum.txt"))
	if False in [ checksumdict[x]["ok"] for x in checksumdict ]:
		sys.stderr.write("\nERROR: following files had mismatching md5-hashes: {}\n--> PLEASE DELETE AND TRY AGAIN!\n".format("\n\t-".join([ checksumdict[x] for x in checksumdict if checksumdict[x]["ok"] == False])))
	sys.stderr.write("\npublication-reference dataset successfully downloaded! To use it as a reference dataset for mdmcleaner runs, specify it in a (local) mdmcleaner.config file via the following command:\n"
	                  "mdmcleaner.py set_configs --db_basedir {}\n".format(os.path.abspath(args.outdir)))
	
	
	
	

def main(args, configs): #todo: make option to read targetdir from configfile or to WRITE targetdir to configfile when finished
	if args.get_pub_data:
		get_publication_set(args, configs)
	else:
		getNprepare_dbdata_nonncbi(args.outdir, verbose=args.verbose, settings=configs.settings)

if __name__ == '__main__':
	sys.stderr.write("testing gtdb_download\n")
	download_gtdb_stuff(sourcedict = gtdb_source_dict, targetfolder=".", verbose=True, maxtries=3)
	# ~ import argparse
	# ~ myparser = argparse.ArgumentParser(prog=os.path.basename(sys.argv[0]), description = "Download and prepare GTDB- and SILVA-data for MDMcleaner")
	# ~ myparser.add_argument("-o", "--outdir", action = "store", dest = "outdir", default = "./db/gtdb", help = "target directory for gtdb-/silva-data. may not be the current working directory. Default = './db/gtdb'")
	# ~ myparser.add_argument("--verbose", action = "store_true", dest = "verbose", default = False, help = "verbose output (download progress etc)") #todo: finish implementing
	# ~ myparser.add_argument("--quiet", action = "store_true", dest = "quiet", default = False, help = "quiet mode (suppress any status messages except Errors and Warnings)") #todo: implement
	# ~ args = myparser.parse_args()
	# ~ main(args)

