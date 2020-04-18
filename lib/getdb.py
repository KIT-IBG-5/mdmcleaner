#!/usr/bin/env python
import os
import sys



'''
module for downloading and parsing the ncbi taxonomy databases for bin refiner
work in progress
things to consider:
	-maybe use SQLite for databases (check if that uses less memory)
	-maybe just use already existing python SQLite implementation for ncbi TAxonomy: taxadb?
		--> if that is not "overkill" for this rather streamlined purpose
	-or instead try to streamline own dictionary by minimizing use of strings as values/keys 
		- turn accession numbers into ASCII or UTF-8 or similar encoded character strings (should use less RAM?)
		- replace taxlevel strings with integers representing each level (1= Domain, 2=Phylum, 3=Class, 5=Order, 6=Family, 7=Genus, 8=Species, ignore subclasses, subfamilies etc and strains for now)
	---> implement these different options and compare;
		A.) RAM usage
		B.) size of the dictionary pickle file on he hard disk
'''

# TODO: change all stderr messages to logger statements!!! make verbosity adjustable!!!

# hardcoded for now: where does NCBI currently store the stuff? May change:
# also, only considering proteins for now...
taxdmp_source = "https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz"
taxdmpnew_source = "https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/new_taxdump/new_taxdump.tar.gz" #ncbi is trying out a new format for taxdump-files. Might as well prepare for those also already...
acc2taxid_source = "ftp://ftp.ncbi.nlm.nih.gov/pub/taxonomy/accession2taxid/prot.accession2taxid.gz"
acc2taxiddead_source = "ftp://ftp.ncbi.nlm.nih.gov/pub/taxonomy/accession2taxid/dead_nucl.accession2taxid.gz"
dbsource_dict = { 	"taxdmp" : taxdmp_source, \
					"taxdmp_new" : taxdmpnew_source, \
					"acc2taxid" : acc2taxid_source, \
					"acc2taxid_dead" : acc2taxiddead_source }



def openfile(infilename, filemode = "rt"): #proably move this to a "basics" module since other steps need it too?
	if infilename.endswith(".gz"):
		import gzip
		filehandle = gzip.open(infilename, filemode)
	else:
		filehandle = open(infilename, filemode)
	return filehandle 

def checkmd5(infile, md5file):
	import hashlib #for comparing md5-checksums of downloaded databases
	
	#calculate hash of downloaded file:
	blocksize = 2**20 #chunks to read file in
	with open(infile, "rb") as f:
		filehash = hashlib.md5()
		while True:
			data = f.read(blocksize)
			if not data:
				break
			filehash.update(data)
	md5ist = filehash.hexdigest()
	
	#read hash from md5-checkfile:
	with open(md5file, "r") as m:
		md5soll = m.readline().split()[0]
	
	sys.stderr.write("\n comparing md5checksums: ")
	sys.stderr.flush()
	
	return md5ist == md5soll

def download_db(dbtype, targetdir="."):
	import urllib.request
	import tarfile
	max_attempts = 3 #maximum number of times to try to re-download the database if md5sums don't match
	
	if dbtype not in dbsource_dict: #usage help
		raise KeyError("\nERROR: unknown dbtype '{}'!\nMust be one of {}\n".format(dbtype, ", ".join(dbsource_dict.keys())))
	tempfilelist = [ os.path.join(targetdir, "delmetemp_{}".format(os.path.basename(dbsource_dict[dbtype]))) ]

	#attempt to download, and verify md5sums
	attempt_counter = 0
	while True:
		attempt_counter += 1
		sys.stderr.write("\ndownloading {} (attempt Nr. {})\n".format(tempfilelist[0], attempt_counter))
		urllib.request.urlretrieve(dbsource_dict[dbtype], os.path.join(targetdir, tempfilelist[0]), reporthook)
		urllib.request.urlretrieve(dbsource_dict[dbtype] + ".md5", os.path.join(targetdir, tempfilelist[0] + ".md5"))
		if checkmd5(tempfilelist[0], tempfilelist[0] +  ".md5"):
			break
		sys.stderr.write("--> md5sums don't match! (need to try again)\n")
		sys.stderr.flush()
		assert attempt_counter <= 3, "exceeded number of tries to download database. Please check connection!\n"
	sys.stderr.write(" --> md5sums match! Successful!\n")
	sys.stderr.flush()
	
	#after file is successfully downloaded: IF it is a taxdump.tar.gz, we need some specific files from that tarball
	if dbtype.startswith("taxdmp"):
		tf = tarfile.open(tempfilelist[0], "r:gz")
		wantedfiles = ["nodes.dmp", "names.dmp"]
		for wantedfile in wantedfiles:
			tf.extract(wantedfile, targetdir)
		tf.close()
		for df in [ tempfilelist[0], tempfilelist[0] + ".md5" ]:
			os.remove(df)
		tempfilelist = [ os.path.join("targetdir", wf) for wf in wantedfiles ]
	return tempfilelist

def reporthook(blockcount, blocksize, totalsize): #todo make this a class with a updatable totalcounter shared between instances
	perc = (blockcount * blocksize) / totalsize 
	sys.stderr.write("\r-->Downloaded {:.0%}".format(perc))
	
	sys.stderr.flush()
#TODO: actually read these downloaded temfiles and create a database

def read_nodesdmp(infilehandle):
	delim="\t|\t" # convoluted delimiter, but that's what ncbi taxdump uses...
	taxdict =  {}
	infile = openfile(infilehandle)
	for line in infile:
		tokens = line.rstrip().split("\t|\t")
		tax_id = int(tokens[0])
		parent = int(tokens[1])
		

def test_download():
	sys.stderr.write("hi")
	sys.stderr.flush()
	db="taxdmp"
	download_db(db)
	sys.stderr.write("\nfinished\n")
	sys.stderr.flush()
	
test_download()
