#!/usr/bin/env python
""" 
downloading, parsing and modifying taxonomy databases for binrefiner.py.

Uses the following files of the ncbi taxonomy database
	- "nodes.dmp" (from "taxdump.tar.gz" or "new_taxdump.tar.gz")
	- "names.dmp" (from "taxdump.tar.gz" or "new_taxdump.tar.gz")
	- "prot.accession2taxid" (and, optionally, "dead_prot.accession2taxid.gz")
To decrease memory consumption, taxonomix ranks are converted to integer-indices.
These indices can be mapped back and forth using the dictionaries "getdb.rank2index" and/or "getdb.index2rank"
"""

import os
import sys



'''
module for downloading and parsing the ncbi taxonomy databases for bin refinerwork in progress
things to consider:
	-maybe use SQLite for databases (check if that uses less memory)
	-maybe just use already existing python SQLite implementation for ncbi TAxonomy: taxadb?
		--> if that is not "overkill" for this rather streamlined purpose
	-or instead try to streamline own dictionary by minimizing use of strings as values/keys 
		- turn accession numbers into ASCII or UTF-8 or similar encoded character strings (should use less RAM?)
		- replace taxlevel strings with integers representing each level (10= Domain, 20=Phylum, 30=Class, 50=Order, 60=Family, 70=Genus, 80=Species, ignore subclasses, subfamilies etc and strains for now)
	---> implement these different options and compare;
		A.) RAM usage
		B.) size of the dictionary pickle/compressed JSON file on he hard disk
'''

# TODO: change all stderr messages to logger statements!!! make verbosity adjustable!!!

# hardcoded for now: where does NCBI currently store the stuff? May change:
# also, only considering proteins for now...
ftp_source_taxdmp = "https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz"
ftp_source_taxdmpnew= "https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/new_taxdump/new_taxdump.tar.gz" #ncbi is trying out a new format for taxdump-files. Might as well prepare for those also already...
ftp_source_acc2taxid= "ftp://ftp.ncbi.nlm.nih.gov/pub/taxonomy/accession2taxid/prot.accession2taxid.gz"
ftp_source_acc2taxiddead= "ftp://ftp.ncbi.nlm.nih.gov/pub/taxonomy/accession2taxid/dead_prot.accession2taxid.gz"
_dbsource_dict = { 	"taxdmp" : ftp_source_taxdmp, \
					"taxdmp_new" : ftp_source_taxdmpnew, \
					"acc2taxid" : ftp_source_acc2taxid, \
					"acc2taxid_dead" : ftp_source_acc2taxiddead }

rank2index = { "no rank" : 0, \
				"ignored rank" : 1, \
				"superkingdom" : 10, \
				"phylum" : 20, \
				"class" : 30, \
				"order" : 40, \
				"genus" : 50, \
				"species" : 60 } # using increments of 10 in case i want to use the indermediate ranks (e.g. subfamily) at some later point also

index2rank = { rank2index[key] : key for key in rank2index }

def openfile(infilename, filemode = "rt"): #proably move this to a "basics" module since other steps need it too?
	""" a convenience function that will be moved to another more general module later"""
	if infilename.endswith(".gz"):
		import gzip
		filehandle = gzip.open(infilename, filemode)
	else:
		filehandle = open(infilename, filemode)
	return filehandle 

def download_db(dbtype, targetdir="."):
	""" 
	downloads from a specific set of ncbi taxonomy database files.
	dbtype should be one of: 
		- "taxdmp" (for the nodes.dmp and names.dmp files)
		- "taxdmp_new" (a new experimental format at ncbi. Might as well prepare for that already
		- "acc2taxid" (for assigning protein accession numbers to taxids)
		- "acc2taxid_dead" (for assigning old, now deleted protein accession numbers to taxids)
	the ftp download paths used are given in "getdb.ftp_source_<dbtype>", respectively.
	"""
	import urllib.request
	import tarfile
	
	##### some encapsulated subfunctions: ##################################
	def reporthook(blockcount, blocksize, totalsize): #todo make this a class with a updatable totalcounter shared between instances
		perc = (blockcount * blocksize) / totalsize 
		sys.stderr.write("\r-->Downloaded {:.0%}".format(perc))
		
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
	##### end of encapsulated subfunctions #################################
	
	sys.stderr.flush()
	
	max_attempts = 3 #maximum number of times to try to re-download the database if md5sums don't match
	
	if dbtype not in _dbsource_dict: #usage help
		raise KeyError("\nERROR: unknown dbtype '{}'!\nMust be one of {}\n".format(dbtype, ", ".join(_dbsource_dict.keys())))
	tempfilelist = [ os.path.join(targetdir, "delmetemp_{}".format(os.path.basename(_dbsource_dict[dbtype]))) ]

	#attempt to download, and verify md5sums
	attempt_counter = 0
	while True:
		attempt_counter += 1
		sys.stderr.write("\ndownloading {} (attempt Nr. {})\n".format(tempfilelist[0], attempt_counter))
		urllib.request.urlretrieve(_dbsource_dict[dbtype], os.path.join(targetdir, tempfilelist[0]), reporthook)
		urllib.request.urlretrieve(_dbsource_dict[dbtype] + ".md5", os.path.join(targetdir, tempfilelist[0] + ".md5"))
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


#TODO: actually read these downloaded temfiles and create a database

def read_taxdict_from_dmp(download_dir = "."):
	"""
	Creates a taxonomy lookup dictionary from downloaded or provided "nodes.dmp" and "names.dmp" files.
	The files "nodes.dmp" and "names.dmp" should be located in "<download_dir>"
	Returns a dictionary for assigning taxids to scientific names and for tracing the full taxonomic lineage down to the rank "superkingdom".
	""" 
	#TODO: This should NOT return a dictionary. Instead it should create a Json or yaml or pickle or other suitable file for creating a taxdb object later
	#should just return the filename forthis 
	##### some encapsulated subfunctions: ##################################
	def read_nodesdmp(nodefile):
		taxdict =  {}
		infile = openfile(nodefile)
		for line in infile:
			tokens = line.rstrip().split("\t|\t") # convoluted delimiter, but that's what ncbi taxdump uses...
			taxid = int(tokens[0])
			parent = int(tokens[1])
			rank =  rank2index.get(tokens[3], 1) # ranks other than the major ranks will considered 1 = "ignored rank"
			taxname = None
			taxdict[taxid] = { "parent" : parent, \
								"rank" : rank, \
								"taxname" : taxname }
		infile.close()
		return taxdict
		
	def read_namesdmp(namesfile, taxdict):
		infile = openfile(namesfile)
		for line in infile:
			tokens = line.rstrip().rstrip("\t|").split("\t|\t") # convoluted delimiter, but that's what ncbi taxdump uses...
			taxid = int(tokens[0])
			taxname = tokens[1]
			taxname_class = tokens[3]
			if taxname_class != "scientific name":
				#print("---{}---".format(taxname_class))
				continue
			#print("{} -> yes".format(taxid))
			#print("--> {}".format(taxname))
			assert taxdict[taxid]["taxname"] == None, "{} There are different scientific names for the same taxid. Need to look into this and choose which one to take...".format(taxid)
			taxdict[taxid]["taxname"] = taxname
		return taxdict
	##### end of encapsulated subfunctions ################################
	taxdict = read_nodesdmp(os.path.join(download_dir, "nodes.dmp"))
	taxdict = read_namesdmp(os.path.join(download_dir, "names.dmp"), taxdict)
	infile.close()
	return taxdict

def read_accdict_from_file(acc2taxidfile): 
	### FORGT THIS! TOO LARGE FOR DICT
	### INSTEAD MAKE A FUNCTION THAT CREATES AN SIMPLIFIED TABLE FILE SORTED ALPHABETICALLY BY ACC
	### THEN SEARCH THAT USING BINARY SEARCH
	### ALTERNATIVELY TRY OUT SQLite??? (if faster)

#def create_sorted_acc2taxid_lookup(acc2taxidfile)
	'''
	for creating an alphabetically sorted list of acessions and their corrsponding taxids
	for use with binary search
	'''
#	pass

#def create_sqlite_acc2taxid_db(acc2taxidfile):
	'''
	for trying out the alternative SQLite method
	'''
#	pass
	
class taxdb(object):
	def __init__(self, acc2taxid_lookupfile, taxdbfile):
		#self.taxdict = self.read_taxddbfile(taxdbfile) #todo: write tis!
		self.acc_lookup_handle = openfile(acc2taxid_lookupfile) #keep in mind that this function will be moved to another module
		self.acc_lookup_handle_filesize = self.acc_lookup_handle.seek(0,2) #jump to end of file and give bytesize (alternative to "os.path.getsize()" which hopefully should also work with compressed files)
		#reminer to self: do NOT use compressed acc2taxid_lookupfile. It increases time for binary search ca 200x!

	def acc2taxid(self, queryacc,start = 0):
		#for using binary search on a simple (maybe compressed?) sorted textfile
		start = 0
		stop = self.acc_lookup_handle_filesize
		subjectacc = None
		
		while start < stop: 
			currentpos = int((start + stop) / 2)
			self.acc_lookup_handle.seek(currentpos, 0)

			if currentpos != start:
				self.acc_lookup_handle.readline()
			tokens =  self.acc_lookup_handle.readline().strip().split("\t")
			subjectacc = tokens[0]
			subjecttaxid = tokens[1]
			if subjectacc == queryacc:
				return subjecttaxid, start
			if subjectacc > queryacc:
				stop = currentpos
			else:
				start = currentpos
		return none, 0
	
	def acclist2taxiddict(self, queryacclist): #may still reconsider whether i actually want to output a dictionary or not
		start = 0
		outdict = {}
		for queryacc in sorted(queryacclist): #VERY important that queries are sored as well!
			outdict[queryacc], start = self.acc2taxid(queryacc, start)
		return outdict 
		
def _test_modulewide_globals(): # to reassure me that global variables within modules really function as intended. Left in here to remind me of that until this whole thing is completed
	print("hi can you see this")
	print(rank2index)
	print("hope you do...")

def _test_download():
	sys.stderr.write("hi")
	sys.stderr.flush()
	db="taxdmp"
	download_db(db)
	sys.stderr.write("\nfinished\n")
	sys.stderr.flush()

if __name__ == '__main__':	
	_test_download()
