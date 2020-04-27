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
import time
import traceback

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

# hardcoded for now: where does NCBI currently store the stuff? If that changes, edit the ftp-adresses here:
# also, binrefiner only considers proteins for now, but including other accession-dbs also for compatability with KRONA
ncbi_ftp_server = "ftp.ncbi.nlm.nih.gov"
#ftp_adress_taxonomy = "https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/" #https (urlib + urlib2 version)
ftp_adress_taxonomy = "pub/taxonomy/"
ftp_adress_accessiondbs = "{}accession2taxid/".format(ftp_adress_taxonomy)
ftp_source_taxdmp = "{}taxdump.tar.gz".format(ftp_adress_taxonomy)
ftp_source_taxdmpnew = "{}/new_taxdump/new_taxdump.tar.gz".format(ftp_adress_taxonomy) #ncbi is trying out a new format for taxdump-files. Might as well prepare for those also already...
ftp_source_prot_acc2taxid = "{}prot.accession2taxid.gz".format(ftp_adress_accessiondbs)
ftp_source_prot_acc2taxiddead= "{}dead_prot.accession2taxid.gz".format(ftp_adress_accessiondbs)
ftp_source_nucl_acc2taxid= "{}nucl.accession2taxid.gz".format(ftp_adress_accessiondbs)
ftp_source_nucl_acc2taxiddead= "{}dead_nucl.accession2taxid.gz".format(ftp_adress_accessiondbs)

testfile= "pub/taxonomy/accession2taxid/dead_prot.accession2taxid.gz"
_dbsource_dict = { 	"taxdmp" : ftp_source_taxdmp, \
					"taxdmp_new" : ftp_source_taxdmpnew, \
					"prot_acc2taxid" : ftp_source_prot_acc2taxid, \
					"prot_acc2taxid_dead" : ftp_source_prot_acc2taxiddead, \
					"nucl_acc2taxid" : ftp_source_nucl_acc2taxid, \
					"nucl_acc2taxid_dead" : ftp_source_nucl_acc2taxiddead }

rank2index = { "no rank" : 0, \
				"ignored rank" : 1, \
				"superkingdom" : 10, \
				"phylum" : 20, \
				"class" : 30, \
				"order" : 40, \
				"genus" : 50, \
				"species" : 60 } # using increments of 10 in case i want to use the indermediate ranks (e.g. subfamily) at some later point also

index2rank = { rank2index[key] : key for key in rank2index }


taxdb_outfilebasename = "taxonomy_br.json" #using the same as krona but exporting it into a faster to read format
acc2taxid_outfilebasename = "all.accession2taxid.sorted" #using the same as krona


def openfile(infilename, filemode = "rt"): #proably move this to a "basics" module since other steps need it too?
	""" a convenience function that will be moved to another more general module later"""
	if infilename.endswith(".gz"):
		import gzip
		filehandle = gzip.open(infilename, filemode)
	else:
		filehandle = open(infilename, filemode)
	return filehandle 

def  _download_db(dbtype, targetdir="."):
	""" 
	downloads from a specific set of ncbi taxonomy database files.
	dbtype should be one of: 
		- "taxdmp" (for the nodes.dmp and names.dmp files)
		- "taxdmp_new" (a new experimental format at ncbi. Might as well prepare for that already
		- "prot_acc2taxid" (for assigning protein accession numbers to taxids)
		- "prot_acc2taxid_dead" (for assigning old, now deleted protein accession numbers to taxids)
	the ftp download paths used are given in "getdb.ftp_source_<dbtype>", respectively.
	"""
	import urllib.request #TODO: switch to urllib2 for better error/timeout-handling?
	import tarfile
	import time
	
	##### some encapsulated subfunctions: ##################################
	def reporthook(blockcount, blocksize, totalsize): #todo make this a class with a updatable totalcounter shared between instances
		if blockcount % 100 == 0:
			totaldownload = (blockcount * blocksize)
			gb = 1000**3
			mb = 1000**2
			currentstatgb = float(totaldownload)/(gb)
			totalsizegb = float(totalsize)/(gb)
			perc = totaldownload / totalsize
			sys.stderr.write("\r-->Downloaded {:.3f} GB of {:.3f} GB ({:.1%})".format(currentstatgb, totalsizegb, perc))
		
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
		raise KeyError("\nERROR: unknown dbtype '{}'!\nMust be one of {}\n".format(dbtype, ", ".join(sorted(_dbsource_dict.keys()))))
	tempfilelist = [ os.path.join(targetdir, "delmetemp_{}".format(os.path.basename(_dbsource_dict[dbtype]))) ]

	#attempt to download, and verify md5sums
	attempt_counter = 1
	while True:
		if not (os.path.exists(tempfilelist[0]) and os.path.exists(tempfilelist[0] + ".md5")):
			sys.stderr.write("\ndownloading {} (attempt Nr. {})\n".format(_dbsource_dict[dbtype], attempt_counter))
			start = time.time()
			urllib.request.urlretrieve(_dbsource_dict[dbtype], os.path.join(targetdir, tempfilelist[0]), reporthook)
			end = time.time()
			sys.stderr.write("\ndownload took {:.1f} hours\n".format((end-start)/3600))
			urllib.request.urlretrieve(_dbsource_dict[dbtype] + ".md5", os.path.join(targetdir, tempfilelist[0] + ".md5"))
			
		else:
			sys.stderr.write("\n{} was apparently already downloaded...\n".format(os.path.basename(_dbsource_dict[dbtype])))
		if checkmd5(tempfilelist[0], tempfilelist[0] +  ".md5"):
			break
		sys.stderr.write("--> but md5sums don't match! (need to try again)\n")
		sys.stderr.flush()
		#remove faulty files:
		os.remove(tempfilelist[0])
		os.remove(tempfilelist[0] + ".md5")
		assert attempt_counter < max_attempts, "exceeded number of tries to download database. Please check connection!\n"
		attempt_counter += 1
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

def _download_db2(dbtype, targetdir="."): #test ftp-module --> works. TODO: rename to download_db
	import ftplib
	#TODO: ensure that there is a check for already downloaded files before redownloading them
	gb = 1000**3
	max_attempts = 10
	
	if dbtype not in _dbsource_dict: #usage help
		raise KeyError("\nERROR: unknown dbtype '{}'!\nMust be one of {}\n".format(dbtype, ", ".join(sorted(_dbsource_dict.keys()))))
	
	ftpfilelist = (_dbsource_dict[dbtype], _dbsource_dict[dbtype] + ".md5")
	outfilenamelist = tuple((os.path.join(targetdir, os.path.basename(f)) for f in ftpfilelist )) #can't add prefix "delmetemp_" because that changes the md5sum
	#outfile = openfile(outfilename, "wb")
	#outmd5file = openfile(outfilename + ".md5", "w")

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
	
	def mycallback(data):
		nonlocal total_blocks_downloaded
		nonlocal total_bytes_downloaded
		total_blocks_downloaded += 1
		currentblocksize = len(data) #hope that works like this with binary data...
		total_bytes_downloaded += currentblocksize
		total_bytes_downloaded_gb = float(total_bytes_downloaded) / gb
		perc = total_bytes_downloaded_gb / final_size_gb
		if total_blocks_downloaded % 100 == 0:
			sys.stderr.write("\r-->Downloaded {:.3f} GB of {:.3f} GB ({:.1%})".format(total_bytes_downloaded_gb, final_size_gb, perc))
		outfile.write(data)

	current_attempts = 1
	success = False
	sys.stderr.write("connecting to {}\n".format(ncbi_ftp_server))
	while current_attempts < max_attempts and success == False:
		
		#outfilelist = [ openfile(f, "wb") for f in outfilenamelist ] #moved file creation to the dowload part, so that it can be skipped if file already exists.
		#TODO: add a check wether file not oly already exists, but is actually not older than the ftp version?
		total_bytes_downloaded = 0
		total_blocks_downloaded = 0
		try:
			sys.stderr.flush()
			sys.stderr.write("\tattempt no. {}\n".format(current_attempts))
			ftp_conn = ftplib.FTP(ncbi_ftp_server, user="anonymous")
			ftp_conn.sendcmd("Type i") # set connection to binary
		except Exception as e:
			sys.stderr.write("Error during connection attempt:")
			sys.stderr.write("\n{}\n{}\n".format(e, traceback.print_exc()))
			current_attempts += 1
			sys.stderr.write("\nwill try again in two seconds\n")
			ftp_conn.close()
			time.sleep(2) #apparently i get a name resolution error if i try to reconnect immediately
			continue
		try:
			for ftpfile, outfilename in zip(ftpfilelist, outfilenamelist):
				if not os.path.exists(outfilename):
					outfile = openfile(outfilename, "wb")
					sys.stdout.write("\n{} + {} \n".format(ncbi_ftp_server, ftpfile))
					final_size_gb = ftp_conn.size(ftpfile) / gb
					ftp_conn.retrbinary("RETR " + ftpfile, mycallback)
					outfile.close()
				else:
					sys.stderr.write("outfile already exists: {} --> skipping for now\n".format(outfilename))
		except EOFError as e:
			sys.stderr.write("\n{}\n".format(e, traceback.print_exc()))
			current_attempts += 1
			sys.stderr.write("\nERROR: Got an 'EOFError'. Connection must have suddenly been broken:\n")
			ftp_conn.close()
			outfile.close()
			time.sleep(2)
			continue			
		except Exception as e: #TOdo: find a way to pick download up where it aborted last
			sys.stderr.write("\n{}\n".format(e, traceback.print_exc()))
			current_attempts += 1
			sys.stderr.write("Unexpected Error during download attempt:")
			sys.stderr.write("\nwill try again in 2 seonds\n")
			ftp_conn.close()
			outfile.close()
			time.sleep(2)
			continue
		ftp_conn.quit()
		if checkmd5(outfilenamelist[0], outfilenamelist[1]):
			success = True
			break
		else:
			sys.stderr.write("\MD5checksums do not match! will try to download again\n")
			assert False, "no i wont" #something is going wrong. have to find out what
			for f in outfilenamelist:
				os.remove(f) #delete erroneous outfilescat
			current_attempts += 1
			time.sleep(2)
			continue
	assert success, "\nFTP download of {} failed!\n".format(os.path.join(ncbi_ftp_server, ftpfile))
	ftp_conn.quit()
	
	return outfilenamelist #just handle them all downstream (delete "md5-file" when done. first is always the db, second is only the md5-file




def _download_db3(dbtype, targetdir="."): #fuck it! urllib2 and ftlib don't work for me, so who cares about another dependancy? This at least is easy to use and works
	import wget
	if dbtype not in _dbsource_dict: #usage help
		raise KeyError("\nERROR: unknown dbtype '{}'!\nMust be one of {}\n".format(dbtype, ", ".join(sorted(_dbsource_dict.keys()))))
	
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
		
	max_attempts = 10
	ftpfilelist = ("ftp://" + ncbi_ftp_server + "/" + _dbsource_dict[dbtype], "ftp://" + ncbi_ftp_server + "/" +  _dbsource_dict[dbtype] + ".md5")
	outfilenamelist = tuple((os.path.join(targetdir, "delmetemp_" + os.path.basename(f)) for f in ftpfilelist ))
	current_attempts = 0
	
	while current_attempts < max_attempts:
		success = False
		current_attempts += 1
		for ftpfile, outfilename in zip(ftpfilelist, outfilenamelist):
			if not os.path.exists(outfilename):
				sys.stderr.write("\nnow downloading '{}':\n".format(ftpfile))
				sys.stderr.flush()
				try:
					wget.download(ftpfile, outfilename)
					success = True
				except Exception:
					success = False
					break				
			else:
				sys.stderr.write("\noutfile already exists: {} --> skipping for now\n".format(outfilename))
				success = True
				sys.stderr.flush()
		if success and checkmd5(outfilenamelist[0], outfilenamelist[1]):
			sys.stderr.write("md5sums match! Download successful\n")
			break
		elif success:
			success = False
			sys.stderr.write("\nmd5sums do not match! deleting files and trying again in 2 seconds...\n")
			for i in outfilenamelist:
				os.remove(i)
		else:
			sys.stderr.write("\nsomething went wrong during download! Will try again in 2 seconds...\n")
		time.sleep(2) #wait 2 seconds before next attempt
	assert success, "\nFAILED: Could not download {}\n".format(ftpfilelist[0])
	return outfilenamelist
			
			

def json_taxdb_from_dmp(download_dir = "."):
	"""
	Creates a taxonomy lookup file from downloaded or provided "nodes.dmp" and "names.dmp" files.
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
			rank =  rank2index.get(tokens[3], 1) # ranks other than the 7 major ranks will considered 1 = "ignored rank"
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
		
	def taxdict2json(taxdict, targetdir): #assume that targetdir = downloaddir
		import json
		outdbfilename = od.path.join(targetdir, acc2taxid_outfilebasename)
		with open(outdbfilename, 'wb') as outfile:
			json.dump(taxdict, outfile)
		return outdbfilename
		
	##### end of encapsulated subfunctions ################################
	taxdict = read_nodesdmp(os.path.join(download_dir, "nodes.dmp"))
	taxdict = read_namesdmp(os.path.join(download_dir, "names.dmp"), taxdict)
	infile.close()
	return taxdict2json(taxdict, download_dir) #is smaller than kronas db. herfore do not attempt to build kronas db

def json_taxdb_from_kronadb(kronadb_dir = "."):
	ktdb_basename = "taxonomy.tab"
	kronadb = os.path.join(kronadb_dir, ktdb_basename)
	infile = openfile(kronadb)
	for line in infile:
		pass #finish this sometime
		

def download_accessiondb(targetdir=".", dbmode = "minimal", enforce=False): #TODO: implement filecheck and ask for user-feedback (skippable with "-f" or "-y" argument) if it looks as if a krona-db file could be overwritten
	acc2taxid_outfilename = os.path.join(targetdir, acc2taxid_outfilebasename) #change to whatever krona is using
	if os.path.isfile(acc2taxid_outfilename):
		warning = "ATTENTION: if you are planning to use/update a database also intended for use with Krona, it is recommended to use Kronas \"updateAccessions.sh\" script instead!\n"
		if enforce != True:
			raise IOError("\nERROR: {} already exists and \"-f\" argument not set. will not overwrite!\n".format(acc2taxid_outfilename) + warning)
		else:
			sys.stderr.write("\nARNING: {} already exists. --> will overwrite it!\n".format(acc2taxid_outfilename) + warning)
	downloaded_filelist = []
	if dbmode == "minimal":
		wishlist = ["prot_acc2taxid", "prot_acc2taxid_dead"]
	elif dbmode == "krona":
		wishlist = [ "prot_acc2taxid", "prot_acc2taxid_dead", "nucl_acc2taxid" , "nucl_acc2taxidn_dead" ]# if you want to update the db in a way that it is usable by KRONAtools also, is also needs the nucleotide accessions...
	for w in wishlist:
		downloaded_filelist.append(_download_db2(w, targetdir))
	sys.stderr.write("finished downloading accessiondbs for selection: '{}'".format(dbmode))
	#extract and sort them:
	_create_sorted_acc2taxid_lookup(downloaded_filelist, acc2taxid_outfilename)
	return acc2taxid_outfilename
	

def _create_sorted_acc2taxid_lookup(acc2taxidfilelist, acc2taxid_outfilename):
	'''
	for creating an alphabetically sorted list of acessions and their corrsponding taxids
	for use with binary search
	TODO: try making this presort multiple input files in paralall (if multiple cpus are specified)
	'''
	import subprocess
	#import shlex #allow string splitting accoring to Shell -like syntax
	presortcmd = "zcat {infile} | cut -f 2,3| grep -v accession | sed 's#\.[0-9]*##'| sort > {outfile}"
	finalsortcmd = "sort -m {filelist} > {finaldb}"
	tempfilelist = []
	for f in acc2taxidfilelist:
		if f.endswith(".gz"):
			tempfile = f[:-3]
		tempfile = tempfile + ".sorted"
		sys.stderr.write("extracting and presorting {}\n".format(f))
		sout, serr = subprocess.Popen(["bash", "-c", presortcmd.format(infile=f, outfile=tempfile)], stdout=subprocess.PIPE).communicate()
		if serr != None:
			raise RuntimeError("...extraction exited with Error:\n{}\n".format(serr))
		os.remove(f)
		os.remove(f + ".md5")
		tempfilelist.append(tempfile)
	sys.stderr.write("combining sorted tempfiles\n")
	sout,serr = subprocess.Popen(["bash", "-c", finalsortcmd.format(filelist=" ".join(tempfilelist), finaldb=acc2taxid_outfilename)], stdout=subprocess.PIPE).communicate()
	if serr != None:
			raise RuntimeError("...extraction exited with Error:\n{}\n".format(serr))
	#now clean up:
	sys.stderr.write("removing temporary downloads")
	for f in tempfilelist:
		os.remove(f)

class taxdb(object):
	def __init__(self, acc2taxid_lookupfile, taxdbfile):
		#self.taxdict = self.read_taxddbfile(taxdbfile) #todo: write tis!
		self.acc_lookup_handle = openfile(acc2taxid_lookupfile) #keep in mind that this function will be moved to another module
		self.acc_lookup_handle_filesize = self.acc_lookup_handle.seek(0,2) #jump to end of file and give bytesize (alternative to "os.path.getsize()" which hopefully should also work with compressed files)
		self.taxdict = None
		#reminer to self: do NOT use compressed acc2taxid_lookupfile. It increases time for binary search ca 200x!

	def read_taxdb(self, taxdbfile):#needs to be json format. Krona taxdbs need to be converted to this format first, using the kronadb2json function above
		import json
		infile = openfile(taxdbfile)
		self.taxdict = json.load(infile)
		infile.close()
	
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
		for queryacc in sorted(queryacclist): #VERY important that queries are sorted as well!
			outdict[queryacc], start = self.acc2taxid(queryacc, start)
		return outdict 
		
def _test_modulewide_globals(): # to reassure me that global variables within modules really function as intended. Left in here to remind me of that until this whole thing is completed
	print("hi can you see this")
	print(rank2index)
	print("hope you do...")

def _test_download():
	sys.stderr.write("_test_dowload\n")
	sys.stderr.flush()
	download_accessiondb(".", "minimal", False)
	sys.stderr.write("\nfinished\n")
	sys.stderr.flush()

def _test_download2():
	sys.stderr.write("_test_dowload2\n")
	sys.stderr.flush()
	_download_db3("taxdmp")
	
	
def _test_makeaccdb():
	sys.stderr.write("\ntest_makeaccdb\n")
	sys.stderr.flush()
	acc2taxidfilelist = [	"delmetemp_prot.accession2taxid.gz" , \
							"delmetemp_dead_prot.accession2taxid.gz" ]
	acc2taxid_outfilename = "sorted_shit_yeah.db"
	_create_sorted_acc2taxid_lookup(acc2taxidfilelist, acc2taxid_outfilename)
	sys.stderr.write("\nFINISHED\n")
if __name__ == '__main__':	
	_test_download2()
