#!/usr/bin/env python
""" 
downloading, parsing and modifying taxonomy databases for binrefiner.py.

Uses the following files of the ncbi taxonomy database
	- "nodes.dmp" (from "taxdump.tar.gz" or "new_taxdump.tar.gz")
	- "names.dmp" (from "taxdump.tar.gz" or "new_taxdump.tar.gz")
	- "prot.accession2taxid" (and, optionally, "dead_prot.accession2taxid.gz")
To decrease memory consumption, taxonomix ranks are converted to integer-indices (although sonly minimal amout of RAM saved this way).
Taxids were als supposed to be represented as intergers, but this was incompatible with json, slowed down database construction and saved only little RAM --> dropped
These indices can be mapped back and forth using the dictionaries "getdb.rank2index" and/or "getdb.index2rank"
"""

import os
import sys
import time
import traceback

from misc import openfile
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
				"ignored rank" : -1, \
				"superkingdom" : 10, \
				"phylum" : 20, \
				"class" : 30, \
				"order" : 40, \
				"family" : 50, \
				"genus" : 60, \
				"species" : 70 } # using increments of 10 in case i want to use the indermediate ranks (e.g. subfamily) at some later point also
				#rank "root" does not exist (as previously plannes. Instead checking for taxid=1 (= root)

index2rank = { rank2index[key] : key for key in rank2index }


taxdb_outfilebasename = "taxonomy_br.json.gz" #different from krona. additional wasteful "children" info per node (may get rid of it again after perfecting LCA). saved as compressed file because uncompression does not increase loading time much but decreases file size a LOT
acc2taxid_outfilebasename = "all.accession2taxid.sorted" #using the same as krona
lcawalkdb_outfilebasename = "lcawalkdb_br.db"

def  _download_db(dbtype, targetdir="."):
	""" 
	downloads from a specific set of ncbi taxonomy database files.
	dbtype should be one of: 
		- "taxdmp" (for the nodes.dmp and names.dmp files)
		- "taxdmp_new" (a new experimental format at ncbi. Might as well prepare for that already
		- "prot_acc2taxid" (for assigning protein accession numbers to taxids)
		- "prot_acc2taxid_dead" (for assigning old, now deleted protein accession numbers to taxids)
		- "nucl_acc2taxid" (optinal)
		- "nucl_acc2taxid_dead"
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
			
			

def lca_and_json_taxdb_from_dmp(download_dir = "."):
	"""
	Creates a taxonomy lookup file from downloaded or provided "nodes.dmp" and "names.dmp" files.
	The files "nodes.dmp" and "names.dmp" should be located in "<download_dir>"
	Returns a dictionary for assigning taxids to scientific names and for tracing the full taxonomic lineage down to the rank "superkingdom".
	""" 
	#TODO: This should NOT return a dictionary. Instead it should create a Json or yaml or pickle or other suitable file for creating a taxdb object later
	#should just return the filename forthis 
	##### some encapsulated subfunctions: ##################################
	def read_nodesdmp(nodefile): #test:keep track of parent AND children per node at the same time. Wasteful for RAM but easier for LCA...
		taxdict =  {}
		lca_walk_tree = {} #instead of wastefully adding the children info for each node directly to the taxonomy dict that are actually only needed to create the paht_lits for LCA later on, just make a second temporary dict for that 
		infile = openfile(nodefile)
		for line in infile:
			tokens = line.rstrip().split("\t|\t") # convoluted delimiter, but that's what ncbi taxdump uses...
			taxid = tokens[0]
			parent = tokens[1]
			#assert tokens[2] in rank2index, "\nNO OFFICIAL RANK: -{}- for taxid -{}-. maybe you meant this: -{}-".format(tokens[2], taxid, tokens[3])
			rank =  rank2index.get(tokens[2], -1) # ranks other than the 7 major ranks will considered 1 = "ignored rank"			
			taxname = None
			taxdict[taxid] = { "parent" : parent, \
								"rank" : rank, \
								"taxname" : taxname }
			if taxid not in lca_walk_tree: #necessary to make sure terminal leaves also have entry in lca_walk_tree
				lca_walk_tree[taxid] = {	"level" : None, \
											"children" : [] } 
			
			if parent in lca_walk_tree: #adding children info now also
				if parent == "1" and taxid == "1": #no need to link root to itself (would lead to "maximum recursion error during eulers walk)
					continue
				lca_walk_tree[parent]["children"].append(taxid)
			elif not (parent == "1" and taxid == "1"): #no need to link root to itself (would lead to "maximum recursion error during eulers walk)
				lca_walk_tree[parent] = { 	"level" : None, \
											"children" : [taxid] }
		infile.close()
		return taxdict, lca_walk_tree
		
	def read_namesdmp(namesfile, taxdict):
		infile = openfile(namesfile)
		for line in infile:
			tokens = line.rstrip().rstrip("\t|").split("\t|\t") # convoluted delimiter, but that's what ncbi taxdump uses...
			taxid = tokens[0]
			taxname = tokens[1]
			taxname_class = tokens[3]
			if taxname_class != "scientific name":
				#print("---{}---".format(taxname_class))
				continue
			#print("{} -> yes".format(taxid))
			#print("--> {}".format(taxname))
			#print(type(taxdict))
			#print(len(taxdict))
			#print(taxdict[1])
			assert taxdict[taxid]["taxname"] == None, "{} There are different scientific names for the same taxid. Need to look into this and choose which one to take...".format(taxid)
			taxdict[taxid]["taxname"] = taxname
		infile.close()
		return taxdict

	def get_level(taxdict, taxid, level = 0): #get level (=depth) of taxon within ncbi taxonomy-tree. Required for LCA later on...
		level += 1
		if taxid == "1":
			return level
		elif taxdict[taxid]["parent"] in taxdict:
			return get_level(taxdict, taxdict[taxid]["parent"], level)
		else:
			return None
	
	def add_levelinfo(taxdict, lca_walk_tree):
		for taxid in lca_walk_tree:
			lca_walk_tree[taxid]["level"] = get_level(taxdict, taxid)
		return lca_walk_tree
			
	def build_lca_db(lca_walk_tree, targetdir):
		"""
		iterates through lca_walk_tree from root to each branch in the form of a "eulers walk". 
		stores the visited nodes of this eulers walk in a list. stores corresponding node-depths in a second list
		thoses lists form the actual lookup-table for lca queries
		"""
		def walk(lca_walk_tree, walk_list, depth_list, currnode = "1"): #simple attempt to get eulers walk across taxdict
			walk_list.append(currnode)
			# ~ print("=========")
			# ~ print(currnode)
			# ~ print(lca_walk_tree[currnode])
			depth_list.append(lca_walk_tree[currnode]["level"])
			for child in sorted(lca_walk_tree[currnode]["children"]):
				#print("{} --> {}".format(currnode, child))
				walk_list, depth_list = walk(lca_walk_tree, walk_list, depth_list, child)
				walk_list.append(currnode)
				depth_list.append(lca_walk_tree[currnode]["level"])
			return walk_list, depth_list
		
		walk_list = []
		depth_list = []
		walk_list, depth_list = walk(lca_walk_tree, walk_list, depth_list)
		
		outfilename = os.path.join(targetdir, lcawalkdb_outfilebasename)
		outfile = openfile(outfilename, "wt")
		outfile.write("\t".join(walk_list) + "\n")
		outfile.write("\t".join([str(l) for l in depth_list]))
		outfile.close()
		return outfilename

	def taxdict2json(taxdict, targetdir): #assume that targetdir = downloaddir
		import json
		outdbfilename = os.path.join(targetdir, taxdb_outfilebasename)
		outfile = openfile(outdbfilename, 'wt')
		json.dump(taxdict, outfile)
		outfile.close()
		return outdbfilename

	# ~ def taxdict2yaml(taxdict, targetdir): #assume that targetdir = downloaddir #optional in cas i decide to switch to yaml
		# ~ import yaml
		# ~ outdbfilename = os.path.join(targetdir, taxdb_outfilebasename + ".yaml")
		# ~ with open(outdbfilename, 'w') as outfile:
			# ~ yaml.dump(taxdict, outfile)
		# ~ return outdbfilename
		
	##### end of encapsulated subfunctions ################################
	sys.stderr.write(" reading nodes.dmp\n")
	taxdict, lca_walk_tree = read_nodesdmp(os.path.join(download_dir, "nodes.dmp"))
	#print(lca_walk_tree)
	sys.stderr.write(" reading names.dmp\n")
	taxdict = read_namesdmp(os.path.join(download_dir, "names.dmp"), taxdict)
	#print(lca_walk_tree)
	sys.stderr.write(" adding level-info\n")
	lca_walk_tree = add_levelinfo(taxdict, lca_walk_tree)
	#print(lca_walk_tree)
	sys.stderr.write(" created! now only have to write to file...\n")
	taxdictjson_file = taxdict2json(taxdict, download_dir)
	print("taxdict")
	#print(taxdict['375451'])
	print("----------")
	print("lca_walk_tree")
	print(lca_walk_tree)
	lca_paths_file = build_lca_db(lca_walk_tree, download_dir)
	return taxdictjson_file, lca_paths_file 

def json_taxdb_from_kronadb(kronadb):
	raise Exception("This function does not exist yet")
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
		if os.path.exists(f + ".md5"):
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
	def __init__(self, acc2taxid_lookupfile, taxdbfile = None, lca_pathsfile = None):
		#self.taxdict = self.read_taxddbfile(taxdbfile) #todo: write tis!
		self.acc_lookup_handle = openfile(acc2taxid_lookupfile)
		self.acc_lookup_handle_filesize = self.acc_lookup_handle.seek(0,2) #jump to end of file and give bytesize (alternative to "os.path.getsize()")
		if taxdbfile != None:
			try:
				self.read_taxdb(taxdbfile)
				if lca_pathsfile == None:
					self.read_lca_paths(os.path.join(os.path.dirname(taxdbfile), lcawalkdb_outfilebasename))
				else:
					self.read_lca_paths(os.path.join(lca_pathsfile))
			except Exception as e: #TODO: replace this with the specific exception throen when there was actually an problem parsing the file as json
				sys.stderr.write("\n{}\n".format(e, traceback.print_exc()))
				sys.stderr.write("\nperhabs the taxdbfile is not in json-format? Assuming a krona-taxonomydb and trying to covert it to json\n")
				self.taxdbfile = json_taxdb_from_kronadb(taxdbfile)
				self.read_taxdb(self.taxdbfile)
		else:
			self.taxdict = None
	
	def read_lca_paths(self, lca_pathsfile):
		infile = openfile(lca_pathsfile)
		#only two lines are recognized. if there is anything else, it will be ignored
		self.walk_list = infile.readline().strip().split("\t")
		self.depth_list = [ int(l) for l in infile.readline().strip().split("\t") ]
		
	
	def get_lca(self, taxA, taxB):
		indexA = self.walk_list.index(taxA)
		indexB = self.walk_list.index(taxB)
		walk_slice = self.walk_list[min([indexA, indexB]):max([indexA, indexB])]
		depth_slice = self.depth_list[min([indexA, indexB]):max([indexA, indexB])]
		print("Walk slice:")
		print(walk_slice)
		slice_tuples = sorted([ (w,d) for w,d in zip(walk_slice, depth_slice) ], key = lambda x:x[1])
		lca = slice_tuples[0][0]
		print("LCA of '{}' & '{}' is '{}'".format(taxA, taxB, lca))
			
	def read_taxdb(self, taxdbfile):#needs to be json format. Krona taxdbs need to be converted to this format first, using the kronadb2json function above
		import json
		# ~ def keystoint(x): #sumb json format does not allow integers as keys (WHYYY=!) and converts them to strings during dump (WHYY?). doing this in the hope this saves RAM. otherwise just go with stupid strings...
			# ~ replace_dict = {}
			# ~ for k, v in x.items():
				# ~ if type(v) == dict:
					# ~ replace_dict[int(k)] = v
				# ~ else:
					# ~ replace_dict[k] = v
			# ~ return replace_dict
		infile = openfile(taxdbfile)
		# ~ self.taxdict = json.load(infile, object_hook = keystoint) #converting to int does NOT save much RAM BUT costs a LOT of speed!
		self.taxdict = json.load(infile)
		#print(list(self.taxdict.keys()))
		infile.close()
	
	def acc2taxid(self, queryacc,start = 0):
		#for using binary search on a simple sorted textfile #reminer to self: do NOT use compressed acc2taxid_lookupfile. It increases time for binary search ca 200x!
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
	
	def taxid2taxpath(self, taxid, fullpath = True, unofficials = True): #may skip the outformat and return all levels as tuples (taxname, taxid, rank). MAy change fullpath default to False AFTER i checked how to best deal with "unofficial candidate phyla"
		# ~ id_candidate_phyla = 1783234
		# ~ id_bacteria_incertae_sedis = 2323 
		assert self.taxdict != None, "\nError in taxid2taxpath: you must provide a taxdb-file\n"
		assert isinstance(fullpath, bool), "\nError in taxid2taxpath: 'fullpath' must be either True or False\n"
		assert isinstance(unofficials, bool), "\nError in taxid2taxpath: 'unofficials' must be either True or False\n"

		taxpath = []
		# ~ is_candidate_phylum = False
		# ~ is_incertae_sedis = False
		
		official_phylum_level_set = False
		placeholder_phylum = None
		placeholder_phylum_listindex = None
		#print("0000000000000")
		#print(taxid)
		#print(type(taxid))
		#sys.stderr.write("\n______\n{}\n".format(taxid))
		#loopcounter = 0
		while int(taxid) > 1: #assuming ALL taxpaths lead down to "root" (taxid=1); otherwise implement a maximum iteration counter
			#loopcounter += 1
			tax = self.taxdict[taxid]
			taxname = tax["taxname"]
			taxrank = tax["rank"]
			sys.stderr.write("   --> taxname: {} --> rank: {} \n".format(taxname, taxrank)) # somehow all ranks turn out to be 0?? look into that!
			sys.stderr.flush()
			taxparent = tax["parent"]
			#workaround for candidate phyla unrecognized by ncbi taxonomy (when limiting filtering to major ranks such as phylum):
			#will probably drop this here and integrate it in the LCA portion instead, because most candidate phyla actually have an official "phylum rank" and i just want to make sure the "candidate phyla" info is not lost, when the LCA ends up below that rank
			if  unofficials and "Candidate phyla" in taxname and not official_phylum_level_set:
				# ~ is_candidate_phylum = True
				placeholder_phylum = "Candidate phylum " + taxpath[-1][0] #use the lowest level clade in ncbi as "placeholder" for phylum
				placeholder_taxid = taxpath[-1][1]
			if unofficials and "incertae sedis" in taxname and not official_phylum_level_set: 
				# ~ is_incertae_sedis = True
				if placeholder_phylum == None:
					placeholder_phylum = "{} {}".format(taxname, taxpath[-1][0]) #use the lowest level clade in ncbi as "placeholder" for phylum
				else:
					placeholder_phylum = "{} {}".format(taxname, placeholder_phylum)
					placeholder_taxid = taxpath[-1][1]
				if placeholder_phylum_listindex != None: #make sure that only the lowest level instance of "incertae sedis" is interpreted as phylum
					taxpath.pop(placeholder_phylum_listindex)
				placeholder_phylum_listindex = len(taxpath)
				taxpath.append((placeholder_phylum, placeholder_taxid, 20))	
			if taxrank == 20: #safeguard to make sure only one phylum-level entry is in taxpath even when looking at current "inofficials", in case ncbi Taxonomy changes e.g. in regard to "Bacteria candidate phyla"
				phylum_level_set = True
				if placeholder_phylum_listindex != None: #if an placeholder-phylum was set BUT now we find an official ncbi-taxonomy-recognized phylum, delete the placeholder
					taxpath.pop(placeholder_phylum_listindex)
			#end of workaround for candidate phyla. may likely drop the above portion here, and instead adapt it for the LCA portion later? ALthough this is probably mostly used for filtering anyay. so would fit better here...?
			taxpath.append((taxname, taxid, taxrank)) 	
			#if taxparent == 620:
			#	sys.stderr.write("\n\nWTF: ({}, {}, {})".format(taxname, taxid, taxrank)) 		
			taxid = taxparent

			
		if fullpath:
			return list(reversed(taxpath))
		else:
			return list(reversed([ t for t in taxpath if t[2] > 0 ])) #only ranks with indices larger than zero == the 7 official ranks
			

################################################################		

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
	
def _test_lookup():#check how fast accession2taxid lookup actually is for different numbers of input data
	import time
	sys.stderr.write("\ncreating database_object...\n")
	db = taxdb("sorted_shit_yeah.db")
	for infilename in ["100.acc", "1000.acc", "100000.acc"]:
		infile=openfile(infilename)
		acclist=[acc.strip() for acc in infile]
		infile.close()
		sys.stderr.write("\n\n{}\n{}\n".format("-"*20, infilename))
		start_time = time.time()
		taxiddict = db.acclist2taxiddict(acclist)
		stop_time = time.time()
		sys.stderr.write("  this took {} seconds ---\n".format(stop_time - start_time))
		outfile = openfile(infilename + "_taxids.tsv", "wt")
		outfile.write("{}".format("\n".join(["{}\t{}".format(a, taxiddict[a]) for a in acclist] )))
		outfile.close()

def _test_taxpath():#assumes "nodes.dmp" and "names.dmp" "sorted_shit_yeah.db" and "100.acc" are in current working directory
	import time
	sys.stderr.write("\ncreating json_taxdb from dmp-files\n")
	start_time = time.time()
	lca_and_json_taxdb_from_dmp()
	stop_time = time.time()
	sys.stderr.write("  --> Done. This took {} seconds\n".format(stop_time - start_time))

	sys.stderr.write("\nrcreating db-object and readig taxdb from json_file\n")	
	start_time = time.time()	
	db = taxdb("sorted_shit_yeah.db", taxdbfile = taxdb_outfilebasename)
	stop_time = time.time()
	sys.stderr.write("  --> Done. This took {} seconds ---\n".format(stop_time - start_time))		
	
	accfile = "100.acc"
	infile = openfile(accfile)
	acclist = [acc.strip() for acc in infile]
	infile.close()
	
	sys.stderr.write("\ngetting taxids for sample-accessions\n")
	start_time = time.time()
	taxiddict = db.acclist2taxiddict(acclist)
	stop_time = time.time()
	sys.stderr.write("  --> Done. This took {} seconds ---\n".format(stop_time - start_time))	
	
	sys.stderr.write("\nNow getting full taxpath for each taxid\n")
	outlines = []
	start_time = time.time()
	for acc in acclist:
		#print(acc)
		#print(taxiddict[acc])
		fuck = db.taxid2taxpath(taxiddict[acc])
		#print("-"*60)
		#print(fuck)
		#print("====")
		pathstring="\t".join([";".join([str(y) for y in x]) for x in fuck] )
		outlines.append("{}\t{}".format(acc, pathstring))
	stop_time = time.time()
	sys.stderr.write("  --> Done. This took {} seconds ---\n".format(stop_time - start_time))
	
	sys.stderr.write("\nwriting results to file...\n")
	outfile = openfile("testfulltaxpath.out.tsv", "wt")
	outfile.write("\n".join(outlines))

def test_opentaxdbspeed():
	sys.stderr.write("\nrcreating db-object and readig taxdb from json_file\n")	
	start_time = time.time()	
	db = taxdb("dummy.tab", taxdbfile = taxdb_outfilebasename)
	stop_time = time.time()
	sys.stderr.write("  --> Done. This took {} seconds ---\n".format(stop_time - start_time))		

def test_lcawalk():
	import time
	sys.stderr.write("\ncreating json_taxdb from dmp-files\n")
	start_time = time.time()
	taxdictjson_file, lca_paths_file = lca_and_json_taxdb_from_dmp()
	stop_time = time.time()
	sys.stderr.write("  --> Done. This took {} seconds\n".format(stop_time - start_time))
	sys.stderr.write("\nrcreating db-object and readig taxdb from json_file\n")	
	start_time = time.time()	
	db = taxdb(taxdictjson_file, taxdbfile = taxdb_outfilebasename)
	stop_time = time.time()
	sys.stderr.write("  --> Done. This took {} seconds ---\n".format(stop_time - start_time))		
	tax1 = sys.argv[1]
	tax2 = sys.argv[2]
	db.get_lca(tax1,tax2)

if __name__ == '__main__':	
	#_test_lookup()
	#~ _test_taxpath()
	#~ test_opentaxdbspeed()
	test_lcawalk()
