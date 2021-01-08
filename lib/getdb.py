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
# TODO: add flags specifying progress, for better being able to resume where it last stopped

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
#TODO: ADD accession2taxig.FULL.gz (which is now new)!


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
				#rank "root" does not exist (as previously planned. Instead checking for taxid=1 (= root)

index2rank = { rank2index[key] : key for key in rank2index }


ncbi_taxdb_outfilebasename = "ncbi_taxonomy_br.json.gz" #different from krona. additional wasteful "children" info per node (may get rid of it again after perfecting LCA). saved as compressed file because uncompression does not increase loading time much but decreases file size a LOT
gtdb_taxdb_outfilebasename = "gtdb_taxonomy_br.json.gz"
ncbi_acc2taxid_outfilebasename = "all.accession2taxid.sorted" #using the same as krona
gtdb_acc2taxid_outfilebasename = "gtdb_all.accession2taxid.sorted"
ncbi_lcawalkdb_outfilebasename = "ncbi_lcawalkdb_br.db"
gtdb_lcawalkdb_outfilebasename = "gtdb_lcawalkdb_br.db"


#TODO: actually read these downloaded temfiles and create a database
def calculate_filehash(infile): #TODO probably move to misc.py
	"""
	calculate md5-hash of a file
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

def _download_ncbidb3(dbtype, targetdir="."): #TODO: TEMPORARY! should be replaced by the version used to download the silva, refseq and gtdb datasets
	import wget
	if dbtype not in _dbsource_dict: #usage help
		raise KeyError("\nERROR: unknown dbtype '{}'!\nMust be one of {}\n".format(dbtype, ", ".join(sorted(_dbsource_dict.keys()))))
	
	def checkmd5(infile, md5file):
		import hashlib #for comparing md5-checksums of downloaded databases
		
		#calculate hash of downloaded file:
		md5ist = calculate_filehash(infile)
		
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
			
	# ~ def taxdict2yaml(taxdict, targetdir): #assume that targetdir = downloaddir #optional in cas i decide to switch to yaml
		# ~ import yaml
		# ~ outdbfilename = os.path.join(targetdir, ncbi_taxdb_outfilebasename + ".yaml")
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
	taxdictjson_file = dict2jsonfile(taxdict, os.path.join(download_dir, ncbi_taxdb_outfilebasename))
	print("taxdict")
	#print(taxdict['375451'])
	print("----------")
	print("lca_walk_tree")
	print(lca_walk_tree)
	lca_paths_file = build_lca_db(lca_walk_tree, os.path.join(download_dir, ncbi_lcawalkdb_outfilebasename))
	return taxdictjson_file, lca_paths_file 

def build_lca_db(lca_walk_tree, outfilename, startingnode = "1"): #this default starting node only works for ncbi based taxonomy. need to pass other staring node when using gtdb
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
	walk_list, depth_list = walk(lca_walk_tree, walk_list, depth_list, startingnode)
	
	#outfilename = os.path.join(targetdir, ncbi_lcawalkdb_outfilebasename)
	outfile = openfile(outfilename, "wt")
	outfile.write("\t".join(walk_list) + "\n")
	outfile.write("\t".join([str(l) for l in depth_list]))
	outfile.close()
	return outfilename

def dict2jsonfile(taxdict, outdbfilename): 
	import json
	outfile = openfile(outdbfilename, 'wt')
	json.dump(taxdict, outfile)
	outfile.close()
	return outdbfilename

def jsonfile2dict(jsonfile):#needs to be json format. Krona taxdbs need to be converted to this format first, using the kronadb2json function above
	import json
	infile = openfile(jsonfile)
	outdict = json.load(infile)
	infile.close()
	return outdict
		
def json_taxdb_from_kronadb(kronadb):
	raise Exception("This function does not exist yet")
	infile = openfile(kronadb)
	for line in infile:
		pass #todo: finish this sometime
		

def download_ncbiaccessiondb(targetdir=".", dbmode = "minimal", enforce=False): #TODO: implement filecheck and ask for user-feedback (skippable with "-f" or "-y" argument) if it looks as if a krona-db file could be overwritten
	acc2taxid_outfilename = os.path.join(targetdir, ncbi_acc2taxid_outfilebasename) #change to whatever krona is using
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
		downloaded_filelist.append(_download_db3(w, targetdir)) #todo: switched from download_db2 to download_db3. change to universal download-function (gtdb-stuff) and make sure it works
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
					self.read_lca_paths(os.path.join(os.path.dirname(taxdbfile), ncbi_lcawalkdb_outfilebasename))
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
		self.taxdict = jsonfile2dict(taxdbfile)
	
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

def _test_download3():
	sys.stderr.write("_test_dowload3\n")
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
	db = taxdb("sorted_shit_yeah.db", taxdbfile = ncbi_taxdb_outfilebasename)
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
	db = taxdb("dummy.tab", taxdbfile = ncbi_taxdb_outfilebasename)
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
	db = taxdb(taxdictjson_file, taxdbfile = ncbi_taxdb_outfilebasename)
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
