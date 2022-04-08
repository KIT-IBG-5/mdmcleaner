#!/usr/bin/env python
#TODO: CHANGE THIS WHOLE THING TO USE DOWNLOAD_FUNCTIONS DESIGNED FOR read_gtdb_taxonomy.py!!!!!

""" 
downloading, parsing and modifying ncbi reference and taxonomy databases for MDMcleaner.

Uses the following files of the ncbi taxonomy database
	- "nodes.dmp" (from "taxdump.tar.gz" or "new_taxdump.tar.gz")
	- "names.dmp" (from "taxdump.tar.gz" or "new_taxdump.tar.gz")
	- "prot.accession2taxid" (and, optionally, "dead_prot.accession2taxid.gz")
To decrease memory consumption, taxonomix ranks are converted to integer-indices (although sonly minimal amout of RAM saved this way). #TODO: this has changed.Check and correct this
Taxids were als supposed to be represented as intergers, but this was incompatible with json, slowed down database construction and saved only little RAM --> dropped
These indices can be mapped back and forth using the dictionaries "getdb.rank2index" and/or "getdb.index2rank" #TODO: this has changed.Check and correct this
"""
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


taxdb_outfilebasename = "ncbi_taxonomy_br.json.gz" #different from krona. additional wasteful "children" info per node (may get rid of it again after perfecting LCA). saved as compressed file because uncompression does not increase loading time much but decreases file size a LOT #todo: "_br" still stands for "binrefiner". change that!
acc2taxid_outfilebasename = "all.accession2taxid.sorted" #using the same as krona #todo: "_br" still stands for "binrefiner". change that!
lcawalkdb_outfilebasename = "ncbi_lcawalkdb_br.db"



def _download_ncbidb3(dbtype, targetdir="."): #TODO: TEMPORARY! should be replaced by the version used to download the silva, refseq and gtdb datasets
	import wget
	from mdmcleaner import misc
	if dbtype not in _dbsource_dict: #usage help
		raise KeyError("\nERROR: unknown dbtype '{}'!\nMust be one of {}\n".format(dbtype, ", ".join(sorted(_dbsource_dict.keys()))))
	
	def checkmd5(infile, md5file):
		import hashlib #for comparing md5-checksums of downloaded databases
		
		#calculate hash of downloaded file:
		md5ist = misc.calculate_md5hash(infile)
		
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
	taxdictjson_file = dict2jsonfile(taxdict, os.path.join(download_dir, taxdb_outfilebasename))
	print("taxdict")
	#print(taxdict['375451'])
	print("----------")
	print("lca_walk_tree")
	print(lca_walk_tree)
	lca_paths_file = build_lca_db(lca_walk_tree, os.path.join(download_dir, lcawalkdb_outfilebasename))
	return taxdictjson_file, lca_paths_file 

def download_ncbiaccessiondb(targetdir=".", dbmode = "minimal", enforce=False): #TODO: implement filecheck and ask for user-feedback (skippable with "-f" or "-y" argument) if it looks as if a krona-db file could be overwritten
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
		downloaded_filelist.append(_download_db3(w, targetdir)) #todo: switched from download_db2 to download_db3. change to universal download-function (gtdb-stuff) and make sure it works
	sys.stderr.write("finished downloading accessiondbs for selection: '{}'".format(dbmode))
	#extract and sort them:
	_create_sorted_acc2taxid_lookup(downloaded_filelist, acc2taxid_outfilename)
	return acc2taxid_outfilename
