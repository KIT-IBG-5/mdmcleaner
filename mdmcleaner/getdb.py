#!/usr/bin/env python


import os
import sys
import time
import traceback
from mdmcleaner import misc
import subprocess
from mdmcleaner.misc import openfile
import locale #hopefully this ensures that sorts and string comparisons behave like during lookuptable-generation
locale.setlocale(locale.LC_ALL, "C") #hopefully this ensures that sorts and string comparisons behave like during lookuptable-generation

'''
module for downloading and parsing the ncbi taxonomy databases for MDMcleaner work in progress
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

rank2index = { "no rank" : 0, \
				"superkingdom" : 10, \
				"phylum" : 20, \
				"class" : 30, \
				"order" : 40, \
				"family" : 50, \
				"genus" : 60, \
				"species" : 70, \
				"ignored rank" : -1} # using increments of 10 in case i want to use the indermediate ranks (e.g. subfamily) at some later point also
				#rank "root" does not exist (as previously planned. Instead checking for taxid=1 (= root)
				#todo: conflicting naming scheme "superkingdom"  and "domain". Fix this!
				#todo: switched order of keys (moved "ignored rank" from second, to last position). Keep an eye out on whether this breaks something somewehere(although relying on exact key-order in dicts WOULD be DUMB!)!

index2rank = { rank2index[key] : key for key in rank2index }

#todo: split gtdb protblastdbs into several subdbs that can be blasted in parallel (check out if faster first)
#todo: find common names for gtdb and ncbi dbs to simplify things
#todo: find actual names for ncbi dbs
dbfiles = { "gtdb" : {	"protblastdbs" : ["gtdbplus_protdb.dmnd"], \
						"nucblastdbs" : ["concat_refgenomes", "SILVA_138.1_SSURef_NR99_tax_silva", "SILVA_138.1_LSURef_NR99_tax_silva"] ,\
						"ssu_nucblastdbs" : ["concat_refgenomes", "SILVA_138.1_SSURef_NR99_tax_silva"], \
						"lsu_nucblastdbs" : ["concat_refgenomes", "SILVA_138.1_LSURef_NR99_tax_silva"], \
						"genome_nucblastdbs" : ["concat_refgenomes"], \
						"mdmdbs" : ["gtdb_all.accession2taxid.sorted", "gtdb_taxonomy_br.json.gz", "gtdb_lcawalkdb_br.db"] },\
			"ncbi" : {  "protblastdbs" : ["nr"], \
						"nucblastdbs" : ["nt"], \
						"ssu_nucblastdbs" : ["nt"], \
						"lsu_nucblastdbs" : ["nt"], \
						"genome_nucblastdbs" : ["nt"], \
						"mdmdbs" : [ "ncbi_accession2taxid", "ncbi_taxonomy_br.json.gz", "ncbi_lcawalkdb_br.db"] } }


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

def dict2jsonfile(taxdict, outdbfilename): #TODO: move this to misc?
	import json
	outfile = openfile(outdbfilename, 'wt')
	json.dump(taxdict, outfile)
	outfile.close()
	return outdbfilename

def jsonfile2dict(jsonfile):#needs to be json format. Krona taxdbs need to be converted to this format first, using the json_taxdb_from_kronadb function below  #TODO: move this to misc?
	import json
	#print("loading {}".format(jsonfile))
	infile = openfile(jsonfile)
	outdict = json.load(infile)
	infile.close()
	return outdict
		
def json_taxdb_from_kronadb(kronadb):
	raise Exception("This function does not exist yet")
	infile = openfile(kronadb)
	for line in infile:
		pass #todo: finish this sometime

def _create_sorted_acc2taxid_lookup(acc2taxidfilelist, acc2taxid_outfilename): 
	'''
	for creating an alphabetically sorted list of acessions and their corrsponding taxids
	for use with binary search
	TODO: try making this presort multiple input files in paralall (if multiple cpus are specified)
	'''
	import subprocess
	# ~ import time
	sys.stderr.write("\n\tcreating {}...\n".format(acc2taxid_outfilename))
	# ~ start = time.time() #todo: for debugging
	#import shlex #allow string splitting accoring to Shell -like syntax
	#presortcmd = "zcat {infile} | cut -f 2,3| grep -v accession | sed 's#\.[0-9]*##'| sort > {outfile}" #using shell commands probably way faster than anything i can do in pure python
	#NOTE: the presortcmd partially renames the accession (remoes version number). This may have been necessary for ncbi, but does not work for gtdb data. 
	#      so am removing this for now (will find a way later to flexibly deal with it for ncbi data)
	#TODO: find out if removing verison numbers from accessions actually necessary for ncbi data. Find a flexible workaround that works for ncbi AND gtdb data
	presortcmd = "zcat {infile} | cut -f 2,3| grep -v accession | grep -v -P '^$'| env LC_ALL=C sort > {outfile}" #IMPORTANT! removes a sed substitution in accession field. check  if still works for nucleotide basts.  additional note: using shell commands probably way faster than anything i can do in pure python
	#ALSO IMPORTANT: added removal of empty lines in the above command. This is because some input files contain empty lines. This resulted in empty lines in the sorted lookupfile, resulting in broken lookups. This is fixed now.
	#ALSO IMPORTANT: added "env LC_ALL=C before the sort command, in a desperate attempt to ensure same localse seetings for creating and using the accession index. Sort behaves differently based on locale settings and that can cause problems
	finalsortcmd = "env LC_ALL=C sort -m {filelist} > {finaldb}"
	tempfilelist = []
	for f in acc2taxidfilelist:
		if f.endswith(".gz"):
			tempfile = f[:-3]
		tempfile = tempfile + ".sorted"
		sys.stderr.write("\textracting and presorting {}\n".format(f))
		sout, serr = subprocess.Popen(["bash", "-c", presortcmd.format(infile=f, outfile=tempfile)], stdout=subprocess.PIPE).communicate()
		if serr != None:
			raise RuntimeError("...extraction exited with Error:\n{}\n".format(serr))
		os.remove(f) 
		if os.path.exists(f + ".md5"):
			os.remove(f + ".md5")
		tempfilelist.append(tempfile)
	sys.stderr.write("\tcombining sorted tempfiles: {}\n".format(", ".join(tempfilelist)))
	sout,serr = subprocess.Popen(["bash", "-c", finalsortcmd.format(filelist=" ".join(tempfilelist), finaldb=acc2taxid_outfilename)], stdout=subprocess.PIPE).communicate()
	if serr != None:
			raise RuntimeError("...extraction exited with Error:\n{}\n".format(serr))
	#now clean up:
	# ~ end = time.time()
	# ~ print("this took : {}".format(end - start)) #todo: for debugging
	# ~ sys.stdout.flush()	#todo: for debugging
	sys.stderr.write("\tremoving temporary downloads")
	for f in tempfilelist:
		os.remove(f)
	return acc2taxid_outfilename

class taxdb(object):
	def __init__(self, configs): # acc2taxid_lookupfile, taxdbfile = None, lca_pathsfile = None): #todo: create and read a "config file" to get file-locations from
		#self.taxdict = self.read_taxddbfile(taxdbfile) #todo: write tis!
		self.set_db_attributes(configs)
		self.check_db_folder(configs)
		self.read_db_versions()
		self.acc_lookup_handle = misc.openfile(self.acc2taxid_lookupfile)
		self.acc_lookup_handle_filesize = self.acc_lookup_handle.seek(0,2) #jump to end of file and give bytesize (alternative to "os.path.getsize()")
		self.read_taxdb(self.taxdbfile)
		self.read_lca_paths(self.lca_pathsfile)

	def set_db_attributes(self, configs, db_basedir = None):
		from mdmcleaner import read_gtdb_taxonomy #todo: after reimplementing optional ncbi taxonomy, put both into the same module and import THAT here
		if db_basedir != None:
			self.db_basedir = db_basedir
		else:
			self.db_basedir = configs.settings["db_basedir"][0]
		self.db_type = configs.settings["db_type"][0]
		self.dbpath = os.path.join(self.db_basedir, self.db_type)
		 
		self.acc2taxid_lookupfile = os.path.join(self.dbpath, dbfiles[self.db_type]["mdmdbs"][0])
		self.taxdbfile = os.path.join(self.dbpath, dbfiles[self.db_type]["mdmdbs"][1])
		self.lca_pathsfile = os.path.join(self.dbpath, dbfiles[self.db_type]["mdmdbs"][2])
		
		self.nucdbs_ssu_rRNA = [os.path.join(self.dbpath, x) for x in dbfiles[self.db_type]["ssu_nucblastdbs"]]
		self.nucdbs_lsu_rRNA = [os.path.join(self.dbpath, x) for x in dbfiles[self.db_type]["lsu_nucblastdbs"]]
		self.nucdbs_genome = [os.path.join(self.dbpath, x) for x in dbfiles[self.db_type]["genome_nucblastdbs"]]
		self.nucdbs_all = [os.path.join(self.dbpath, x) for x in dbfiles[self.db_type]["nucblastdbs"]]
		
		self.protdbs_all = [os.path.join(self.dbpath, x) for x in dbfiles[self.db_type]["protblastdbs"]]
		
		self.versionfile = os.path.join(self.dbpath, "DB_versions.txt")
		self.final_progress_marker = os.path.join(self.dbpath, "progress_step{}.json".format(read_gtdb_taxonomy._progress_steps["finished"][-1]))#todo: make a "get_progress_marker_filename()" function in read_gtdb_taxonomy.py?

	def check_db_folder(self, configs):
		if not (os.path.exists(self.dbpath) and os.path.isdir(self.dbpath)):
			if self.db_basedir.endswith(self.db_type) and os.path.exists(self.db_basedir):
				sys.stderr.write(	"\n"+"!"*150+"\n"
									"WARNING: It seems that 'db_basedir' is specified up to the 'gtdb' subfolder!\n"
					                "        in fact you should only specify the directory CONTAINING the 'gtdb' subfolder (and then specify 'gtdb' as 'db_type')\n"
									"        assuming for now that the specified path to 'gtdb' is correct, but if it is you should correct your configuration with the following command\n"
									"        for local config file:\n"
									"            mdmcleaner.py set_configs --db_basedir {actualpath}\n"
									"        for global config file:\n"
									"            sudo mdmcleaner.py set_configs -s global --db_basedir {actualpath}\n".format(actualpath=os.path.dirname(self.db_basedir)) )
				sys.stderr.write(	"!"*150+"\n")
				self.set_db_attributes(configs, db_basedir = os.path.dirname(self.db_basedir))
			else:
				sys.exit("\n\nERROR: can't find specified database folder '{dbpath}'!\nyou may need to download and create the database with 'mdmcleaner.py download_db -o {dbpath}'".format(dbpath=self.db_basedir))
		if not os.access(self.dbpath, os.R_OK):
			sys.exit("\n\nERROR: insufficient permissions to read from db_dir: '{}'\n".format(self.dbpath))
		if not os.path.exists(self.final_progress_marker):
			sys.exit("\n\nERROR: Database creation is not complete! Please run 'mdmcleaner.py download_db -o {dbpath}' to finish it!\n".format(dbpath=self.dbpath))
		for f in [self.acc2taxid_lookupfile, self.taxdbfile, self.lca_pathsfile, self.versionfile]:
			assert os.path.exists(f), "\n\nERROR: can't find file '{}' in '{}'! Reference database is not complete! Please run 'mdmcleaner.py download_db -o {dbpath}' to finish it!\n".format(os.path.basename(f), self.db_basedir)
	
	def read_db_versions(self):
		self.versions = {}
		with misc.openfile(self.versionfile) as infile:
			for line in infile:
				tokens = line.strip().split("=")
				if len(tokens) != 2:
					continue
				db = tokens[0].strip()
				version = tokens[1].strip()
				self.versions[db] = version			
		self.print_db_versions()
	
	def print_db_versions(self):
		sys.stderr.write("\n\tDatabase versions:\n")
		for db in self.versions:
			sys.stderr.write("\t\t{} = {}\n".format(db, self.versions[db]))
		sys.stderr.write("\n")	
		
	def read_lca_paths(self, lca_pathsfile):
		infile = misc.openfile(lca_pathsfile)
		#only two lines are recognized. if there is anything else, it will be ignored
		self.walk_list = infile.readline().strip().split("\t")
		self.depth_list = [ int(l) for l in infile.readline().strip().split("\t") ]
		infile.close()
		
	def get_strict_pairwise_lca(self, taxA, taxB):
		indexA = self.walk_list.index(taxA)
		indexB = self.walk_list.index(taxB)
		walk_slice = self.walk_list[min([indexA, indexB]):max([indexA, indexB])+1]
		depth_slice = self.depth_list[min([indexA, indexB]):max([indexA, indexB])+1]
		# ~ print("taxA : {}, taxB : {}".format(taxA, taxB))
		# ~ print("Walk slice:")
		# ~ print(walk_slice)
		# ~ print("Depth_slice")
		# ~ print(depth_slice)
		# ~ print("--------------")
		slice_tuples = sorted([ (w,d) for w,d in zip(walk_slice, depth_slice) ], key = lambda x:x[1])
		lca = slice_tuples[0][0]
		#print("LCA of '{}' & '{}' is '{}'".format(taxA, taxB, lca))
		return lca
			
	def read_taxdb(self, taxdbfile):#needs to be json format. Krona taxdbs need to be converted to this format first, using the kronadb2json function above
		self.taxdict = jsonfile2dict(taxdbfile)

	def _delmegetline(self, searchterm):# only for debugging
		getlinecmd = ["grep", "-n", searchterm, self.acc2taxid_lookupfile]
		myproc = subprocess.run(getlinecmd, stdout = subprocess.PIPE, stderr = subprocess.PIPE, text = True)
		linenumber = myproc.stdout.split(":")[0]
		return linenumber
	
	def acc2taxid(self, queryacc,start = 0):		
		#for using binary search on a simple sorted textfile #reminer to self: do NOT use compressed acc2taxid_lookupfile. It increases time for binary search ca 200x!
		import locale
		locale.setlocale(locale.LC_ALL, "C")
		stop = self.acc_lookup_handle_filesize
		# ~ finalstop = stop #todo: only for debugging
		# ~ totallinecount = 35602626 #todo: only for debugging
		subjectacc = None
		
		while start < stop: 
			currentpos = int((start + stop) / 2)
			self.acc_lookup_handle.seek(currentpos, 0)

			if currentpos != start:
				self.acc_lookup_handle.readline()
			tokens =  self.acc_lookup_handle.readline().strip().split("\t")
			# ~ print("{}  -  {}   /   {}".format(start, stop, finalstop))
			# ~ startline=(start/finalstop) * totallinecount
			# ~ stopline=(stop/finalstop) * totallinecount
			# ~ currentline = (currentpos/finalstop) * totallinecount
			# ~ print("--> currentpos: {}".format(currentpos))
			# ~ print("---> currentline: {}".format(self._delmegetline(tokens[0])))
			# ~ print(tokens)
			# ~ print("queryacc: {}".format(queryacc))
			# ~ print("startline={}, stopline={}, currline={}".format(startline, stopline, currentline))
			# ~ print("---")
			if tokens == [""]:
				sys.stderr.write("\nWARNING: do not recognize accession '{}'\n".format(queryacc))
				break
			subjectacc = tokens[0]
			subjecttaxid = tokens[1]
			if locale.strcoll(subjectacc, queryacc) == 0:
				return subjecttaxid, start
			if locale.strcoll(subjectacc, queryacc) >= 1:
				stop = currentpos
			else:
				start = currentpos
		return None, 0
	
	def acclist2taxiddict(self, queryacclist): #may still reconsider whether i actually want to output a dictionary or not
		start = 0
		outdict = {}
		for queryacc in sorted(queryacclist): #VERY important that queries are sorted as well!
			outdict[queryacc], start = self.acc2taxid(queryacc, start)
		return outdict
	
	def isnot_bacteria(self, taxid): #todo. double check if this also still works when using ncbi data!
		tp = self.taxid2taxpath(taxid)
		if len(tp) >= 2 and tp[1][0] == "Bacteria": #todo: make sure this simplified index-based lookup still works if using taxpaths that include non-official ranks (but probably should, at least for "root" and "domain/superkingdom"-levels
			return False
		return True
	
	def isnot_archaea(self, taxid): #todo. double check if this also still works when using ncbi data!
		tp = self.taxid2taxpath(taxid)
		if len(tp) >= 2 and tp[1][0] == "Archaea": #todo: make sure this simplified index-based lookup still works if using taxpaths that include non-official ranks (but probably should, at least for "root" and "domain/superkingdom"-levels
			return False
		return True	
		
	def is_viral(self, taxid): #todo. double check if this also still works when using ncbi data!
		tp = self.taxid2taxpath(taxid)
		if len(tp) >= 1 and tp[0][0] in ["Viruses", "Virus"]: #todo: make sure this simplified index-based lookup still works if using taxpaths that include non-official ranks (but probably should, at least for "root" and "domain/superkingdom"-levels
			return True
		return False
	
	def is_eukaryote(self, taxid):
		tp = self.taxid2taxpath(taxid)
		if len(tp) >= 1 and tp[1][0] == "Eukaryota": #todo: make sure this simplified index-based lookup still works if using taxpaths that include non-official ranks (but probably should, at least for "root" and "domain/superkingdom"-levels
			return True
		return False	

	def get_domain_phylum(self, taxid):
		'''
		returns a tuple containing only the domain and phylum designation associated with the given taxid (if any)
		'''
		domain_phylum = [None, None]
		tp = self.taxid2taxpath(taxid)[1:3]
		for i in range(len(tp)):
			domain_phylum[i] = tp[i][0]
		return tuple(domain_phylum)
		

	def get_specific_taxlevel_subtaxid(self, taxid, taxlevel="domain", returnvalue="taxid"):
		"""
		accepts a taxid and a taxlevel (one of ["root", "domain", "phylum", "class", "order", "family", "genus", "species"]) as input.
		looks up the taxpath for the given taxid and returns the ("sub-")taxid or taxname at the given taxlevel (default="domain"), if present.
		If the taxpath does not contain the given taxlevel (e.g. when requesting a "subtaxid" at genus level, but giving a taxid of a phylum), None is returned
		"""
		from mdmcleaner import lca
		assert taxlevel in lca.taxlevels, "\nERROR: taxlevel must be one of {}\n".format(lca.taxlevels[:2] + ["superkingdom"] + lca.taxlevels[2:])
		assert returnvalue in ["taxid", "taxname"], "ERROR: returnvalue mist bei either \"taxid\" or \"taxname\"\n"
		outtaxid = None
		taxpath = self.taxid2taxpath(taxid)
		if taxlevel == "domain":
			taxlevel = "superkingdom"
		if taxlevel == "root":
			outtaxid = taxpath[0][1]
			outtaxname = taxpath[0][0]
		else:
			for t in taxpath[1:]:
				outtaxid = t[1]
				outtaxname = t[0]
				if self.taxid2taxlevel(outtaxid) == taxlevel:
					break
				outtaxid, outtaxname = None, None
		if returnvalue == "taxid":
			return outtaxid
		elif returnvalue == "taxname":
			return outtaxname
		
		
	def taxid2taxname(self, taxid):
		return self.taxdict[taxid]["taxname"]
	
	def taxid2pathstring(self, taxid):
		"""returns semicolon-seperated string-representation of the full taxpath for a given taxid
		"""
		if taxid:
			return ";".join([ p[0] if type(p) == tuple else p.taxname for p in self.taxid2taxpath(taxid)]) #unnecessarily flexible, as taxid2taxpath returns tuples, not namedtuples. But may avoid confugion in the furture (when i switch to ONE sytem for storing taxpaths)


	def taxid2taxpath(self, taxid, fullpath = True, unofficials = False): #may skip the outformat and return all levels as tuples (taxname, taxid, rank). MAy change fullpath default to False AFTER i checked how to best deal with "unofficial candidate phyla"
		#print("Hi, taxid2taxpath here. I got this: '{}'".format(taxid))
		def notroot(taxid):
			if taxid == "root":
				return False
			else:
				try:
					if int(taxid) <= 1:
						return False
				except ValueError:
					return True
			return True
		# todo: switch tuples to namedtuples
		"""
		Returns a list of tuples representing the taxonomic Path for a given Taxid.
		Each tuple represents a single taxonomic level and consists of three values: ( taxon name [String] , taxon ID [String] , rank [Integer])
		The returned path will include all ranks,including minor intermediate ranks such as "Subfamily" if fullpath == True, otherwise it will consistonly of the 7 major ranks.
		When using the ncbi taxonomic system, set unofficials==True to include unofficial ranks such as Candidate phyla or "Incertae sedis" taxa. 
		When using the gtdb taxonomy, set unofficials==False.
		"""
		
		#todo: veryfy that this works for ncbi as well as for gtdb taxonomy-dbs
		# ~ id_candidate_phyla = 1783234
		# ~ id_bacteria_incertae_sedis = 2323 
		assert self.taxdict != None, "\nError in taxid2taxpath: you must provide a taxdb-file\n"
		assert isinstance(fullpath, bool), "\nError in taxid2taxpath: 'fullpath' must be either True or False\n"
		assert isinstance(unofficials, bool), "\nError in taxid2taxpath: 'unofficials' must be either True or False\n"
		if not taxid:
			return
		taxpath = []
		# ~ is_candidate_phylum = False
		# ~ is_incertae_sedis = False
		
		official_phylum_level_set = False
		placeholder_phylum = None
		placeholder_phylum_listindex = None
		while notroot(taxid): #assuming ALL taxpaths lead down to "root" (taxid=1); otherwise implement a maximum iteration counter
			tax = self.taxdict[taxid]
			taxname = tax["taxname"]
			taxrank = tax["rank"]
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
			taxid = taxparent

			
		if fullpath:
			return list(reversed(taxpath))
		else:
			return list(reversed([ t for t in taxpath if t[2] > 0 ])) #only ranks with indices larger zero == the 7 official ranks


			
	def taxid2taxlevel(self, taxid):
		"""
		supposed to return the taxlevel (domain, phylum, etc) of a given taxid, if possible
		"""
		rankindex = self.taxdict[taxid]["rank"]
		return index2rank[rankindex]

	def taxids2contradicting_taxpaths(self, taxA, taxB):
		"""
		takes two taxids, and checks if their taxpaths contradict each other
		if the taxpaths contradict, it returns the major contradicting taxlevel/rank, or it returns "minor rank below <last matching official taxlevel>" if contradicting at minor rank
		returns None if no contradiction is found.
		unequal lengths of taxpaths are NOT considered a contradicion (taxpaths will only be compared up to the shorter length of the two taxpaths)
		"""
		taxpathA = self.taxid2taxpath(taxA)
		taxpathB = self.taxid2taxpath(taxB)
		return self.contradicting_taxpaths(taxpathA, taxpathB)

	def contradicting_taxpaths(self, taxpathA, taxpathB):
		"""
		takes two taxpaths and checks for contradictions
		if the taxpaths contradict, it returns the major contradicting taxlevel/rank, or it returns "minor rank below <last matching official taxlevel>" if contradicting at minor rank
		returns None if no contradiction is found.
		unequal lengths of taxpaths are NOT considered a contradicion (taxpaths will only be compared up to the shorter length of the two taxpaths)
		"""
		found_contradiction = False
		last_common_rankindex = None
		for ta, tb in zip(taxpathA, taxpathB):
			if ta != tb:
				found_contradiction = True
				if ta[2] < 0 and tb[2] < 0: 
					continue #if contradicting taxlevel is an minor rank, such as "subfamily"
				else:
					return index2rank(min([ x for x in [ta[2], tb[2]] if x >= 0 ])) 
			else:
				last_common_rankindex = ta[2]
		if found_contradiction:
			return "minor rank below {}".format(index2rank(last_common_rankindex))
		
	def _gtdb_refseq_or_silva(self, acc):
		"""
		determines which db an accession-nr was from (gtdb, silva or refseq), based on the mdmcleaner-internal formatting/style of the accession-nr
		only tested for the gtdb-based workflow (which includes eukaryotic and viral protein sequences from RefSeq)
		primarily used to determine gtdb/silva refDB.ambiguities within the mdmcleaner workflow.
		accepts an accession-id as input, returns "silva" for SILVA-style accessions, "gtdb" for gtdb-style and "refseq_eurvircat" for refseq. returns None if none of the patterns match
		"""
		import re
		slv_pattern = ("silva", "^\w+\.\d+\.\d+$")
		gtdb_contig_pattern = ("gtdb", "^\w+\.\d+_\w+\.*\d+$")
		gtdb_genome_pattern = ("gtdb", "^(RS|GB)_GC[AF]_\d+\.\d+$")
		refseq_euvircat_pattern = ("refseq", "^[A-Z]{2}_\d+\.\d+$")
		
		for p in [slv_pattern, gtdb_contig_pattern, gtdb_genome_pattern, refseq_euvircat_pattern]:
			if re.search(p[1], acc):
				return p[0]
				
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
	

def test_lcawalk():
	acc2taxid_lookupfile = "tempdir/gtdb_all.accession2taxid.sorted"
	taxdictjson_file = "tempdir/gtdb_taxonomy_br.json.gz"
	lcawalk_file = "tempdir/gtdb_lcawalkdb_br.db"
	db = taxdb(acc2taxid_lookupfile, taxdictjson_file, lcawalk_file)		
	tax1 = sys.argv[1]
	tax2 = sys.argv[2]
	db.get_strict_pairwise_lca(tax1,tax2)

if __name__ == '__main__':	
	#_test_lookup()
	#~ _test_taxpath()
	#~ test_opentaxdbspeed()
	test_lcawalk()
