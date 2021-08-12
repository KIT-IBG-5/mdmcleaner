#!/usr/bin/env python
""" 
creating, reading and parsing of blast files (Diamond and well as blast+).
only creates/handles tabular blast files in "-outfmt 6" format
"""

#todo: assign blasts to contigs-proteins in a nother module.
#todo:   --> in that module create "contigs" classes, that combine all blast hits (and corresponding marker-gene status) for each contig
#todo:		--> in that module combine contig-objects of an assembly dataset into a dataset-object

#todo: use differend score cutoffs for protein/unknown-contig-blasts and 16S/23S-blasts (für proteins accept hits down to 50% of best score, for 16S/23S accept hits only down to 90-95%% of best score)

import sys
import os
assert sys.version_info >= (3, 7), "This module requires python 3.7+! You, however, are running '{}'".format(sys.version)
from collections import namedtuple
import misc
hit = namedtuple('hit', 'accession taxid identity score')

from misc import openfile, run_multiple_functions_parallel

def read_lookup_table(lookup_table, table_format = None): #should be able to read tsv, csv and gff. should output a simple dict with locus_tags as keys and contignames as values
	#consider putting this under "lookup_handler.py"? Probably naaaw...
	#start of subfunctions
	def _read_anything(lookup_table):
		global table_format
		linelist = []
		infile = openfile(lookup_table)
		firstline = True
		for line in infile:
			if line.strip().startswith("#"):
				if firstline and line.startswith("##gff-version"):
					table_format = "gff"
				if line == "##FASTA":
					table_format = "gff"
					break 
				continue
			linelist.append(line)
		return linelist
	
	def _guess_format(linelist):
		"""
		format is assumed to be tsv, if splitting by tab yields exactly two columns per line
		format is assumed to be csv, if splitting by comma yields exactly two columns per line
		otherwise format is assumed to be gff
		"""
		no_tab_fields = linelist[0].split("\t")
		no_comma_fields = linelist[0].split(",")
		if no_tab_fields == 9:
			return "gff"
		if no_tab_fields == 2:
			return "tsv"
		if no_tab_comma_fields == 2:
			return "csv"
		else:
			raise IOError("Can't recognize format of lookup_table! It should be either gff OR a simple tab- or comma seperated table with TWO columns \
			(locus_tag in the first, and contig_name in the second column)")
	
	def _parse_gff(linelist): #there may be some customized gffs out there (e.g. only "CDS"-lines grepped). Does not harm to include them also...
		lookup_dict = {}
		for line in linelist:
			tokens = line.split("\t")
			contig = tokens[0]
			attributes = tokens[9].split[";"]
			for att in attributes:
				if att.startswith("ID="):
					locustag = att[3:]
					break
			lookup_dict[locustag] = contig
		return lookup_dict
		
	def _parse_simpletable(linelist, seperator="\t"):
		lookup_dict = {}
		for line in linelist:
			line.split(seperator)[0] = locustag
			line.split(seperator)[0] = contig
			lookup_dict[locustag] = contig
		return lookup_dict

	#end of subfunctions
	
	possible_table_formats = ["gff", "tsv", "csv", None]
	assert table_format in possible_table_formats, "ERROR: tableformat must be one of {}".format(",".join([str(x) for x in possible_table_formats]))
	linelist = _read_anything(lookup_table)
	if table_format == None:
		table_format = _guess_format(linelist)
	if table_format == "gff":
		return _parse_gff(linelist)
	if table_format("tsv"):
		return _parse_simpletable("\t")
	if table_format("csv"):
		return _parse_simpletable(",")


class blastdata(object): #todo: define differently for protein or nucleotide blasts (need different score/identity cutoffs)
	_blasttsv_columnnames = {"query" : 0, "subject" : 1, "ident" : 2, "alignlen" : 3, "qstart": 6, "qend" : 7, "sstart" : 8, "ssend" : 9, "evalue" : 10, "score" : 11, "contig" : None, "stype" : None, "taxid": None} #reading in all fields, in case functionality is added later that also uses the sstat send etc fields
	# ~ _nuc_stypelist = ["ssu_rRNA", "lsu_rRNA"]
	# ~ _prot_stypelist = ["total", "bac_marker", "prok_marker", "arc_marker"]
	def __init__(self, *blastfiles, max_evalue = None, min_ident = None, score_cutoff_fraction = 0.75, keep_max_hit_fraction = 0.5, keep_min_hit_count = 2, continue_from_json = False, auxilliary = False, seqtype=None): #todo: add a min_score #todo: change score_cutoff_fraction to 0.75?
		#methods:
		#	method 1: filter by query-genes. For each gene remove all hits with score below <score_cutoff_fraction> (default = 0.5) of maximum score for that gene
		#			  from these keep the <keep_max_hit_fraction> of hits (default = 0.5), but keep at least <keep_min_fraction> (default = 2) in every case if possible
		#			  later classify each gene by simple lca
		# note: for 16S/23S the score should be more than 1000 and the identity should be above 90??
		# ~ if from_pickle and len(blastfiles) == 1: #if from_pickle is set, that means input shuld NOT be a list of blastfiles, but a single blastdata-object-pickle-file
			# ~ print("blastdata object already exists! loading from pickle")
			# ~ self.unpickleyourself(blastfiles[0])
		assert seqtype in ["nuc", "prot", None], "\nERROR: seqtype must be either 'nuc', 'prot' or None (if unknown)\n"
		self.seqtype = seqtype #todo: do something with this (e.g. set some cutoff)
		self.min_ident = min_ident
		self.max_evalue = max_evalue
		self.score_cutoff_fraction = score_cutoff_fraction
		self.keep_max_hit_fraction = keep_max_hit_fraction
		self.keep_min_hit_count = keep_min_hit_count
		if continue_from_json:
			assert len(blastfiles) == 1, "ERROR: if supposed to create blastdata-object from json, you can provide only one input file (the json)"
			self.from_json(blastfiles[0])
		else:
			self.blastlinelist = []
			for bf in blastfiles:
				self.read_blast_tsv(bf, max_evalue = max_evalue, min_ident = min_ident) #TODO: Note: not passing bindata-object right here. if it needs to be looped through later anyway (after condensing the list) it makes more sense to do all further assignments later at that point
			self.blastlinelist = [ dict(t) for t in {tuple(bl.items()) for bl in self.blastlinelist} ] #remove duplicate blast hits that apparently turn up in gtdb and silva dbs ... alternatively i could simply only blast vs silva, but that could miss some potentially incorrectly called rna genes...
			self.blastlinelist = self.sort_blastlines_by_gene()
			self.filter_hits_per_gene(score_cutoff_fraction, keep_max_hit_fraction, keep_min_hit_count)
	
	# ~ def _prot_or_nuc(self, stype):
		# ~ stype_query = stype.strip().split()[0]
		# ~ if stype_query in self._nuc_stypelist:
			# ~ return "nuc"
		# ~ elif stype_query in self._prot_stypelist:
			# ~ return "prot"
		
			
	def filter_hits_per_gene(self, score_cutoff_fraction = 0.75, keep_max_hit_fraction = 0.5, keep_min_hit_count = 2): #todo: allow additional filter settings for rRNA data (e.g. filter by identity not score)
		"""
		assumes the blastlinelist is already sorted by decreasind score!
		"""
		assert 0<= score_cutoff_fraction <= 1, "score_cutoff_fraction must lie in the range 0.0 - 1.0 !"
		assert 0< keep_max_hit_fraction <= 1, "score_cutoff_fraction must be larger than 0.0, and lower or equal to 1.0!"
		previous_query = None
		query_best_score = 0
		#blindex = 0
		query_templist = []
		filtered_blastlinelist = []
		for currline in self.blastlinelist: 
			#steps: store all hits for a specific query in a seperate templist (if above cutoff_fraction). 
			#when the next query is reached check howmany entries are in templist. write obly the top <keep_max_fraction> to new_blastlinelist (but at least <keep_minhit_count>, if possible)
			#currline = self.blastlinelist[blindex]
			if currline["query"] == previous_query:
				#print(" --> {} < ({} * {} ) ? --> {} --> {}!".format(currline["score"], query_best_score, score_cutoff_fraction, (query_best_score * score_cutoff_fraction), currline["score"] < (query_best_score * score_cutoff_fraction)))
				if currline["score"] < (query_best_score * score_cutoff_fraction):
					continue
				query_templist.append(currline)
			else:
				#print(currline["query"])
				if len(query_templist) > 0:
					filtered_blastlinelist.extend(query_templist[:max(int(len(query_templist)*keep_max_hit_fraction), keep_min_hit_count)])
				query_templist = [currline]
				query_best_score = currline["score"]
			previous_query = currline["query"]
		if len(query_templist) > 0: #after finished looping, append the final gene also...
			filtered_blastlinelist.extend(query_templist[:max(int(len(query_templist)*keep_max_hit_fraction), keep_min_hit_count)])
		self.blastlinelist = filtered_blastlinelist
	
	def add_info_to_blastlines(self, bindata_obj, taxdb_obj = None):
		import time #todo: remove this later
		print("add_info_to_blastlines old version")
		starttime = time.time()
		for i in range(len(self.blastlinelist)):
			self.blastlinelist[i]["contig"] = bindata_obj.marker2contig(self.blastlinelist[i]["query"])
			self.blastlinelist[i]["stype"] = bindata_obj.markerdict[self.blastlinelist[i]["query"]]["stype"]
			# ~ if self.seqtype != None:
				# ~ if self._prot_or_nuc(self.blastlinelist[i]["stype"]) != self.seqtype:
					# ~ self.seqtype = "mixed"
			# ~ else:
				# ~ self.seqtype = self._prot_or_nuc(self.blastlinelist[i]["stype"])
			if taxdb_obj != None:
				# ~ print("--{}--".format(self.blastlinelist[i]))
				# ~ import pdb; pdb.set_trace()
				self.blastlinelist[i]["taxid"] = taxdb_obj.acc2taxid(self.blastlinelist[i]["subject"])[0]
		endtime = time.time()
		print("\nthis took {} seconds\n".format(endtime - starttime))

	def from_json(self, jsonfilename):
		self.blastlinelist = misc.from_json(jsonfilename)

	def to_json(self, jsonfilename):
		misc.to_json(self.blastlinelist, jsonfilename)
	# ~ def add_info_to_blastlines(self, bindata_obj, taxdb_obj = None):
		# ~ import time #todo: remove this later
		# ~ print("add_info_to_blastlines NEW version")
		# ~ starttime = time.time()
		# ~ if taxdb_obj != None:
			# ~ acc2taxiddict = taxdb_obj.acclist2taxiddict(list({bl["subject"] for bl in self.blastlinelist})) 
		# ~ for i in range(len(self.blastlinelist)):
			# ~ self.blastlinelist[i]["contig"] = bindata_obj.marker2contig(self.blastlinelist[i]["query"])
			# ~ self.blastlinelist[i]["stype"] = bindata_obj.markerdict[self.blastlinelist[i]["query"]]["stype"]
			# ~ if taxdb_obj != None:
				# ~ self.blastlinelist[i]["taxid"] = acc2taxiddict.get(self.blastlinelist[i]["subject"])
		# ~ endtime = time.time()
		# ~ print("\nthis took {} seconds\n".format(endtime - starttime))
	
	def sort_blastlines_by_gene(self):
		return sorted(self.blastlinelist, key = lambda x: (x["contig"], x["query"], -x["score"])) #sorts hits first increasingly by contig and query-name (not the same in case of rRNA genes), then decreasingly by score	

	def sort_blastlines_by_contig(self):
		return sorted(self.blastlinelist, key = lambda x: (x["contig"], -x["score"])) #sorts hits first increasingly by contig and query-name (not the same in case of rRNA genes), then decreasingly by score	
			
	def get_best_hits_per_gene(self, keep_max_best_hit_fraction = 1.0, keep_min_hit_count = 2): #max_best_hitfraction and keep_min_hit_count are not supposed to be actually used. should already be taken care of by "filter_hits_by_gene". Just keeping option open to filter again using different cutoffs	#todo: combine with below function. add keyword to return contig or gene
		#print("OI!")
		assert 0< keep_max_best_hit_fraction <= 1, "score_cutoff_fraction must be larger than 0.0, and lower or equal to 1.0!"
		previous_query = None
		#blindex = 0
		query_templist = []
		filtered_blastlinelist = []
		for currline in self.sort_blastlines_by_gene():
			if currline["query"] == previous_query:
				query_templist.append(currline)
			else:
				if len(query_templist) > 0:
					yield previous_query, [ hit(accession=q["subject"], taxid=q["taxid"], identity=q["ident"], score=q["score"]) for q in query_templist[:max(int(len(query_templist)*keep_max_best_hit_fraction), keep_min_hit_count)]]
				query_templist = [currline]
				query_best_score = currline["score"]
			previous_query = currline["query"]
		if len(query_templist) > 0:
			yield previous_query, [ hit(accession=q["subject"], taxid=q["taxid"], identity=q["ident"], score=q["score"]) for q in query_templist[:max(int(len(query_templist)*keep_max_best_hit_fraction), keep_min_hit_count)]]
	
	def get_best_hits_per_contig(self, keep_max_best_hit_fraction = 1.0, keep_min_hit_count = 20): #todo: combine with above function. add keyword to return contig or gene
		#hit = namedtuple('hit', 'accession taxid identity score')
		assert 0<  keep_max_best_hit_fraction <= 1, "score_cutoff_fraction must be larger than 0.0, and lower or equal to 1.0!"
		previous_contig = None
		#blindex = 0
		contig_templist = []
		filtered_blastlinelist = []
		for currline in self.sort_blastlines_by_contig: 
			#steps: store all hits for a specific query in a seperate templist (if above cutoff_fraction). 
			#when the next query is reached check howmany entries are in templist. write obly the top <keep_max_fraction> to new_blastlinelist (but at least <keep_minhit_count>, if possible)
			#currline = self.blastlinelist[blindex]
			if currline["contig"] == previous_contig:
				contig_templist.append(currline)
			else:
				if len(query_templist) > 0:
					yield previous_contig, [ hit(accession=q["subject"], taxid=q["taxid"], identity=q["ident"], score=q["score"]) for q in query_templist[:max(int(len(query_templist)* keep_max_best_hit_fraction), keep_min_hit_count)]]
				query_templist = [currline]
				query_best_score = currline["score"]		
			previous_contig = currline["contig"]

	def read_blast_tsv(self, infilename, max_evalue = None, min_ident = None, bindata_obj = None): #todo: test performace/memory-usage difference when using dataclasses or namedtuples instead of dictionaries for blast entries
		"""
		reads in each line from input blast files as dictionary. Stored all lines in a common list "blastlinelist".
		already roughly filters by maximal evalue or minimum identity if correpsonding arguments are set.
		already associated query contig and markertype to each query marker if a bindata-object is passed
		will not assign taxids at this point yet, because this list is likely to be greatly reduced in later steps and assigning taxids is a relatively slow process  
		"""
		def string_or_int_or_float(teststring): #todo: consider moving this to misc.py (if ever needed somewhere else?)
			try:
				testfloat = float(teststring)
			except (TypeError, ValueError):
				return teststring
			else:
				try:
					testint = int(teststring)
				except (TypeError, ValueError):
					return testfloat
				else:
					return testint

		infile = openfile(infilename)
		print("NOWHANDLINGFILE: {}".format(infilename))
		for line in infile:
			tokens = line.strip().split("\t")
			bl = { x : string_or_int_or_float(tokens[self._blasttsv_columnnames[x]]) if type(self._blasttsv_columnnames[x]) == int else None for x in self._blasttsv_columnnames}
			if max_evalue and bl["evalue"] > max_evalue: #bl.evalue > max_evalue:
				continue
			if min_ident and bl["ident"] < min_ident: #bl.ident < min_ident:
				continue
			if bindata_obj != None:
				bl["contig"] = bindata.marker2contig(bl["query"])
				bl["stype"] = bindata.markerdict[bl["query"]]["stype"]
			# ~ print("*"*20)
			# ~ print(bl)
			# ~ print("*"*20)
			self.blastlinelist.append(bl)
		

	def get_contradicting_tophits(self, markernames, db, cutoff, markerlevel): #only checking tax-classification contradictions on phylumlevel or below!
		import lca, getdb
		def get_besthit_info(db, blastline, blastlineindex):
			subject = blastline["subject"]
			taxid = blastline["taxid"]
			domain, phylum = None, None
			taxpath = db.taxid2taxpath(taxid)
			if len(taxpath) >= 2:
				domain = db.taxid2taxpath(taxid)[1][0]
				if len(taxpath) >= 3:
					phylum = db.taxid2taxpath(taxid)[2][0]
			identity = blastline["ident"]
			score = blastline["score"]
			return {"subject": subject, "taxid": taxid, "domain": domain, "phylum": phylum, "identity" : identity, "score": score, "blastlineindex" : blastlineindex}
			
		# ~ tempmarkerdict = { mn : { "tophit": None, "tophit_blastline" : None, "contrahit": None, "contrahit_blastline" : None, "lca_taxlevel" : None} for mn in markernames }
		inconsistencies = { mn: { "lca": None, "best_hit" : None, "contradiction_level" : None, "best_contradiction" : None, "above_cutoff": False} for mn in markernames}   
		tempsinglemarkerlcas = []
		for mn in markernames:
			templines = sorted([ self.blastlinelist[x] for x in range(len(self.blastlinelist)) if self.blastlinelist[x]["query"] == mn ], key =  lambda x: -x["score"])
			tophit = templines[0]
			inconsistencies[mn]["lca"] = lca.strict_lca(db, blasthitlist = _blastlines2blasthits(templines)) #todo: is redundant to repeat this here. Should instead save individuyl marker-lcas when first doing them!
			inconsistencies[mn]["best_hit"] = get_besthit_info(db, tophit, 0) 
			templca = templines[0]["taxid"]
			currentline = 1
			while currentline < len(templines):
				templca = db.get_strict_pairwise_lca(templca, templines[currentline]["taxid"])
				templcataxlevel = db.taxid2taxlevel(templca)
				if templcataxlevel not in lca.taxlevels[3:] + [ "ignored rank"]:  # refdb-inconsistency if lca is "root" (= "no rank" at the moment, unfortunately), "domain" (aka "superkingdom") or "phylum" (=lca.taxlevels[:3]) despite high average identities
					inconsistencies[mn]["best_contradiction"] = get_besthit_info(db, templines[currentline], currentline)
					if templines[currentline]["ident"] >= cutoff:
						inconsistencies[mn]["above_cutoff"] = True
					contradictionlevelindex = list(getdb.rank2index.keys()).index(templcataxlevel) + 1 #can be only "domain", "phylum" or None, but this way allows adding more levels in the future easier
					actual_contradictionlevel = lca.taxlevels[contradictionlevelindex]
					inconsistencies[mn]["contradiction_level"] = actual_contradictionlevel
					break
				currentline += 1
		
		tempsinglemarkerlcas = [ inconsistencies[mn]["lca"] for mn in markernames if inconsistencies[mn]["lca"] != None ]			
		outtaxpath, top2_contras, top2_contras_avidents = lca.weighted_lca(db, blasthitlist=tempsinglemarkerlcas, taxlevel=markerlevel, return_contradicting_top2 = True)
		
		print("\n"+"-"*50)
		print("looking at refdb_inconsisentcies:")
		print("single marker lcas")
		print("+"*50)
		for mn in inconsistencies:
			if inconsistencies[mn]["best_contradiction"] == None:
				print("{} : no contradiction on or below phylum-level".format(mn))
				continue
			print(mn)
			print("  best hit: {}".format(inconsistencies[mn]["best_hit"]))
			print("  best contra: {}".format(inconsistencies[mn]["best_contradiction"]))
			print(" contra_level: {}, above_cutoff : {}".format(inconsistencies[mn]["contradiction_level"], inconsistencies[mn]["above_cutoff"]))
			print("___")
		print("weighted_lca")
		if top2_contras != None:
			print([val for pair in zip(top2_contras, top2_contras_avidents) for val in pair])
		else:
			print("no (further) contradiction on weighted_lca level")
			print(outtaxpath)
		
				#next steps: while resulting taxlevel > phylum: make strict_lca of rising fraction of blastlies. return the taxid of the fist blastline that returns a strict lca of < phylum (with average id >= 80)
	# ~ def pickleyourself(self, filename): #todo: does not work yet! figure out why!
		# ~ try:
			# ~ import cPickle as pickle #this way the faster cPickle gets used IF available, but the slower standard pickle gets used otherwise
		# ~ except:
			# ~ import pickle
		# ~ f = open(filename, 'wb')
		# ~ pickle.dump(self.__dict__, f, 2)
		# ~ f.close()
	
	# ~ def unpickleyourself(self, filename):#todo: does not work yet! figure out why!
		# ~ try:
			# ~ import cPickle as pickle #this way the faster cPickle gets used IF available, but the slower standard pickle gets used otherwise
		# ~ except:
			# ~ import pickle
		# ~ f = open(filename, 'rb')
		# ~ tmp_dict = pickle.load(f)
		# ~ f.close()          
		# ~ self.__dict__.update(tmp_dict)
				
def read_blast_tsv(infilename, max_evalue = None, min_ident = None, dbobj = None, bindata_obj = None):
	""" 
	returns a a list of dictionaries, each representing the data in a blast line
	if max_evalue or min_ident are set, elines above or below these cutoffs are ignored
	each blast will later be assigned to a contig and to a "subject-type" (=stype), which can be either of ["16S", "23S", "univ_marker", "proc_marker", "bact_marker", "arch_marker", "other"]
	"""
	#note to self: wanted to use namedtuples here, but namedtuples won't work well here, because i may need them to be mutable.
	columninfos = {"query" : 0, "subject" : 1, "ident" : 2, "alignlen" : 3, "qstart": 6, "qend" : 7, "sstart" : 8, "ssend" : 9, "evalue" : 10, "score" : 11, "contig" : None, "stype" : None, "taxid": None} #Since python 3.7 all dicts are now ordered by default (YAY!). So the order of keys is maintained, without having to keep a seperate "orderlist". --> THEREFORE:  Ensure/Enforce python 3.7+ !!!
	infile = openfile(infilename)
	blastlinelist = []
	for line in infile:
		tokens = line.strip().split("\t")
		bl = { x : tokens[columninfos[x]] if type(columninfos[x]) == int else None for x in columninfos}
		if bindata_obj != None:
			bl[contig] = bindata.marker2contig(bl["query"])
		if max_evalue and bl["evalue"] > max_evalue: #bl.evalue > max_evalue:
			continue
		if min_ident and bl["ident"] < min_ident: #bl.ident < min_ident:
			continue
		blastlinelist.append(bl)
	infile.close()
	return blastlinelist

def _blastlines2blasthits(blastlinelist):
	return [ hit(accession=q["subject"], taxid=q["taxid"], identity=q["ident"], score=q["score"]) for q in blastlinelist ]

def __add_contigs2blasthits_later(blastlinelist, parsetype = "prodigal", lookup_table = None): #todo: probably obsolete...
	"""
	assigns contigs to blast hits
	parsetype must be one of ["prodigal", "rnammer", "barrnap", "lookup"]
	if parsetype is set to "lookup", a lookup_table must also be provided
	"""
	#start of subfunctions
	def _parse_prodigal(blastlinelist):
		#sys.stderr.write("\nPARSING PRODIGAL STYLE\n")
		import re
		pattern = re.compile("_\d+$") #this pattern only works for standard naming scheme of standard prodigal output
		for blindex in range(len(blastlinelist)):
			assert re.search(pattern, blastlinelist[blindex]["query"]), "\nERROR: protein {} does not seem to be named according to prodigal standards. You should supply a lookup-table\n".format(blastlinelist[blindex]["query"])
			contig = re.sub(pattern, "", blastlinelist[blindex]["query"])
			blastlinelist[blindex]["contig"] = contig
			#sys.stderr.write("\n{} --> {}\n".format(blastlinelist[blindex]["query"], contig))
		return blastlinelist
	
	def _parse_lookup_table(blastlinelist, lookup_table):
		for blindex in range(len(blastlinelist)):
			contig = lookup_table[blastlinelist[blindex]["query"]]
			blastlinelist[blindex]["contig"] = contig
		return blastlinelist
		
	def _parse_rnammer(blastlinelist):
		import re
		patternprefix = re.compile("^\d+S_rRNA_") #this pattern only works for standard naming scheme of standard rnammer output
		patternsuffix = re.compile("_\d+-\d+_DIR[+-]$") #this pattern only works for standard naming scheme of standard rnammer output
		assert re.search(patternprefix, blastlinelist[blindex]["query"]) and re.search(patternsuffix, blastlinelist[blindex]["query"]), "\nERROR: Query {} does not seem to be named according to rnammer standards. You should supply a lookup-table\n".format(blastlinelist[blindex]["query"])
		for blindex in range(len(blastlinelist)):
			contig = re.sub(patternprefix, "", blastlinelist[blindex]["query"])
			contig = re.sub(patternsuffix, "", contig)
			blastlinelist[blindex]["contig"] = contig
		return blastlinelist
	
	def _parse_barrnap(blastlinelist):
		import re
		patternprefix = re.compile("^\d+S_rRNA::") #this pattern only works for standard naming scheme of standard barrnap output
		patternsuffix = re.compile(":\d+-\d+_\([+-]\)$") #this pattern only works for standard naming scheme of standard barrnap output
		assert re.search(patternprefix, blastlinelist[blindex]["query"]) and re.search(patternsuffix, blastlinelist[blindex]["query"]), "\nERROR: Query {} does not seem to be named according to barrnap standards. You should supply a lookup-table\n".format(blastlinelist[blindex]["query"])
		for blindex in range(len(blastlinelist)):
			contig = re.sub(patternprefix, "", blastlinelist[blindex]["query"])
			contig = re.sub(patternsuffix, "", contig)
			blastlinelist[blindex]["contig"] = contig
		return blastlinelist
		
	#end of subfunctions
	assert parsetype in ["prodigal", "rnammer", "barrnap", "lookup"], "Do not recognize parsetype {}".format(parsetype)
	if parsetype == "lookup":
		assert lookup_table != None, "if parsetype is set to \"lookup\", a corresponding lookup_table must also be provided"
		return _parse_lookup_table(blastlinelist, lookup_table)
	if parsetype == "prodigal":
		return _parse_prodigal(blastlinelist)
	if parsetype == "rnammer":
		return _parse_rnammer(blastlinelist)
	if parsetype == "barrnap":
		return _parse_barrnap(blastlinelist)
	
def __add_stype2blasthits_later(blastlinelist, markersetdict): #todo: probably obsolte... #markersetdict should be derived from a different module "fastahandler.py" that reads the markerfastas and returns either the sequences, or just the locus_tags/fasta_headers.
	#"markersetlist" should be a list of sets. sets should be ordered in decreasing hierarchy (16S first, then 23S then markerprots then totalprots)
	#each set contains the locus_tags of the markers used in the blast
	for blindex in range(len(blastlinelist)):
		for ms in markersetdict:
			if blastlinelist[bl]["query"] in markersetdict[ms]:
				blastlinelist[bl]["stype"] = ms
	return blastlinelist

def __add_taxid2blasthits_later(blastlinelist, db_obj): #todo: probably obsolte...
	for blindex in range(len(blastlinelist)):
		blastlinelist[bl]["taxid"] = db_obj.acc2taxid(blastlinelist[bl]["taxid"]["subject"])
	return blastlinelist

def run_single_blast(query, db, blast, outname, threads = 1): 
	import subprocess
	assert os.path.basename(blast) in ["blastn", "blastp", "blastx"]
	if type(query) == str and os.path.isfile(query):
		inputarg = None
	elif type(query) == list:
		inputarg =  "\n".join([record.format("fasta") for record in query])
		query = "-"
	else:
		raise IOError("\nERROR: don't recognize query argument\n")
	blastcmd = subprocess.run([blast, "-query", query, "-db", db, "-evalue", "1e-10",\
							   "-outfmt", "6", "-num_threads", str(int(threads)), "-out", outname + ".tmp"], \
							  input = inputarg, stdout = subprocess.PIPE, stderr = subprocess.PIPE, text = True)
	try:
		blastcmd.check_returncode()
	except Exception:
		sys.stderr.write("\nAn error occured during blastp run with query '{}'\n".format(query))
		sys.stderr.write("{}\n".format(blastcmd.stderr))
		raise RuntimeError
	return outname

def run_single_blastp(query, db, blast, outname, threads = 1):
	assert os.path.basename(blast) == "blastp" 
	return run_single_blast(query, db, blast, outname, threads)

def run_single_blastn(query, db, blast, outname, threads = 1):
	assert os.path.basename(blast) == "blastn"
	return run_single_blast(query, db, blast, outname, threads)

def run_single_blastx(query, db, blast, outname, threads = 1):
	assert os.path.basename(blast) == "blastx"
	return run_single_blast(query, db, blast, outname, threads)
	
def run_single_diamondblastp(query, db, diamond, outname, threads = 1): #TODO: currently not setting "--tmpdir" & "--parallel-tmpdir" here! figure something out if this turns out to be problematic on hpc systems
	#TODO: add a maxmem arguemt that states how much memory can be used. use this to determine optimal blocksize and chunks for more efficient blasting. BLOCKSIZE=INT(MEMORY/6) CHUNKS=4/2/1 IF MEMORY >=12/24/48
	import subprocess
	blastcmd = subprocess.run([diamond, "blastp", "--query", query, "--db", db, "--evalue", "1e-10",\
							   "--outfmt", "6", "--threads", str(int(threads)), "--out", outname + ".tmp"], \
							   stdout = subprocess.PIPE, stderr = subprocess.PIPE, text = True)
	try:
		print(" ".join(blastcmd.args))
		blastcmd.check_returncode()
	except Exception:
		sys.stderr.write("\nAn error occured during diamond blastp run with query '{}'\n".format(query))
		sys.stderr.write("{}\n".format(blastcmd.stderr))
		raise RuntimeError
	return outname
	
def _run_any_blast(query, db, path_appl, outname, threads, force = False):
	#TODO: fix stupid problem that blast (unfortunately) cannot handle gzipped query files. solution read in compressed fastas, pipe records to blast-functions via stdin (TODO: add this function also to diamond (already present in ncbi-blast)
	appl = os.path.basename(path_appl)
	#print("\nblasting --> {}".format(outname))
	if os.path.isfile(outname) and force == False:
		sys.stderr.write("\nWARNING: blast result file '{}' already exists and 'force' not set to True --> skipping this blast\n".format(outname))
		return outname
	assert appl in ["blastp", "blastn", "diamond"], "\nError: unknown aligner '{}'\n".format(appl)
	_command_available(command=path_appl)
	if appl == "blastp":
		outname = run_single_blastp(query, db, path_appl, outname, threads)
	elif appl == "blastn":
		outname = run_single_blastn(query, db, path_appl, outname, threads)
	elif appl == "diamond":
		run_single_diamondblastp(query, db, path_appl, outname, threads)
	os.rename(outname + ".tmp", outname)
	return outname

def _distribute_threads_over_jobs(total_threads, num_jobs): # to distribute N threads over M groups as evenly as possible, put (N/M) +1 in (N mod M) groups, and (N/M) in the rest
	##TODO: test using misc.run_multiple_functions_parallel() for this instead! DELETE THIS IF MISC VERSION WORKS!
	if num_jobs == 0:
		sys.stderr.write("nothing to do...\n")
		return [], total_threads
	if num_jobs < total_threads:
		more_threads_list = [ (total_threads / num_jobs) + 1 ] * (total_threads % num_jobs) #these should go to the lower priority markers, as there will be more of them
		fewer_threads_list = [ (total_threads / num_jobs) ] * (num_jobs - len(more_threads_list))
		return fewer_threads_list + more_threads_list, num_jobs
	return [1] * num_jobs, total_threads

def run_multiple_blasts_parallel(basic_blastarg_list, outbasename, total_threads): #basic_blastarg_list = list of tuples such as [(query1, db1, blast1), (query2, db2, blast2),...])
	#TODO: test using misc.run_multiple_functions_parallel() for this instead! DELETE THIS IF MISC VERSION WORKS!
	if len(basic_blastarg_list) == 0:
		sys.stderr.write("nothing to blast...")
	from multiprocessing import Pool
	thread_args, no_processes = _distribute_threads_over_jobs(total_threads, len(basic_blastarg_list))
	arglist = []
	for i in range(len(basic_blastarg_list)):
		if type(basic_blastarg_list[i][0]) == list:
			print("blast input is a list with {} entries".format(len(basic_blastarg_list[i][0])))
			querystring = "auxblasts"
		else:
			querystring = os.path.basename(basic_blastarg_list[i][0])
		outname = "{prefix}_{counter:03d}{appl}_{query}_vs_{db}.tsv".format(prefix = outbasename, counter = i, \
					appl = os.path.basename(basic_blastarg_list[i][2]), query = querystring, \
					db = os.path.basename(basic_blastarg_list[i][1]))
		arglist.append(tuple(bba for bba in basic_blastarg_list[i]) + (outname, thread_args[i]))
	masterblaster = Pool(processes = no_processes)
	outfile_list = masterblaster.starmap(_run_any_blast, arglist)
	#print("finished blasting all")
	#print(outfile_list)
	masterblaster.close()
	masterblaster.join()
	return outfile_list

def make_diamond_db(infasta, outfilename, diamond = "diamond", threads = 1):
	import subprocess
	_command_available(command="diamond")
	mkdbcmd = [ diamond, "makedb", "--in", infasta, "--db", outfilename, "--threads", str(threads) ]
	mkdbproc = subprocess.run(mkdbcmd, stdout = subprocess.PIPE, stderr = subprocess.PIPE, text = True)
	sys.stderr.write("\n creating diamond-db \"{}\"...".format(outfilename))
	sys.stderr.flush()
	try:
		mkdbproc.check_returncode()
		sys.stderr.write("...finished!\n")
		sys.stderr.flush()
	except Exception:
		sys.stderr.write("\nAn error occured while creating '{}'\n".format(outfilename))
		sys.stderr.write("{}\n".format(mkdbproc.stderr))
		raise RuntimeError
	return outfilename	

def make_blast_db(infasta, outfilename, makeblastdb="makeblastdb", db_type="nucl", threads = 1): #threads  argument is ignored. Is only there to work with multiprocessing function
	import subprocess
	#todo: combine with "make_blast_db_from_gz" below: check if infasta ends with .gz. If yes: pipe from zcat, otherwise run makeblastdb directly
	assert db_type in ["nucl", "prot"],  "dbtype has to be either 'nucl' or 'prot'"
	sys.stderr.write("\n-->creating blastdb {}...\n".format(outfilename))
	sys.stderr.flush()
	mkdbcmd = [ makeblastdb, "-in", infasta, "-title", outfilename, "-out", outfilename, "-dbtype", db_type, "-max_file_sz", "4GB"] #hash_index apparently not possible for very large datasets
	print(" ".join(mkdbcmd))
	mkdbproc = subprocess.run(mkdbcmd, stdout = subprocess.PIPE, stderr = subprocess.PIPE, text = True)
	sys.stderr.write("\n creating blastdb-db \"{}\"...".format(outfilename))
	sys.stderr.flush()
	try:
		mkdbproc.check_returncode()
		sys.stderr.write("...finished!\n")
		sys.stderr.flush()
	except Exception:
		sys.stderr.write("\nAn error occured while creating '{}'\n".format(outfilename))
		sys.stderr.write("{}\n".format(mkdbproc.stderr))
		raise RuntimeError
	return outfilename	

def make_blast_db_from_gz(infasta, outfilename, makeblastdb="makeblastdb", db_type="nucl", threads = 1): #threads  argument is ignored. Is only there to work with multiprocessing function
	import subprocess
	_command_available(command="makeblastdb")
	# todo: if this works,incorporate it with "make_blast_db" above! check if infasta ends with .gz. If yes: pipe from zcat, otherwise run makeblastdb directly
	sys.stderr.write("\ncreating blastdb {}...\n".format(outfilename))
	sys.stderr.flush()
	if not infasta.endswith(".gz"):
		sys.stderr.write("\ninput is not compressed! using standard make_blast_db function".format(outfilename))
		return make_blast_db(infasta, outfilename, makeblastdb, db_type, threads) #todo: make this the other way round. make it call this function if it finds input to be compressed
	assert db_type in ["nucl", "prot"],  "dbtype has to be either 'nucl' or 'prot'"
	zcatcmd = [ "zcat", infasta]
	mkdbcmd = [ makeblastdb, "-title", outfilename, "-out", outfilename, "-dbtype", db_type, "-max_file_sz", "4GB"] #hash_index apparently not possible for very large datasets
	import time
	start = time.time()
	zcatproc = subprocess.Popen(zcatcmd, stdout = subprocess.PIPE, stderr = subprocess.PIPE, text = True)
	print(" ".join(mkdbcmd))
	sys.stdout.flush()
	mkdbproc = subprocess.Popen(mkdbcmd, stdin = zcatproc.stdout, stdout = subprocess.PIPE, stderr = subprocess.PIPE, text = True)
	stdout, stderr = mkdbproc.communicate()
	mkdbproc.wait()
	end = time.time()
	print("piping it this way took : {}".format(end - start))
	sys.stdout.flush()
	return outfilename
	#todo: compare how long it take to pipe this via python subprocess and ow long to do the same on the bash shell


def get_blast_combinations(dblist, querylist, blast = "blastn"):
	import itertools #todo move up
	acceptable_blasttools = ["blastn", "blastp", "blastx"]
	assert os.path.basename(blast) in acceptable_blasttools, "Error: 'blast' must be one of {}".format(acceptable_blasttools) #maybe remove this check or add a workaround to use diamond whith this also?
	return [ blasttuple + (blast,) for blasttuple in list(itertools.chain(*list(zip(querylist, permu) for permu in itertools.permutations(dblist, len(querylist)))))] #todo:Only works as long as dblist is longer or equal to querylist... find a way to ensure/workaround this

def _command_available(command="diamond"): #todo: move this to misc, maybe?
	'''
	for checking if diamond or blast are actually installed in the current environment
	'''
	import subprocess
	try:
		subprocess.call([command],stderr=subprocess.PIPE, stdout=subprocess.PIPE)
	except FileNotFoundError:
		raise FileNotFoundError("\nError: can't run command '{}'. Maybe it is not installen in PATH, or you do not have the right conda environment activated?\n")
	else:
		return True
	return False
	
###### Testing functions below:
def test_prodigal_blasts():
	blastfile = sys.argv[1]
	sys.stderr.write("\nreading blast file \"{}\n".format(blastfile))
	blastlinelist = read_blast_tsv(blastfile)
	for bl in blastlinelist[:5]:
		#sys.stderr.write("\n{}\n".format(str(bl)))
		#sys.stderr.write("\n{}\n".format(type(bl)))
		sys.stderr.write("{}\n".format("\t".join([str(bl[x]) for x in bl])))
		sys.stderr.flush()
	sys.stderr.write("\n-------------------------\n")
	sys.stderr.write("\nadding contig_info:\n")
	blastlinelist = contigs2blasthits(blastlinelist)
	for bl in blastlinelist[:5]:
		#sys.stderr.write("\n{}\n".format(str(bl)))
		#sys.stderr.write("\n{}\n".format(type(bl)))
		sys.stderr.write("{}\n".format("\t".join([str(bl[x]) for x in bl])))
		sys.stderr.flush()
	sys.stderr.write("DONE (for now...")

def test_multiblastjobs():
	print("test_multiblastjobs")
	query_list = ["query1.fasta", "query2.fasta", "query3.fasta"]
	db = "/opt/db/blast/nr"
	tool = "blastp"
	jobtuplelist = []
	for query in query_list:
		jobtuplelist.append(("blasthandler", "_run_any_blast", {"query":query, "db":db, "path_appl": tool, "outname":"{}_{}_vs_{}.tab.tmp".format(tool, os.path.basename(query), os.path.basename(db))}))
	run_multiple_functions_parallel(jobtuplelist, 8)
	
def main():
	""" for testing purposes only"""
	#test_prodigal_blasts()
	print("main")
	blastdatabase = "/home/ww5070/delme3/dbdir/gtdb/concat_refgenomes"
	import pdb; pdb.set_trace()
	# ~ test_multiblastjobs()
	print("FINISHED")
####### End of test functions

if __name__ == '__main__':
	main()



