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
from mdmcleaner import misc
hit = namedtuple('hit', 'accession taxid identity score')
supported_outfmts = ["6", "6 std", "6 std qlen", "6 std qlen slen"] 

from mdmcleaner.misc import openfile, run_multiple_functions_parallel

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
			
class blastdata_baseobject(object): #todo: define differently for protein or nucleotide blasts (need different score/identity cutoffs)
	_blasttsv_columnnames = {"query" : 0, "subject" : 1, "ident" : 2, "alignlen" : 3, "qstart": 6, "qend" : 7, "sstart" : 8, "send" : 9, "evalue" : 10, "score" : 11, "qlen" : 12, "slen": 13, "contig" : None, "stype" : None, "taxid": None} #reading in all fields, in case functionality is added later that also uses the sstat send etc fields
	# ~ _nuc_stypelist = ["ssu_rRNA", "lsu_rRNA"]
	# ~ _prot_stypelist = ["total", "bac_marker", "prok_marker", "arc_marker"]
	def __init__(self, *args, max_evalue = None, min_ident = None, score_cutoff_fraction = 0.75, keep_max_hit_fraction = 0.5, keep_min_hit_count = 2,  auxilliary = False, seqtype=None, blacklist=None):
		#methods:
		#	method 1: filter by query-genes. For each gene remove all hits with score below <score_cutoff_fraction> (default = 0.5) of maximum score for that gene
		#			  from these keep the <keep_max_hit_fraction> of hits (default = 0.5), but keep at least <keep_min_fraction> (default = 2) in every case if possible
		#			  later classify each gene by simple lca
		# note: for 16S/23S the score should be more than 1000 and the identity should be above 90??
		# ~ if from_pickle and len(blastfiles) == 1: #if from_pickle is set, that means input shuld NOT be a list of blastfiles, but a single blastdata-object-pickle-file
			# ~ print("blastdata object already exists! loading from pickle")
			# ~ self.unpickleyourself(blastfiles[0])
		assert seqtype in ["nuc", "prot", None], "\nERROR: seqtype must be either 'nuc', 'prot' or None (if unknown)\n"
		assert blacklist == None or type(blacklist) == set, "\nERROR: blacklist must be of type == set\n"
		self.seqtype = seqtype #todo: do something with this (e.g. set some cutoff)
		self.min_ident = min_ident
		self.max_evalue = max_evalue
		self.score_cutoff_fraction = score_cutoff_fraction
		self.keep_max_hit_fraction = keep_max_hit_fraction
		self.keep_min_hit_count = keep_min_hit_count
		if blacklist == None:
			self.blacklist = set()
		else:
			self.blacklist = blacklist
	
	# ~ def _prot_or_nuc(self, stype):
		# ~ stype_query = stype.strip().split()[0]
		# ~ if stype_query in self._nuc_stypelist:
			# ~ return "nuc"
		# ~ elif stype_query in self._prot_stypelist:
			# ~ return "prot"
				 
	
	def filter_hits_per_gene(self, score_cutoff_fraction = 0.75, keep_max_hit_fraction = 0.5, keep_min_hit_count = 2): #todo: allow additional filter settings for rRNA data (e.g. filter by identity not score)
		"""
		assumes the blastlinelist is already sorted by decreasing score!
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

	def filter_blasthits_by_cov_and_ident(self, *, mincov=50.0, minident=90.0, filterbylen = "min", combine_hsps = True):
		'''
		meant for filtering potential refDB contaminations. 
		'mincov' and 'minident' should be given as percentages between 0-100.
		alignment-length filtering can be done based on "query" or "subject" length or the "min" or "max" of those values 
 
		'''
		def get_comparelen(line, filterbylen):
			if filterbylen == "min":
				return min(line["qlen"], line["slen"])
			elif filterbylen == "max":
				return max(line["qlen"], line["slen"])
			else:
				return line[filterbylen]
		
		filterbylen_options = ["min","max","query","subject"]
		assert filterbylen in filterbylen_options, "\n\nERROR: 'filterbylen' must be one of {}\n\n".format(filterbylen_options)

		# ~ import pdb; pdb.set_trace()

		if combine_hsps:
			outlines = []
			query_subject_dict = {}
			for line in self.blastlinelist:
				lkey = (line["query"],line["subject"])
				if lkey in query_subject_dict:
					query_subject_dict[lkey].append(line)
				else:
					query_subject_dict[lkey] = [line]
			for lkey in query_subject_dict.keys():
				sorted_line_sublist = sorted(query_subject_dict[lkey], key = lambda x : (min([x["qstart"],x["qend"]]), -max([x["qstart"], x["qend"]])))
				while len(sorted_line_sublist) > 1:
					a = sorted_line_sublist[0]
					b = sorted_line_sublist[1]
					if max([a["qstart"], a["qend"]]) > min([b["qstart"], b["qend"]]):
						if a["score"] < b["score"]:
							sorted_line_sublist.pop(0)
						else:
							sorted_line_sublist.pop(1)
					else:
						a["ident"] = sum([a["ident"]*a["alignlen"],b["ident"]*b["alignlen"]])/sum([a["alignlen"],b["alignlen"]])
						a["alignlen"] += b["alignlen"]
						a["qend"] = b["qend"]
						a["sstart"] = None #not bothering about alignment orientatione here...(also this could serve to mark this as a "virtual HSP")
						a["send"] = None #not bothering about alignment orientatione here...(also this could serve to mark this as a "virtual HSP")
						a["evalue"] = sum([a["evalue"],b["evalue"]])/2
						a["score"] += b["score"]
						sorted_line_sublist[0] = a
						sorted_line_sublist.pop(1)				
				outlines += sorted_line_sublist
				# ~ import pdb; pdb.set_trace()
			self.blastlinelist = sorted(outlines, key = lambda x:x["score"], reverse=True)
		# ~ import pdb; pdb.set_trace()
		filteredlist = [ line for line in self.blastlinelist if ((line["alignlen"]/get_comparelen(line, filterbylen))*100 >= mincov and line["ident"] >= minident) ] 
		# the following was dropped. delete it when absolutely sure not implementing it again
		# ~ if filterbyoverlap == True: #if the contigs don't overlap completely, make sure the aligning part is not just some conserbed region in the middle of the contigs, but positionoed so that it at least looks like a possible overlap
			# ~ i = 0
			# ~ while i < (len(filteredlist)):
				# ~ if filteredlist[i]["alignlen"] < min(filteredlist[i]["qlen"], filteredlist[i]["slen"]) * 0.95: #allowing difference for up to 5% of contig
					# ~ qmin = min(filteredlist[i]["qstart"],filteredlist[i]["qend"] )
					# ~ qmax = max(filteredlist[i]["qstart"],filteredlist[i]["qend"])
					# ~ smin = min(filteredlist[i]["sstart"],filteredlist[i]["send"])
					# ~ smax = max(filteredlist[i]["sstart"],filteredlist[i]["send"])
					
					# ~ if (qmin > end_mismatch_allowance and qmax < filteredlist[i]["qlen"] - end_mismatch_allowance) or (smin > end_mismatch_allowance and smax < filteredlist[i]["slen"] - end_mismatch_allowance):
						# ~ filteredlist.pop(i)
						# ~ continue
				# ~ i += 1
		self.blastlinelist = filteredlist
		# ~ import pdb; pdb.set_trace()
	
	def add_info_to_blastlines(self, bindata_obj = None, taxdb_obj = None, verbose=True): #todo: add verbose and/or quiet argument
		# ~ import time #todo: remove this later
		if verbose:
			sys.stderr.write("\tadding info to blastlines (old version)\n")
		# ~ starttime = time.time()
		for i in range(len(self.blastlinelist)):
			if bindata_obj != None:
				# ~ import pdb; pdb.set_trace()
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
		# ~ endtime = time.time()
		# ~ print("\nthis took {} seconds\n".format(endtime - starttime))

	def from_json(self, jsonfilename):
		self.blastlinelist = misc.from_json(jsonfilename)

	def to_json(self, jsonfilename):
		misc.to_json(self.blastlinelist, jsonfilename)
		
	# ~ def add_info_to_blastlines(self, bindata_obj, taxdb_obj = None):
		# ~ import time #todo: remove this later
		# ~ sys.stderr.write("\tadding info to blastlines (NEW version)")
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
	
	def sort_blastlines_by_gene(self, contig=None):
		'''
		returns a list of blastlines, sorted first increasingly by contig-id and query-id, and then decreasingly by score.
		if contig is != None: only blastlines corresponding to the specified contig are processed. Otherwise all blastlines are processed.
		'''
		if contig != None:
			return sorted(self.get_blastlines_for_contig(contig), key = lambda x: (x["contig"], x["query"], -x["score"]))
		return sorted(self.blastlinelist, key = lambda x: (x["contig"], x["query"], -x["score"])) #sorts hits first increasingly by contig and query-name (not the same in case of rRNA genes), then decreasingly by score	

	def sort_blastlines_by_contig(self):
		'''
		returns a list of blastlines, sorted first increasingly by contig-id, and then decreasingly by score.
		'''
		
		return sorted(self.blastlinelist, key = lambda x: (x["contig"], -x["score"])) #sorts hits first increasingly by contig-name, then decreasingly by score	
			
	def get_best_hits_per_gene(self, keep_max_best_hit_fraction = 1.0, keep_min_hit_count = 2, contig=None): #max_best_hitfraction and keep_min_hit_count are not supposed to be actually used. should already be taken care of by "filter_hits_by_gene". Just keeping option open to filter again using different cutoffs	#todo: combine with below function. add keyword to return contig or gene
		'''
		A generator function that returns gene-identifiers and corresponding blast-hit-tuples for applying strict lca-classification
		"keep_max_best_hit_fraction" and "keep_min_it_count" do not usially need to be set during standard mdmcleaner runs, as these filtering steps should have been applied during earlier analysis steps
		if set, keep-max-hit-fraction should be a floating-pint value between 0-1. 
		optionally a contig name can be specified. only gene-identifiers and corresponding hit-tuples for that specific contig will be returned
		contig is set to 'None' (default setting), then all contigs are processed.
		'''
		#print("OI!")
		assert 0< keep_max_best_hit_fraction <= 1, "score_cutoff_fraction must be larger than 0.0, and lower or equal to 1.0!"
		previous_query = None
		#blindex = 0
		query_templist = []		
		for currline in self.sort_blastlines_by_gene(contig):
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
			if "_" in teststring:
				return(teststring) # quickfix for unforseen problems caused by PEP 515
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
		# ~ print("NOWHANDLINGFILE: {}".format(infilename))
		for line in infile:
			# ~ import pdb; pdb.set_trace()
			tokens = line.strip().split("\t")
			bl = { x : string_or_int_or_float(tokens[self._blasttsv_columnnames[x]]) if (type(self._blasttsv_columnnames[x]) == int and self._blasttsv_columnnames[x] < len(tokens)) else None for x in self._blasttsv_columnnames}
			# ~ print(bl["contig"])
			# ~ if bl["contig"] == "contam_NZ_JAHGVE010000025.1_15" or bl["subject"] in ["GCF_004341205.1_NZ_SLUL01000011.1", "GCA_002434245.1_DJJR01000029.1"]:
				# ~ import pdb; pdb.set_trace()
			# ~ import pdb; pdb.set_trace()
			# ~ bl["query"] = str(bl["query"]) # quickfix for problems with contigs named as numbers# todo: find a more elegant fix!
			# ~ bl["subject"] = str(bl["subject"]) # quickfix for problems with contigs named as numbers# todo: find a more elegant fix!
			if max_evalue and bl["evalue"] > max_evalue: #bl.evalue > max_evalue:
				continue
			if min_ident and bl["ident"] < min_ident: #bl.ident < min_ident:
				continue
			if bl["subject"] in self.blacklist:
				# ~ print("' {}' is on blacklist --> ignoring!".format(bl["subject"]))
				continue
			if bindata_obj != None:
				bl["contig"] = bindata_obj.marker2contig(bl["query"])
				bl["stype"] = bindata_obj.markerdict[bl["query"]]["stype"]
			self.blastlinelist.append(bl)

	def filter_blacklist(self, blacklist):
		# ~ sys.stderr.write("\n\t (Re-)filtering blasthits for new blacklist entries\n")
		i=0
		while i < len(self.blastlinelist):
			if self.blastlinelist[i]["subject"] in blacklist:
				self.blastlinelist.pop(i)
				continue
			i += 1
		

	def get_blastlines_for_query(self, queryname):
		'''
		returns all blastlines with a certain value in "query"
		query must be a single string representing a query-name or -accession-number
		'''
		return [ self.blastlinelist[x] for x in range(len(self.blastlinelist)) if self.blastlinelist[x]["query"] == queryname ]

	def pop_blastlines_for_query(self, queryname):
		'''
		returns all blastlines with a certain value in "query" and simultaineously removes them from this blastdata-object
		query must be a single string representing a query-name or -accession-number
		'''
		sublist = []
		index = 0
		while index < len(self.blastlinelist):
			if self.blastlinelist[index]["query"] == queryname:
				sublist.append(self.blastlinelist.pop(index))
			else:
				index += 1
		return sublist
	
	def get_blastlines_for_contig(self, contigname):
		'''
		returns all blastlines with a certain value in "contig"
		query must be a single string representing a contig-name or -accession-number
		'''
		return [ self.blastlinelist[x] for x in range(len(self.blastlinelist)) if self.blastlinelist[x]["contig"] == contigname ]		

	def get_contradicting_tophits(self, markernames, db, cutoffs, markerlevel): #only checking tax-classification contradictions on phylumlevel or below!
		#todo: distinguish between:
		#				- probable actual contaminations (in either refdb or sagmag) DONE
		#				- only weakly supported contradictions (should only happen in some cases where contradicitons scored JUST high enough to factor in the weighted LCA) DONE
		#				- likely known gtdb taxon-problems (Firmicutes_A, Firmicuted_B etc, also some stuff in archaea #scrapped that.
		#				- contradicting taxonomies silva/gtdb (this would require recognizing whether an individual entry comes from silfǘa or gtdb)
		#				- rRNA-based silva taxon not backed by genomic data in gtdb
		# 			- keep a lookput for more...
		
		from mdmcleaner import lca, getdb
		cutoff_taxlevels =  ["species", "genus", "order"]
		
		def get_singlehit_domain_phylum_counts(db, markernames, singlemarkerlcadict): #todo: this could be interesting also for general contigdict info? consider putting this in reporting/mdmcleaner...
			from collections import Counter
			lcataxids = [ singlemarkerlcadict[mn]["lca"].taxid for mn in singlemarkerlcadict ]
			domain_phylum_list = [ (db.get_specific_taxlevel_subtaxid(taxid, taxlevel="domain"), db.get_specific_taxlevel_subtaxid(taxid, taxlevel="phylum")) for taxid in lcataxids ]
			domain_phylum_counts = Counter(domain_phylum_list)
			domain_phylum_counts["no_blast_hit"] = len(markernames) - len(domain_phylum_list)
			return domain_phylum_counts
			
			
		def get_besthit_info(db, blastline, blastlineindex):
			subject = blastline["subject"]
			taxid = blastline["taxid"]
			# ~ domain, phylum = None, None
			taxpath = db.taxid2taxpath(taxid)
			domain = db.get_specific_taxlevel_subtaxid(taxid, "domain") #will be None, if taxpath does not include domain
			phylum =  db.get_specific_taxlevel_subtaxid(taxid, "phylum") #will be None, if taxpath does not include phylum
			identity = blastline["ident"]
			score = blastline["score"]
			return {"subject": subject, "taxid": taxid, "domain": domain, "phylum": phylum, "identity" : identity, "score": score, "blastlineindex" : blastlineindex}

		def categorize(outinfo, singlemarkerlcadict):
			#if there is an sm_contradiction and/or weighted_lca_contradiction: --> potential refdb-contamination
			#		--> if contra-identity > species-cutoff --> HIGH indicator (filter out)
			#		--> if species-cutoff > contra-identity > genus-cutoff --> moderate indicator (filter out)
			#		--> if genus-cutoff > contraidentity  --> low indicator (filter out only, if best hit contradicts consenus_lcs)
			#		otherwise: fringe case for current cutoff settings (most likely not a contamination after all, but taxon should be checked nethertheless. Filter out only if best hit contradicts consensus-LCA)
			# if the marker is rRNA, and the blasthits include gtdb-annoations above phylum with high identity BUT equally good silva hits say "None" --> silva/gtdb-conflicting taxomony
			#oterwise, if there is no contradiction, but the lca is still None: --> silva taxon with no representation in GTDB
			amb_type = "" #can be potential refDB-contamination("high indication", "moderate indication", "fringe-case")on sm- and/or weighted-LCA level, silva-gtdb-ambiguity(-->should choose lower marker-level if possible) 
			amb_evidence = ""
			amb_infotext = ""
			if outinfo["sm_representative_contradiction"] != None:
				amb_infotext = "contradictions within single marker LCAs"	
				if outinfo["sm_representative_contradiction"]["bestcontraident"] >= outinfo["cutoffs"][0]:
					amb_type = "potential refDB-contamination [high indication sm-LCA level]"
				elif outinfo["sm_representative_contradiction"]["bestcontraident"] >= outinfo["cutoffs"][1]:
					amb_type = "potential refDB-contamination [moderate indication sm-LCA level]"
				elif outinfo["sm_representative_contradiction"]["bestcontraident"] >= outinfo["cutoffs"][2]:
					amb_type = "potential refDB-contamination [low indication sm-LCA level]"
				else:
					amb_type = "fringe case [sm-LCA level]" #todo: in these cases the script should check whether the best hit agrees with consenus-LCA. Yes --> keep contig, No --> put contig in "potential-refdb-contaminations"
					amb_infotext +=  "; sm-LCA affected by few low-identity cross-phylum/domain hits"
				
				bh_evidence = "'{bhtaxid}'({bhdomphyl}; acc='{bhacc}'; ident={bhident:.2f}%)".format(bhtaxid=outinfo["sm_representative_contradiction"]["besthit_taxid"],\
				bhdomphyl=outinfo["sm_representative_contradiction"]["besthitdomphyl"],bhacc=outinfo["sm_representative_contradiction"]["besthit_subject"], bhident=outinfo["sm_representative_contradiction"]["besthitident"])
				bc_evidence = "'{bctaxid}'({bcdomphyl};acc='{bcacc}'; ident={bcident:.2f}%)".format(bctaxid=outinfo["sm_representative_contradiction"]["bestcontra_taxid"],bcdomphyl=outinfo["sm_representative_contradiction"]["bestcontradomphyl"],bcacc=outinfo["sm_representative_contradiction"]["bestcontra_subject"],\
				bcident=outinfo["sm_representative_contradiction"]["bestcontraident"])
				if outinfo["markerlevel"] in ["lsu_rRNA_tax", "ssu_rRNA_tax", "tsu_rRNA_tax"]:
					bh_accsource = db._gtdb_refseq_or_silva(outinfo["sm_representative_contradiction"]["besthit_subject"])
					bc_accsource = db._gtdb_refseq_or_silva(outinfo["sm_representative_contradiction"]["bestcontra_subject"])
					if bh_accsource != bc_accsource:
						amb_infotext += " and/or contradicting SILVA vs gtdb taxonomies"
						bh_evidence = "{}; {}".format(bh_accsource, bh_evidence)
						bc_evidence = "{}; {}".format(bc_accsource, bc_evidence)
				amb_evidence += " sm_best hit={} ;; sm_best contradiction={}".format(bh_evidence, bc_evidence)

			if outinfo["weighted_lca_top_contradictions"] != None: #todo: I know, I know! Is redundant! I don't have the time to optimize RN. Something for version 1.1...
				amb_infotext+="; contradictions between individual LCAs (weighted LCA)"
				if outinfo["weighted_lca_top_contradictions"]["bestcontraident"] >= outinfo["cutoffs"][0]:
					amb_type += "potential refDB-contamination or chimeric-contig [high indication weighted-LCA level]"
				elif outinfo["weighted_lca_top_contradictions"]["bestcontraident"] >= outinfo["cutoffs"][1]:
					amb_type += "potential refDB-contamination or chimeric-contig [moderate indication weighted-LCA level]"
				elif outinfo["weighted_lca_top_contradictions"]["bestcontraident"] >= outinfo["cutoffs"][2]:
					amb_type += "potential refDB-contamination or chimeric-contig [low indication weighted-LCA level]"	
				else:
					amb_type += "fringe case [weighted-LCA level]" #todo: in these cases the script should check whether the best hit agrees with consenus-LCA. Yes --> keep contig, No --> put contig in "potential-refdb-contaminations"
					amb_infotext +=  "; weighted LCA affected by few low-identity cross-phylum/domain hits"		
				amb_evidence += " weighted_best_tax={}(identity={:.2f}%, score={});;weighted_best_contradiction={}(identity={:.2f}%, score={})".format(outinfo["weighted_lca_top_contradictions"]["besthit_taxid"], outinfo["weighted_lca_top_contradictions"]["besthitident"], outinfo["weighted_lca_top_contradictions"]["besthitscore"],outinfo["weighted_lca_top_contradictions"]["bestcontra_taxid"], outinfo["weighted_lca_top_contradictions"]["bestcontraident"], outinfo["weighted_lca_top_contradictions"]["bestcontrascore"])
				
			if outinfo["sm_representative_contradiction"] == outinfo["weighted_lca_top_contradictions"] == None:
				# ~ hitlines = singlemarkerlcadict[list(singlemarkerlcadict.keys())[0]]["blastlines"] #previous version: looked only blast hits of first marker
				hitlines = []
				for sm in singlemarkerlcadict.keys(): #now: in case tehre are multiple markers (e.g. proteins) get blast hits of all of them and choose the three best based on score (not just simply take the three first hits of the FIRST marker)
					hitlines += singlemarkerlcadict[mn]["blastlines"]
				hitlines = sorted(hitlines, key = lambda x: -x["score"])
				best3hitlines =  hitlines[:3]
				if db.is_eukaryote(best3hitlines[0]["taxid"]):
					amb_infotext = "probable eukaryotic contamination"
					amb_type = False
				elif outinfo["markerlevel"] in ["ssu_rRNA_tax", "lsu_rRNA_tax", "tsu_rRNA_tax"]:
					best_gtdb_hit = None
					best_gtdb_hits = [bl for bl in hitlines if db._gtdb_refseq_or_silva(bl["subject"]) == "gtdb"]
					if len(best_gtdb_hits) > 0:
						best_gtdb_hit = best_gtdb_hits[0]
					if db._gtdb_refseq_or_silva(best3hitlines[0]["subject"]) == "silva" and  best3hitlines[0]["ident"] >= lca.species_identity_cutoffs[outinfo["markerlevel"]] and (best_gtdb_hit == None or (db._gtdb_refseq_or_silva(best_gtdb_hit["subject"]) and best_gtdb_hit["ident"] < lca.species_identity_cutoffs[outinfo["markerlevel"]])):
						amb_type = "unrepresented silva taxon/OTU"
						amb_infotext = "matches a silva taxon/OTU on '{}'-level, which is apparently not represented with any genome data in gtdb".format(outinfo["markerlevel"])
					else:
						amb_type = "gtdb/silva database ambiguity" #todo: in such cases downstream check protein_based markers (marker_prots or totalprots) also
						amb_infotext = "gtdb/silva database ambiguity"
						if best_gtdb_hit:
							amb_infotext += "(possibly conflicting taxonomic placements in silva compared to gtdb)"
							amb_evidence += "best gtdb-hit: acc={},taxid={},ident={:.2f}; ".format(best_gtdb_hit["subject"], best_gtdb_hit["taxid"], best_gtdb_hit["ident"])
						else:
							amb_infotext += "(possibly a novel taxon with hits to silva-OTUs but not to gtdb-representatives)"
				else:
					amb_type = "wtf"
					amb_infotext = "wtf" #todo: look for such cases and try to figure them out if they occur!
				amb_evidence += "overall best three blast hits: {}".format("; ".join(["db={},acc={},inferred_gtdb_taxid={},ident={:.2f}".format(db._gtdb_refseq_or_silva(h["subject"]), h["subject"], h["taxid"], h["ident"]) for h in best3hitlines]))
			#todo: maybe add a key "checktaxa" to outinfo, containing a list of taxa (sm- and weighted LCA best hits) to check against the consensus classification	
			return amb_type, amb_evidence, amb_infotext
				
		outinfo = {"markerlevel" : markerlevel, "cutoffs" : cutoffs, "single_or_weighted_discrep": None, "discrepancy_taxlevel" : None, "amb_type" : None, "amb_infotext" : None, "amb_evidence": None, "sm_representative_contradiction" : None, "sm_contradiction_counts" : None,  "sm_domain_phylum_classification_counts" : None, "weighted_lca_top_contradictions": None}
		singlemarkerlcadict = { mn: { "lca": None, "blastlines" : sorted(self.get_blastlines_for_query(mn), key = lambda x: -x["score"]), "info": {"best_hit" : None, "contradiction_level" : None, "best_contradiction" : None, "above_cutoff": None} }for mn in markernames if len(self.get_blastlines_for_query(mn)) > 0}   
		tempsinglemarkerlcas = []
		lowestcli = 999
		for mn in singlemarkerlcadict:
			# ~ print(mn)
			templines = singlemarkerlcadict[mn]["blastlines"]
			# ~ import pdb; pdb.set_trace()
			tophit = templines[0]
			# ~ breakme = False #todo: debugging only
			# ~ if tophit["contig"] == "PHBV01000266.1": #todo: debugging only
				# ~ breakme = True #todo: debugging only
			singlemarkerlcadict[mn]["lca"] = lca.strict_lca(db, blasthitlist = _blastlines2blasthits(templines)) #todo: is redundant to repeat this here. Should instead save individuyl marker-lcas when first doing them!
			singlemarkerlcadict[mn]["info"]["best_hit"] = get_besthit_info(db, tophit, 0) 
			templca = templines[0]["taxid"]
			currentline = 1
			while currentline < len(templines):
				templca = db.get_strict_pairwise_lca(templca, templines[currentline]["taxid"])
				templcataxlevel = db.taxid2taxlevel(templca)
				if templcataxlevel not in lca.taxlevels[3:] + [ "ignored rank"]:  # refdb-inconsistency if lca is "root" (= "no rank" at the moment, unfortunately), "domain" (aka "superkingdom") or "phylum" (=lca.taxlevels[:3]) despite high average identities
					hitinfo = get_besthit_info(db, templines[currentline], currentline)
					if hitinfo["domain"] == singlemarkerlcadict[mn]["info"]["best_hit"]["domain"] and  hitinfo["phylum"] == singlemarkerlcadict[mn]["info"]["best_hit"]["phylum"]: #e.g. if BOTH hits say "domain=Archaeum; phylum=None", they do not actually contradict each other...
						currentline += 1
						continue
					singlemarkerlcadict[mn]["info"]["best_contradiction"] = hitinfo
					for tli in range(len(cutoff_taxlevels)):
						if templines[currentline]["ident"] >= cutoffs[tli]:
							singlemarkerlcadict[mn]["info"]["above_cutoff"] = cutoff_taxlevels[tli] #="species" if above species-cutoff (strong), "genus" if above genus cutoff(moderate), "order" if above order cutoff (weak), None if below order cutoff (fringe-case)
							break
					contradictionlevelindex = list(getdb.rank2index.keys()).index(templcataxlevel) + 1 #can be only "domain", "phylum" or None, but this way allows relatively easy adding more levels in the future
					actual_contradictionlevel = lca.taxlevels[contradictionlevelindex]
					singlemarkerlcadict[mn]["info"]["contradiction_level"] = actual_contradictionlevel
					if contradictionlevelindex < lowestcli:
						lowestcli = contradictionlevelindex
						outinfo["discrepancy_taxlevel"] = actual_contradictionlevel
					# ~ import pdb; pdb.set_trace()
					break
				currentline += 1
		
		tempsinglemarkerlcas = [ singlemarkerlcadict[mn]["lca"] for mn in singlemarkerlcadict if singlemarkerlcadict[mn]["lca"] != None ]
		sm_contradiction_overview = {} #todo: not used much yet, but already coded in case further analyses should be added in the near future

		for i in singlemarkerlcadict:
			tempdict = {}
			if singlemarkerlcadict[i]["info"]["best_contradiction"] == None:
				continue
			
			outinfo["single_or_weighted_discrep"] = "single"			
			tempdict["besthit_taxid"] = singlemarkerlcadict[i]["info"]["best_hit"]["taxid"]
			tempdict["besthit_subject"] = singlemarkerlcadict[i]["info"]["best_hit"]["subject"]
			tempdict["besthitdomphyl"] = "{},{}".format(singlemarkerlcadict[i]["info"]["best_hit"]["domain"],singlemarkerlcadict[i]["info"]["best_hit"]["phylum"])
			tempdict["besthitident"] = singlemarkerlcadict[i]["info"]["best_hit"]["identity"]
			tempdict["besthitscore"] = singlemarkerlcadict[i]["info"]["best_hit"]["score"]

			tempdict["bestcontra_taxid"] = singlemarkerlcadict[i]["info"]["best_contradiction"]["taxid"]
			tempdict["bestcontra_subject"] = singlemarkerlcadict[i]["info"]["best_contradiction"]["subject"]			
			tempdict["bestcontradomphyl"] = "{},{}".format(singlemarkerlcadict[i]["info"]["best_contradiction"]["domain"],singlemarkerlcadict[i]["info"]["best_contradiction"]["phylum"])
			tempdict["bestcontraident"] = singlemarkerlcadict[i]["info"]["best_contradiction"]["identity"]
			tempdict["bestcontrascore"] = singlemarkerlcadict[i]["info"]["best_contradiction"]["score"]
			
			tempdict["comparevalue"] = tempdict["besthitident"] + tempdict["bestcontraident"] # the marker with the highest overall identity of besthit + bestcontradiction combined is chosen as the representative. If two or more are tied, the one with the highest supportcount is chosen.
			tempdict["supportcount"] = 1
			obc_keytuple = tuple(sorted([tempdict["besthitdomphyl"],tempdict["bestcontradomphyl"]]))
			
			if obc_keytuple in sm_contradiction_overview: 
				if sm_contradiction_overview[obc_keytuple]["comparevalue"] >= tempdict["comparevalue"]:
					sm_contradiction_overview[obc_keytuple]["supportcount"] += 1
					continue
				else:
					tempdict["supportcount"] += sm_contradiction_overview[obc_keytuple]["supportcount"]
			sm_contradiction_overview[obc_keytuple] = tempdict
			
		if len(sm_contradiction_overview) > 0:
			rep_key = sorted(list(sm_contradiction_overview.keys()), key=lambda x: (sm_contradiction_overview[x]["comparevalue"], sm_contradiction_overview[x]["supportcount"]))[0]
			outinfo["sm_representative_contradiction"] = sm_contradiction_overview[rep_key]
			outinfo["sm_contradiction_counts"] = sum([sm_contradiction_overview[x]["supportcount"] for x in sm_contradiction_overview])
				
		outtaxpath, top2_contras, top2_contras_avidents, top2_contras_avscores = lca.weighted_lca(db, blasthitlist=tempsinglemarkerlcas, taxlevel=markerlevel, return_contradicting_top2 = True)
		
		
		if top2_contras != None:
			if outinfo["discrepancy_taxlevel"] == None:
				outinfo["discrepancy_taxlevel"] = db.taxid2taxlevel(top2_contras[0])
			if outinfo["single_or_weighted_discrep"] == "single":
				outinfo["single_or_weighted_discrep"] = "single&weighted"
			else:
				outinfo["single_or_weighted_discrep"] = "weighted"
			outinfo["weighted_lca_top_contradictions"] = {	"besthit_taxid" : top2_contras[0], "besthitident" : top2_contras_avidents[0], "besthitscore" : top2_contras_avscores[0], \
															"bestcontra_taxid" : top2_contras[1], "bestcontraident" : top2_contras_avidents[1], "bestcontrascore" : top2_contras_avscores[1]} #had to correct sorting of taxoptions in weighted lca, but now the order is actually correct (besthit vs best contradiction)
		
		
		amb_type, amb_evidence, amb_infotext = categorize(outinfo, singlemarkerlcadict)
		outinfo["amb_type"] = amb_type
		outinfo["amb_evidence"] = amb_evidence
		outinfo["amb_infotext"] = amb_infotext
		
		outinfo["sm_domain_phylum_classification_counts"] = get_singlehit_domain_phylum_counts(db, markernames, singlemarkerlcadict)
		# ~ if breakme: #todo: debugging only
			# ~ import pdb; pdb.set_trace()
		return outinfo
		
class blastdata_subset(blastdata_baseobject):
	def __init__(self, blastdata_obj,*args, query_id, max_evalue = None, min_ident = None, score_cutoff_fraction = None, keep_max_hit_fraction = None, keep_min_hit_count = None, auxilliary = False, seqtype=None, blacklist=None):
		assert isinstance(blastdata_obj, blastdata), "\nERROR: input must be an object of type 'blastdata', not '{}'\n".format(type(blastdata_obj))
		super().__init__(max_evalue = max_evalue, min_ident = min_ident, score_cutoff_fraction = score_cutoff_fraction, keep_max_hit_fraction = keep_max_hit_fraction, keep_min_hit_count = keep_min_hit_count, auxilliary = auxilliary, seqtype = seqtype, blacklist = blacklist)
		self.blastlinelist = blastdata_obj.pop_blastlines_for_query(query_id)
		self.get_settings_from_blastdata_obj(blastdata_obj)
		
	def get_settings_from_blastdata_obj(self, blastdata_obj): #todo: TEST THIS
		for x in ["max_evalue", "min_ident", "score_cutoff_fraction", "keep_max_hit_fraction", "keep_min_hit_count"]:
			if getattr(self, x) == None:
				setattr(self, x, getattr(blastdata_obj, x))
		 
class blastdata(blastdata_baseobject):
	def __init__(self, *blastfiles, max_evalue = None, min_ident = None, score_cutoff_fraction = 0.75, keep_max_hit_fraction = 0.5, keep_min_hit_count = 2, continue_from_json = False, auxilliary = False, seqtype=None, blacklist=None):
		super().__init__(max_evalue = max_evalue, min_ident = min_ident, score_cutoff_fraction = score_cutoff_fraction, keep_max_hit_fraction = keep_max_hit_fraction, keep_min_hit_count = keep_min_hit_count, auxilliary = auxilliary, seqtype = seqtype, blacklist = blacklist)
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
				
def read_blast_tsv(infilename, max_evalue = None, min_ident = None, dbobj = None, bindata_obj = None): #todo: apaprently obsolete due to blastdata-ojects. delete this function here if that is true
	""" 
	returns a a list of dictionaries, each representing the data in a blast line
	if max_evalue or min_ident are set, elines above or below these cutoffs are ignored
	each blast will later be assigned to a contig and to a "subject-type" (=stype), which can be either of ["16S", "23S", "univ_marker", "proc_marker", "bact_marker", "arch_marker", "other"]
	"""
	#note to self: wanted to use namedtuples here, but namedtuples won't work well here, because i may need them to be mutable.
	columninfos = {"query" : 0, "subject" : 1, "ident" : 2, "alignlen" : 3, "qstart": 6, "qend" : 7, "sstart" : 8, "send" : 9, "evalue" : 10, "score" : 11, "contig" : None, "stype" : None, "taxid": None} #Since python 3.7 all dicts are now ordered by default (YAY!). So the order of keys is maintained, without having to keep a seperate "orderlist". --> THEREFORE:  Ensure/Enforce python 3.7+ !!!
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


def get_contig_from_blastdb(seqid, blastdb, blastdbcmd = "blastdbcmd"): #todo: move this funciton to blasthandler
	'''
	todo: a different functions should also be added to the databaseobject created at getdb, that then calls this blasthandler funciton for the correct blastdb
	seqid may be a single seqid or a list of seqids
	returns a list of Seqrecords
	'''
	import subprocess
	from io import StringIO
	from Bio import SeqIO
	if type(seqid) == list:
		seqid = ",".join(seqid)
	cmd_blastdbcmd = ["blastdbcmd", "-entry", seqid, "-db", blastdb] # todo: add 'blastdbcmd' to dependencies and to configs
	# ~ print("huhuhuhuhu")
	# ~ print(cmd_blastdbcmd)
	# ~ print("hehehehehe")
	p_blastdbcmd = subprocess.run(cmd_blastdbcmd, stdout = subprocess.PIPE, stderr = subprocess.PIPE, text = True)
	try:
		p_blastdbcmd.check_returncode()
	except Exception:
		sys.stderr.write("\nAn error occured when trying to extract seqids {} from blast database {}\n".format(seqid, blastdb))
		sys.stderr.write("{}\n".format(p_blastdbcmd.stderr))
		raise RuntimeError
	fasta_io = StringIO(p_blastdbcmd.stdout)
	return(list(SeqIO.parse(fasta_io, "fasta")))

def _add_contigs2blasthits_later(blastlinelist, parsetype = "prodigal", lookup_table = None): #todo: probably obsolete...
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
	
def _add_stype2blasthits_later(blastlinelist, markersetdict): #todo: probably obsolte... #markersetdict should be derived from a different module "fastahandler.py" that reads the markerfastas and returns either the sequences, or just the locus_tags/fasta_headers.
	#"markersetlist" should be a list of sets. sets should be ordered in decreasing hierarchy (16S first, then 23S then markerprots then totalprots)
	#each set contains the locus_tags of the markers used in the blast
	for blindex in range(len(blastlinelist)):
		for ms in markersetdict:
			if blastlinelist[bl]["query"] in markersetdict[ms]:
				blastlinelist[bl]["stype"] = ms
	return blastlinelist

def _add_taxid2blasthits_later(blastlinelist, db_obj): #todo: probably obsolte...
	for blindex in range(len(blastlinelist)):
		blastlinelist[bl]["taxid"] = db_obj.acc2taxid(blastlinelist[bl]["taxid"]["subject"])
	return blastlinelist

def run_single_blast(query, db, blast, outfmt = "6", outname = None, threads = 1): #todo: add 'max_target_seqs" to optinal arguments (with 500 as default value)
	import subprocess
	assert os.path.basename(blast) in ["blastn", "blastp", "blastx"]
	assert outname, "\n\tERROR: must supply an outfilename!\n"
	assert outfmt in supported_outfmts, "\n\tError: outfmt '{}' not supported! can only be one of : {}".format(outfmt, ", ".join(supported_outfmts))
	if type(query) == str and os.path.isfile(query):
		inputarg = None
	elif type(query) == list:
		inputarg =  "\n".join([record.format("fasta") for record in query])
		query = "-"
	else:
		raise IOError("\nERROR: don't recognize query argument\n")
	blastcmd = subprocess.run([blast, "-query", query, "-db", db, "-evalue", "1e-10",\
							   "-outfmt", outfmt, "-num_threads", str(int(threads)), "-out", outname + ".tmp"], \
							  input = inputarg, stdout = subprocess.PIPE, stderr = subprocess.PIPE, text = True)
	try:
		blastcmd.check_returncode()
	except Exception:
		sys.stderr.write("\nAn error occured during blastp run with query '{}'\n".format(query))
		sys.stderr.write("{}\n".format(blastcmd.stderr))
		raise RuntimeError
	return outname

def run_single_blastp(query, db, blast, outfmt = "6", outname = None, threads = 1):
	assert os.path.basename(blast) == "blastp" 
	return run_single_blast(query, db, blast, outfmt, outname, threads)

def run_single_blastn(query, db, blast, outfmt = "6", outname = None, threads = 1):
	assert os.path.basename(blast) == "blastn"
	return run_single_blast(query, db, blast, outfmt, outname, threads)

def run_single_blastx(query, db, blast, outfmt = "6", outname = None, threads = 1):
	assert os.path.basename(blast) == "blastx"
	return run_single_blast(query, db, blast, outfmt, outname, threads)
	
def run_single_diamondblastp(query, db, diamond, outfmt = "6", outname = None, threads = 1): #TODO: currently not setting "--tmpdir" & "--parallel-tmpdir" here! figure something out if this turns out to be problematic on hpc systems
	#TODO: add a maxmem arguemt that states how much memory can be used. use this to determine optimal blocksize and chunks for more efficient blasting. BLOCKSIZE=INT(MEMORY/6) CHUNKS=4/2/1 IF MEMORY >=12/24/48
	#todo: increase default number of max_target hits stored to ~ 500 (and make it an optional argument)
	import subprocess
	assert outname, "\n\tERROR: must supply an outfilename!\n"
	assert outfmt in supported_outfmts, "\n\tError: outfmt '{}' not supported! can only be one of : {}".format(outfmt, ", ".join(supported_outfmts))
	blastcmd = subprocess.run([diamond, "blastp", "--query", query, "--db", db, "--evalue", "1e-10",\
							   "--outfmt", outfmt, "--threads", str(int(threads)), "--out", outname + ".tmp"], \
							   stdout = subprocess.PIPE, stderr = subprocess.PIPE, text = True)
	try:
		# ~ print(" ".join(blastcmd.args))
		blastcmd.check_returncode()
	except Exception:
		sys.stderr.write("\nAn error occured during diamond blastp run with query '{}'\n".format(query))
		sys.stderr.write("{}\n".format(blastcmd.stderr))
		raise RuntimeError
	return outname

def run_single_diamondblastx(query, db, diamond, outfmt = "6", outname = None, threads = 1): #TODO: currently not setting "--tmpdir" & "--parallel-tmpdir" here! figure something out if this turns out to be problematic on hpc systems
	#TODO: add a maxmem arguemt that states how much memory can be used. use this to determine optimal blocksize and chunks for more efficient blasting. BLOCKSIZE=INT(MEMORY/6) CHUNKS=4/2/1 IF MEMORY >=12/24/48
	import subprocess
	assert outname, "\n\tERROR: must supply an outfilename!\n"
	assert outfmt in supported_outfmts, "\n\tError: outfmt '{}' not supported! can only be one of : {}".format(outfmt, ", ".join(supported_outfmts))
	tempname = None
	import re
	outfmt = re.sub(" std ", " qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore ", outfmt).strip().split() #blast allows "std" keyword in outfmt, diamond doesn't
	if type(query) == list: #list indicates a list of seqrecords
		# ~ print("QUERY IS A LIST")
		tempname="delme_tempinput_mdmcleaner_diamondblastx.faa" # todo: figure out a better solution until a future diamond version hopefully may accept input from stdin...
		misc.write_fasta(query, tempname)
		query = tempname
	# ~ else:
		# ~ print("QUERY IS OF TYPE '{}'".format(type(query)))
	cmdlist = [diamond, "blastx", "--query", query, "--db", db, "--evalue", "1e-10",\
			"--outfmt", *outfmt, "--threads", str(int(threads)), "--out", outname + ".tmp"]
	# ~ print(cmdlist)
	blastcmd = subprocess.run(cmdlist, stdout = subprocess.PIPE, stderr = subprocess.PIPE, text = True)
	try:
		# ~ print(" ".join(blastcmd.args))
		blastcmd.check_returncode()
		if tempname:
			os.remove(tempname)
	except Exception:
		sys.stderr.write("\nAn error occured during diamond blastp run with query '{}'\n".format(query))
		sys.stderr.write("{}\n".format(blastcmd.stderr))
		if tempname:
			os.remove(tempname)
		raise RuntimeError
	return outname

def _run_any_blast(query, db, path_appl, outfmt = "6",outname = None ,threads = 1, force = False):
	#TODO: fix stupid problem that blast (unfortunately) cannot handle gzipped query files. solution read in compressed fastas, pipe records to blast-functions via stdin (TODO: add this function also to diamond (already present in ncbi-blast)
	appl = os.path.basename(path_appl)
	appl_onlypath = os.path.dirname(path_appl)
	#print("\nblasting --> {}".format(outname))
	assert outname, "\n\tERROR: must supply an outfilename!\n" 
	assert outfmt in supported_outfmts, "\n\tError: outfmt '{}' not supported! can only be one of : {}".format(outfmt, ", ".join(supported_outfmts))
	if os.path.isfile(outname) and force == False:
		sys.stderr.write("\n\tWARNING: blast result file '{}' already exists and 'force' not set to True --> skipping this blast\n".format(outname))
		return outname
	assert appl in ["blastp", "blastn", "diamond", "diamond blastp", "diamond blastx"], "\nError: unknown aligner '{}'\n".format(appl)
	_command_available(command=path_appl)
	if appl == "blastp":
		outname = run_single_blastp(query, db, path_appl, outfmt, outname,  threads)
	elif appl == "blastn":
		outname = run_single_blastn(query, db, path_appl, outfmt, outname, threads)
	elif appl in ["diamond", "diamond blastp"]:
		outname = run_single_diamondblastp(query, db, os.path.join(appl_onlypath, "diamond"), outfmt, outname, threads)
	elif appl == "diamond blastx":
		outname = run_single_diamondblastx(query, db, os.path.join(appl_onlypath, "diamond"), outfmt, outname, threads)
	os.rename(outname + ".tmp", outname)
	return outname

def _distribute_threads_over_jobs(total_threads, num_jobs): # to distribute N threads over M groups as evenly as possible, put (N/M) +1 in (N mod M) groups, and (N/M) in the rest
	##TODO: test using misc.run_multiple_functions_parallel() for this instead! DELETE THIS IF MISC VERSION WORKS!
	if num_jobs == 0:
		sys.stderr.write("\tnothing to do...\n")
		return [], total_threads
	if num_jobs < total_threads:
		more_threads_list = [ (total_threads / num_jobs) + 1 ] * (total_threads % num_jobs) #these should go to the lower priority markers, as there will be more of them
		fewer_threads_list = [ (total_threads / num_jobs) ] * (num_jobs - len(more_threads_list))
		return fewer_threads_list + more_threads_list, num_jobs
	return [1] * num_jobs, total_threads

def run_multiple_blasts_parallel(basic_blastarg_list, *, outfmt = "6", outbasename=None, total_threads=1): #basic_blastarg_list = list of tuples such as [(query1, db1, blast1), (query2, db2, blast2),...])
	# ~ import pdb; pdb.set_trace()
	#blasic_blastarg_list should be a list of tuples, each tuple containg NOTHING else but (query, blastdb) in THAT exact order (!) #todo: find a better way to do this. perhaps dicts instead of tuples (and enforce keyword aguments for blast functions)?
	if len(basic_blastarg_list) == 0:
		sys.stderr.write("\tnothing to blast...")
	from multiprocessing import Pool
	thread_args, no_processes = _distribute_threads_over_jobs(total_threads, len(basic_blastarg_list))
	arglist = []
	# ~ print("---")
	# ~ print(outfmt)
	for i in range(len(basic_blastarg_list)):
		# ~ print(basic_blastarg_list[i])
		if type(basic_blastarg_list[i][0]) == list:
			# ~ print("blast input is a list with {} entries".format(len(basic_blastarg_list[i][0])))
			querystring = "auxblasts"
		else:
			querystring = os.path.basename(basic_blastarg_list[i][0])
		outname = "{prefix}_{counter:03d}{appl}_{query}_vs_{db}.tsv".format(prefix = outbasename, counter = i, \
					appl = os.path.basename(basic_blastarg_list[i][2]).replace(" ", "_"), query = querystring, \
					db = os.path.basename(basic_blastarg_list[i][1]))
		arglist.append(tuple(bba for bba in basic_blastarg_list[i]) + (outfmt, outname, thread_args[i]))
		# ~ print(arglist[-1])
	# ~ print(arglist)
	masterblaster = Pool(processes = no_processes)
	# ~ import pdb; pdb.set_trace()
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
	sys.stderr.write("\n\tcreating diamond-db {}...".format(outfilename))
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
	sys.stderr.write("\n\tcreating blastdb {}...".format(outfilename))
	sys.stderr.flush()
	mkdbcmd = [ makeblastdb, "-in", infasta, "-title", outfilename, "-parse_seqids", "-out", outfilename, "-dbtype", db_type, "-max_file_sz", "4GB"] #hash_index apparently not possible for very large datasets
	# ~ print(" ".join(mkdbcmd)) #todo: add verbose option and print mkdbcmd it TRUE
	mkdbproc = subprocess.run(mkdbcmd, stdout = subprocess.PIPE, stderr = subprocess.PIPE, text = True)
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
	if not infasta.endswith(".gz"):
		# ~ sys.stderr.write("\ninput is not compressed! using standard make_blast_db function".format(outfilename))
		return make_blast_db(infasta, outfilename, makeblastdb, db_type, threads) #todo: make this the other way round. make it call this function if it finds input to be compressed
	assert db_type in ["nucl", "prot"],  "dbtype has to be either 'nucl' or 'prot'"
	sys.stderr.write("\n\tcreating blastdb {}...".format(outfilename))
	sys.stderr.flush()
	zcatcmd = [ "zcat", infasta]
	mkdbcmd = [ makeblastdb, "-title", outfilename, "-parse_seqids", "-out", outfilename, "-dbtype", db_type, "-max_file_sz", "4GB"] #hash_index apparently not possible for very large datasets
	# ~ import time
	# ~ start = time.time()
	zcatproc = subprocess.Popen(zcatcmd, stdout = subprocess.PIPE, stderr = subprocess.PIPE, text = True)
	# ~ print(" ".join(mkdbcmd)) #todo: add verbose option and print mkdbcmd it TRUE
	mkdbproc = subprocess.Popen(mkdbcmd, stdin = zcatproc.stdout, stdout = subprocess.PIPE, stderr = subprocess.PIPE, text = True)
	stdout, stderr = mkdbproc.communicate()
	mkdbproc.wait()
	sys.stderr.write("...finished!\n")
	# ~ end = time.time()
	# ~ print("piping it this way took : {}".format(end - start))
	# ~ sys.stdout.flush()
	return outfilename
	#todo: compare how long it take to pipe this via python subprocess and ow long to do the same on the bash shell


def get_blast_combinations(dblist, querylist, blast = "blastn"):
	import itertools #todo move up
	# ~ print("gbc")
	# ~ import pdb; pdb.set_trace()
	querylist = [x for x in querylist if x != ""] #workaround for cases were there is no query file (e.g. no ssu or lsu rRNA found)
	acceptable_blasttools = ["blastn", "blastp", "blastx"]
	assert os.path.basename(blast) in acceptable_blasttools, "Error: 'blast' must be one of {}".format(acceptable_blasttools) #maybe remove this check or add a workaround to use diamond whith this also?
	return [ blasttuple + (blast,) for blasttuple in list(itertools.chain(*list(zip(querylist, permu) for permu in itertools.permutations(dblist, len(querylist)))))] #todo:Only works as long as dblist is longer or equal to querylist... find a way to ensure/workaround this

def _command_available(command="diamond"): #todo: move this to misc, maybe?
	'''
	for checking if diamond or blast are actually installed in the current environment
	'''
	import subprocess
	try:
		subprocess.call([command.split()[0]],stderr=subprocess.PIPE, stdout=subprocess.PIPE)
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
		jobtuplelist.append((".blasthandler", "_run_any_blast", {"query":query, "db":db, "path_appl": tool, "outname":"{}_{}_vs_{}.tab.tmp".format(tool, os.path.basename(query), os.path.basename(db))}))
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



