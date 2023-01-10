#!/usr/bin/env python
from Bio import SeqIO
from mdmcleaner import misc
from mdmcleaner.misc import openfile
import sys
import os
from mdmcleaner import blasthandler
import re
import tempfile
'''
much of this will move to the other modules when finished
'''

sm_contam_pattern = re.compile("potential refDB-contamination \[\w+ indication sm-LCA level\]")

class comparison_hit(object):
	def __init__(self,*, taxid, seqid, domain, phylum, db, markerlevel, settings=None):
		if settings == None:
			self.settings = {x:x for x in ["blastn", "blastp", "diamond"]}
		else:
			self.settings = settings
		self.seqid = seqid
		self.domain = domain
		self.phylum = phylum
		self.db = db
		self.taxid = taxid
		self.markerlevel = markerlevel
		self.taxpath = db.taxid2pathstring(taxid)
		self.seqsource = db._gtdb_refseq_or_silva(seqid)
		self.seqrecord = None
		self.blastlinelist = None
		self.extractdb = None
		self.blastdbs = None
		self.blast = None
		self.blastdata = None
		
	# ~ def __hash__(self): #in case i want to use the objects themselves as dict or hash keys...
		# ~ return hash((self.seqid, self.blast, self.markerlevel))
		
	# ~ def __eq__(self, other):
		# ~ return (self.seqid, self.blast, self.markerlevel) == (other.seqid, other.blast, other.markerlevel)
		
	# ~ def __ne__(self, other): #apparently not necessary anymore (see accepted answer of https://stackoverflow.com/questions/4352244/should-ne-be-implemented-as-the-negation-of-eq), but still implemented here just to be sure...
		# ~ return not self == other
	def return_key(self):
		'''
		returns a tuple of (seqid, blast, markerlevel) for assigning blast hits to the correct objects after running collective blast jobs
		'''
		return (self.seqid, self.blast, self.markerlevel)

	def set_extractdb(self):
		if self.seqsource == "silva":
			if self.markerlevel == "ssu_rRNA_tax":
				self.extractdb = self.db.nucdbs_ssu_rRNA[-1]
			elif self.markerlevel == "lsu_rRNA_tax":
				self.extractdb = self.db.nucdbs_lsu_rRNA[-1]
		elif self.seqsource == "gtdb":
			self.extractdb = self.db.nucdbs_genome[0]
	
	def get_seqrecord(self, blastdbcmd = "blastdbcmd"):
		from mdmcleaner import blasthandler
		self.seqrecord = blasthandler.get_contig_from_blastdb(self.seqid, self.extractdb, blastdbcmd)

	def blast_contigs(self, threads = 1, blacklist=None, outfileprefix = ""):
		# ~ import time
		if self.blast != None:
			basic_blastarglist = [ (self.seqrecord, blastdb, self.blast) for blastdb in self.blastdbs ]
			outbasename = "{}refblast_{}_".format(outfileprefix, self.seqid)
			# ~ start = time.time()
			resultfiles = blasthandler.run_multiple_blasts_parallel(basic_blastarglist, outfmt= "6 std qlen slen", outbasename=outbasename, total_threads=threads)
			# ~ end = time.time()
			# ~ print("blastfiles:\n\t{}\n\t".format("\n\t".join(resultfiles)))
			# ~ print("\nthis blast took {:.4f} seconds\n\n".format(end-start))
			# ~ print(resultfiles)
			self.blastdata = blasthandler.blastdata(*resultfiles, max_evalue = 1e-5, min_ident = 90, score_cutoff_fraction = 0, keep_max_hit_fraction = 1, keep_min_hit_count = 2, continue_from_json = False, auxilliary = False, seqtype=None, blacklist=blacklist) #todo: ignlorelistfile should be changed to ignore list. Can be a list of a filepath. If filepath, read that file as list. if list, use that
			if self.blast.endswith("blastx"):
				self.blastdata.filter_blasthits_by_cov_and_ident(mincov=90, filterbylen=subject)
			else:
				self.blastdata.filter_blasthits_by_cov_and_ident() #todo: stricter identity cutoffs for ssu-rRNA	
			self.blastdata.add_info_to_blastlines(taxdb_obj=self.db, verbose=False)			
			return_category, return_note = self.count_contradictions()
			return {"evaluation": return_category, "markerlevel_checked": [self.markerlevel], "note": return_note}, resultfiles


	def count_contradictions(self): # argument "ch" should be a "comparison_hit"-object
	#blastdata, db, comparison_domain, comparison_phylum, query_acc): #todo. add query-id to required arguments, so that hits to original query can be more safely recognized and ignored! 
		import re
		comparison_domain = re.sub("^d__","", self.domain)
		comparison_phylum = re.sub("^p__","",self.phylum)
		if comparison_domain == "None":
			comparison_domain = None
		if comparison_phylum == "None":
			comparison_phylum = None
		counts = {}
		weights = {}
		best_selfscore = 0
		firstline = True
		for line in self.blastdata.blastlinelist:
			dp_tuple = self.db.get_domain_phylum(line["taxid"])
			if firstline and line["subject"] == self.seqid:
				best_selfscore = line["score"]
				firstline = False
			if dp_tuple in counts:
				counts[dp_tuple] += 1
				weights[dp_tuple] += line["score"]
			else:
				counts[dp_tuple] = 1
				if line["subject"] == self.seqid:
					weights[dp_tuple] = 0 #ignoring scores of hits to self
				else: 
					weights[dp_tuple] = line["score"]
		# ~ print("expected: {}, {}".format(comparison_domain, comparison_phylum))
		# ~ print("best_selfscore: {}".format(best_selfscore))
		# ~ print("actual counts:\n\t{}".format("\n\t".join(["{} : {}".format(key, counts[key]) for key in counts.keys() ])))
		# ~ print("actual weights:\n\t{}".format("\n\t".join(["{} : {}".format(key, weights[key]) for key in weights.keys() ])))
		# ~ print("----------------------------")
		# ~ if self.seqid == "GCA_002842085.1_PHCA01000074.1": #todo:only for debugging
			# ~ import pdb; pdb.set_trace()
		domain_counts = { domain: sum([counts[x] for x in counts.keys() if x[0] == domain]) for domain in set([y[0] for y in counts.keys()]) }
		domain_weights = { domain: sum([weights[x] for x in weights.keys() if x[0] == domain]) for domain in set([y[0] for y in weights.keys()]) }
		# ~ print("domain_counts:\n\t{}".format("\n\t".join(["{} : {}".format(key, domain_counts[key]) for key in domain_counts.keys() ])))
		# ~ try:
		domain_counts_expected = domain_counts[comparison_domain] #todo: only in try_except statement for debugging
		# ~ except Exception:
			# ~ import pdb; pdb.set_trace()
		domain_counts_contradicting = sum([domain_counts[x] for x in domain_counts.keys() if x != comparison_domain])
		# ~ print("domain_counts expected == {}".format(domain_counts_expected))
		# ~ print("domain_counts contradicting == {}".format(domain_counts_contradicting))
		domain_weights_expected = domain_weights[comparison_domain]
		domain_weights_contradicting = sum([domain_weights[x] for x in domain_weights.keys() if x != comparison_domain])
		domain_weights_contradicting_excluding_none = sum([domain_weights[x] for x in domain_weights.keys() if (x != comparison_domain and x != None) ])
		# ~ print("domain_weights:\n\t{}".format("\n\t".join(["{} : {}".format(key, domain_weights[key]) for key in domain_weights.keys() ])))
		# ~ print("domain_weights expected == {}".format(domain_weights_expected))
		# ~ print("domain_weights contradicting == {}".format(domain_weights_contradicting))
		# ~ print("----------------------------")
		phylum_counts = { phylum: sum([counts[x] for x in counts.keys() if x[1] == phylum]) for phylum in set([y[1] for y in counts.keys()]) }
		phylum_weights= { phylum: sum([weights[x] for x in weights.keys() if x[1] == phylum]) for phylum in set([y[1] for y in weights.keys()]) }
		# ~ print("phylum_counts:\n\t{}".format("\n\t".join(["{} : {}".format(key, phylum_counts[key]) for key in phylum_counts.keys() ])))
		# ~ if self.seqid == "AF391990.1.1401":
			# ~ import pdb; pdb.set_trace()
		phylum_counts_expected = phylum_counts[comparison_phylum]
		phylum_counts_contradicting = sum([phylum_counts[x] for x in phylum_counts.keys() if x != comparison_phylum])
		# ~ print("phylum_counts expected == {}".format(phylum_counts_expected))
		# ~ print("phylum_counts contradicting == {}".format(phylum_counts_contradicting))
		phylum_weights_expected = phylum_weights[comparison_phylum]
		phylum_weights_contradicting = sum([phylum_weights[x] for x in phylum_weights.keys() if x != comparison_phylum ])
		phylum_weights_contradicting_excluding_none = sum([phylum_weights[x] for x in phylum_weights.keys() if (x != comparison_phylum and x != None)])
		# ~ print("phylum_weights:\n\t{}".format("\n\t".join(["{} : {}".format(key, phylum_weights[key]) for key in phylum_weights.keys() ])))
		# ~ print("phylum_weights expected == {}".format(phylum_weights_expected))
		# ~ print("phylum_weights contradicting == {}".format(phylum_weights_contradicting))
		# ~ print("phylum_weights contradicting (excluding 'None')== {}".format(phylum_weights_contradicting_excluding_none))
		for pc in [domain_counts, phylum_counts]:
			if len(pc) == 2 and None in pc.keys():
				if "rRNA" in self.markerlevel:
					note = " --> '{}' represents a contradiction between Silva and GTDB taxonomies!".format(self.seqid)
					return "ambiguity", note # SHOULD RETURN "silva_conflict"!  only TEMPORARY FIX for Iss #37! Still need to implement preferrence of gtdb annotation in such cases!
				else:
					note = " --> there appear to be references similar to {} that are not annotated up to phylum level! This is likely an error in the reference DB!".format(self.seqid)
					return "wtf", note
		# if contradicting counts larger than expectedhitcounts AND contradicting weights >= 2x expectedhitcounts: --> contamination
		if domain_counts_contradicting > domain_counts_expected and domain_weights_contradicting >= (domain_weights_expected * 2) and domain_weights_contradicting >= (best_selfscore * 0.2): #the last requirement ensures that random matchtes of very small contigs representing e.g. partial transposases etc do not cause misclassification of very large contigs that match only in those short regions (in other words: very large contigs are only seen as contaminants if other large contigs match)
			note = " ---> '{}' shows indication of being a contamination on domain level! --> ADD TO BLACKLIST!".format(self.seqid)
			return "contamination", note
		# if contradcting counts >=  expectedhitcounts AND contradicting weights >= expectedhitwights: suspicious ambiguity! one may be contamination, but not sure which
		if domain_counts_contradicting >= domain_counts_expected and domain_weights_contradicting >= domain_weights_expected and domain_weights_contradicting >= (best_selfscore * 0.1):
			note = " '{}' has highly similar hits to genomes of other domains! either it or the matching contigs from other domains may be contaminations! --> requires independent evaluation (e.g. blast against refseq, uniprot or nr)".format(self.seqid)
			return "ambiguity", note
		if phylum_counts_contradicting > phylum_counts_expected and phylum_weights_contradicting >= (phylum_weights_expected * 2) and phylum_weights_contradicting >= (best_selfscore * 0.2): #the last requirement ensures that random matchtes of very small contigs representing e.g. partial transposases etc do not cause misclassification of very large contigs that match only in those short regions (in other words: very large contigs are only seen as contaminants if other large contigs match)
			note = " ---> '{}' shows indication of being a contamination on phylum level! --> ADD TO BLACKLIST!".format(self.seqid)
			return "contamination", note
		# if contradcting counts >=  expectedhitcounts AND contradicting weights >= expectedhitwights: suspicious ambiguity! one may be contamination, but not sure which
		if phylum_counts_contradicting >= phylum_counts_expected and phylum_weights_contradicting >= phylum_weights_expected and phylum_weights_contradicting >= (best_selfscore * 0.1):
			note = " '{}' has highly similar hits to genomes of other phyla! either it or the matching contigs in other phyla may be contaminations! --> requires independent evaluation (e.g. blast against refseq, uniprot or nr)".format(self.seqid)
			return "ambiguity", note
		# if contradicting counts < expectedhitcounts OR contrdicting weights < expectedhitweights: actually ok! --> drop ambiguity
		if phylum_counts_contradicting < phylum_counts_expected or phylum_weights_contradicting < phylum_weights_expected or phylum_weights_contradicting < (best_selfscore * 0.1):
			note = "NO sufficient proof for '{}' being a contaminant. assuming it is OK!".format(self.seqid)
			return "OK", note
		note = "there is a case i have not considered, and THIS is it!"
		return "wtf", note
		#todo: add detailed evidence strings to 'contamination' and 'ambiguity' classifications...
		# ~ import pdb; pdb.set_trace()

class comp_refseqprot(comparison_hit): #specifically for eukaryotic and viral references that are currently still only stored as protein sequences...
	def __init__(self,*, taxid, seqid, domain, phylum, db, markerlevel, configs):
		super().__init__(taxid=taxid,seqid=seqid, domain=domain, phylum=phylum, db=db, markerlevel=markerlevel)
		self.seqtype = "prot" #TODO: NOTE: extraction protein sequences from diamond DBs is rather inefficient. so this will be skipped for now, only noted here in case it should be implemented in a later version (if it turns out to be needed, which is unlikely considering only get the accession of one protein per referencecontig this way)
	
	def blast_contigs(self, threads = 1, blacklist=None, outfileprefix = ""):
			pass
			
class comp_prokprotcontig(comparison_hit):
	def __init__(self,*, taxid, seqid, domain, phylum, db, markerlevel, configs):
		super().__init__(taxid=taxid,seqid=seqid, domain=domain, phylum=phylum, db=db, markerlevel=markerlevel)
		self.set_extractdb()
		self.seqtype = "nucl" 
		self.blast = self.settings["blastn"] #except, if the best contradiction is a "refseqprot" = eukaryotic or viral seqeunce --> then change to "diamond blastx"
		self.blastdbs = db.nucdbs_genome
		self.get_seqrecord()
		# ~ self.blast_contigs(threads = configs["threads"], blacklist=None, outfileprefix = "")
	
	def changeblastmethod(self, blasttype):
		assert blasttype in ["blastn", "blastx", "diamond blastx"], "\nERROR: blasttype={}\n".format(blasttype)
		if blasttype in ["blastx", "diamond blastx"]:
			self.blast = self.settings["diamond"] + " blastx"
			self.blastdbs = self.db.protdbs_all
		elif blasttype in ["blastn"]:
			self.blast = self.settings["blastn"]
			self.blastdbs = db.nucdbs_genome
				
class comp_ssurrna(comparison_hit):
	def __init__(self,*, taxid, seqid, domain, phylum, db, markerlevel, configs):
		super().__init__(taxid=taxid,seqid=seqid, domain=domain, phylum=phylum, db=db, markerlevel=markerlevel)
		self.set_extractdb()
		# ~ self.extractdb = db.nucdbs_ssu_rRNA[-1] #todo: may need to terate through genomic db also in some cases?
		self.blastdbs = db.nucdbs_ssu_rRNA
		self.seqtype = "nucl"
		self.blast = self.settings["blastn"]
		self.get_seqrecord()
		# ~ self.blast_contigs(threads = configs["threads"], blacklist=None, outfileprefix = "")
	
class comp_lsurrna(comparison_hit):
	def __init__(self,*, taxid, seqid, domain, phylum, db, markerlevel, configs):
		super().__init__(taxid=taxid,seqid=seqid, domain=domain, phylum=phylum, db=db, markerlevel=markerlevel)
		self.set_extractdb()
		# ~ self.extractdb = db.nucdbs_lsu_rRNA[-1] #todo: may need to terate through genomic db also in some cases?
		self.blastdbs = db.nucdbs_lsu_rRNA
		self.seqtype = "nucl"
		self.blast = self.settings["blastn"]
		self.get_seqrecord()
		# ~ self.blast_contigs(threads = configs["threads"], blacklist=None, outfileprefix = "")

class suspicious_entries(object):
	'''
	meant to replace 'comparison_pair' objects, since "best_hits" and "best_contradictions" are now being reviewed seperately and independently anyway
	'''
	prot_domain_types = ["d__Eukaryota", "None"] #Eukaryotes and Viruses (the latter yield Domain = "None") currently only have protein DBs
	nuc_marker_types = ["ssu_rRNA_tax", "lsu_rRNA_tax"]
	sm_lca_bh_pattern = re.compile("sm_best hit=[a-z; ]{0,7}\'([^\(\);\']+)\'\(([^;,\(\)]+),([^;,\(\)]+);\s*acc=\'([^;,\(\)]+)\'")
	sm_lca_bc_pattern = re.compile("sm_best contradiction=[a-z; ]{0,7}\'([^\(\);\']+)\'\(([^;,\(\)]+),([^;,\(\)]+);\s*acc=\'([^;,\(\)]+)\'")

	
	def __init__(self, db, configs, outbasename = "./"):
		self.db = db
		self.threads = configs.settings["threads"]
		self.configs = configs.settings
		self.blacklist = configs.blacklist #todo: double check that the default value for "blacklist" in configs in set() and NOT None!
		self.blacklist_additions = set()
		self.added2blacklistcount = 0
		self.blastxjobs = {}
		self.blastxdbs = None
		self.seqid2evaluation = {} # key = seqid, value = dictionary --> { "evaluation": a, "markerlevel_checked": [b], "note" : c}
		self.outbasename = outbasename
		self.tempdir = tempfile.mkdtemp(prefix="tempdir_", suffix="_mdmcleaner_refdbcontam", dir=".")
		self.delete_filelist = []
		self.last_checked = []
	
	def last_checked_evaluations(self):
		'''
		returns the most severe evaluation result of the last checked reference-seqids
		'''
		last_evaluations =  [ self.seqid2evaluation[seqid]["evaluation"] for seqid in self.last_checked if seqid in self.seqid2evaluation ]
		# ~ print(last_evaluations)
		if "contamination" in last_evaluations:
			return "contamination"
		elif "ambiguity" in last_evaluations:
			return "ambiguity"
		elif "wtf" in last_evaluations:
			return "wtf"
		for le in last_evaluations:
			assert le == "OK", "\nERROR: evaluation result not accounted for: '{}'\n".format(le)
		return "OK"
		
	def evaluateornot(self, comp, blastxdone=False): #todo: improve this! This just a convoluted solution to make sure multiple mentions of contigs are only re-evaluated if they affect more/different target databases (e.g. "ssu_rRNA" is blasted against genomes AND SSU-databases, but not against LSU databases. So no need to redo genomeblasts on markerlevel "total_prots" or "marker_prots", but may ned to reblast "lsu_rRNA"
		# ~ if comp.seqid in self.blacklist: #unnecessary, because this is already checked in the calling function
			# ~ return #if it was already previously detected as contamination on any level, no need to do more blasts
		blastfilenames=[]
		evaluation = {"evaluation": None, "markerlevel_checked": [None], "note": None} 
		if comp.seqid not in self.seqid2evaluation or (blastxdone and self.seqid2evaluation[comp.seqid] != "contamination"): # todo: convoluted --> simplify. if rRNA dbs searched, no need for additional search of genome-DBs (were already included). but if protein-blast --> search against eukaryotes not in nucleotide-genome-DB --> do again
			if not blastxdone:
				evaluation, blastfilenames = comp.blast_contigs(self.threads, blacklist=self.blacklist, outfileprefix = os.path.join(self.tempdir, ""))
			else:
				return_category, return_note = comp.count_contradictions() #todo: redundant. streamline blastcontigs() and countcontradictions() more
				evaluation = {"evaluation": return_category, "markerlevel_checked": [comp.markerlevel], "note": return_note} 
			self.seqid2evaluation[comp.seqid] = evaluation
			if evaluation["evaluation"] == "contamination":
				self.blacklist.add(comp.seqid)	
				self.blacklist_additions.add(comp.seqid)
		else:
			if "rRNA" in comp.markerlevel and comp.markerlevel not in self.seqid2evaluation[comp.seqid]["markerlevel_checked"]:
				evaluation, blastfilenames = comp.blast_contigs(self.threads, blacklist=self.blacklist, outfileprefix = os.path.join(self.tempdir, ""))
				self.seqid2evaluation[comp.seqid]["markerlevel_checked"].append(comp.markerlevel)
				self.seqid2evaluation[comp.seqid]["note"] += "; {}".format(evaluation["note"])
		self.delete_filelist += blastfilenames
		return evaluation["evaluation"]
	
	def collective_diamondblast(self):
		# ~ print("{} potential diamond blasts".format(len(self.blastxjobs)))
		eval_list = []
		if len(self.blastxjobs) > 0:
			# ~ print([self.blastxjobs[x].seqid for x in self.blastxjobs])
			blastrecords = [self.blastxjobs[x].seqrecord[0] for x in self.blastxjobs if len(self.blastxjobs[x].seqrecord[0]) < 100000] #for now skipping reference contigs larger than 100 kb (takes too long to blastx). TODO: in such cases, search for ribosomal & other markergenes to verify classification!
			sys.stderr.write("\nblasting {} entries with blastx against reference proteins (another {} entries were too long to blastx efficiently\n".format(len(blastrecords), len(self.blastxjobs) - len(blastrecords)))
			if len(blastrecords) == 0:
				return
			basic_blastarglist = [ (blastrecords, blastdb, "diamond blastx") for blastdb in self.blastxdbs ]
			# ~ import pdb; pdb.set_trace()
			import time
			start = time.time()
			resultfiles = blasthandler.run_multiple_blasts_parallel(basic_blastarglist, outfmt= "6 std qlen slen", outbasename= os.path.join(self.tempdir, ""), total_threads=self.threads)
			self.delete_filelist += resultfiles
			end = time.time()
			# ~ print ("THIS BLAST TOOK {:.3f} seconds".format(end-start))
			collective_blastdata = blasthandler.blastdata(*resultfiles, max_evalue = 1e-5, min_ident = 90, score_cutoff_fraction = 0, keep_max_hit_fraction = 1, keep_min_hit_count = 2, continue_from_json = False, auxilliary = False, seqtype=None, blacklist=self.blacklist)
			for x in self.blastxjobs:
				if len(self.blastxjobs[x].seqrecord[0]) < 100000:
					self.blastxjobs[x].blastdata = blasthandler.blastdata_subset(collective_blastdata, query_id = x)
					self.blastxjobs[x].blastdata.add_info_to_blastlines(taxdb_obj=self.blastxjobs[x].db, verbose=False)
					eval_list.append(self.evaluateornot(self.blastxjobs[x], blastxdone = True))
		return eval_list
				
	def parse_evidence(self, amb_evidence, markerlevel): #todo: this should move to a compare_group_object
		return_list = []
		self.last_checked = []
		for p in [self.sm_lca_bh_pattern, self.sm_lca_bc_pattern]:
			patternhit = re.search(p, amb_evidence)
			seqid = patternhit.group(4)
			domain = patternhit.group(2)
			phylum = patternhit.group(3)
			taxid = patternhit.group(1)
			comp = None
			if seqid in self.blacklist:
				continue
			if domain in self.prot_domain_types and markerlevel not in self.nuc_marker_types:
				comp = comp_refseqprot(taxid = taxid, seqid=seqid, domain=domain, phylum=phylum, db=self.db, markerlevel=markerlevel, configs=self.configs)
				
			elif markerlevel == "ssu_rRNA_tax":
				comp = comp_ssurrna(taxid = taxid,seqid=seqid, domain=domain, phylum=phylum, db=self.db, markerlevel=markerlevel, configs=self.configs)
				
			elif markerlevel == "lsu_rRNA_tax":
				comp = comp_lsurrna(taxid = taxid,seqid=seqid, domain=domain, phylum=phylum, db=self.db, markerlevel=markerlevel, configs=self.configs)
			else:
				comp = comp_prokprotcontig(taxid = taxid,seqid=seqid, domain=domain, phylum=phylum, db=self.db, markerlevel=markerlevel, configs=self.configs)
			return_list.append(comp)
		if "prot" in [ x.seqtype for x in return_list ]:
			for x in range(len(return_list)):
				if not isinstance(return_list[x], comp_refseqprot): #TODO: NOTE: extraction protein sequences from diamond DBs is rather inefficient. so this will be skipped for now, only noted here in case it should be implemented in a later version (if it turns out to be needed, which is unlikely considering only get the accession of one protein per referencecontig this way)
					return_list[x].changeblastmethod("diamond blastx")
					self.blastxjobs[return_list[x].seqid] = return_list[x]
					self.blastxdbs = return_list[x].blastdbs #todo: this works as long as there is only one protein-db or if ALWAYS ALL available protein dbs are used (--> for any protein balst, ALWAYS the same set of DBs is used as currently the case). If that changes, make sure this here still works! 
					self.last_checked.append(return_list[x].seqid)
		else:
			for c in return_list:
				self.evaluateornot(c)
				self.last_checked.append(c.seqid)

def read_ambiguity_report(ambiguity_report, configs):
	import os
	from mdmcleaner import getdb
	import re
	
	# ~ outdir, outfileprefix = os.path.split(outbasename)
	# ~ if outdir != "" and not os.path.exists(outdir):
		# ~ os.mkdir(outdir)
	# ~ sm_contam_pattern = re.compile("potential refDB-contamination \[\w+ indication sm-LCA level\]")		
	db = getdb.taxdb(configs)
	suspects = suspicious_entries(db, configs)
	with openfile(ambiguity_report) as infile:
		counter = 0
		for line in infile:
			counter += 1
			tokens = line.strip().split("\t")
			if line.startswith("magsag\t") or len(tokens) == 0:
				continue
			magsag = tokens[0]
			contig = tokens[1]
			markerlevel = tokens[2]			
			ambtype = tokens[6]
			if re.search(sm_contam_pattern, ambtype) == None: #ignore everything other than potential contaminations detected on singlemarker level for now (those exclusively found on weighted LCA level are more indicative for chimeras than refDB contaminations...)
				continue
			amb_evidence = tokens[8]
			sys.stderr.write("\r\tprocessing line {}".format(counter))
			suspects.parse_evidence(amb_evidence, markerlevel)
	sys.stderr.write("\r\tprocessed line {} --> FINISHED!\n".format(counter))
	sys.stderr.write("\nnow running remaining diamond blasts collectively")
	suspects.collective_diamondblast()
	for df in suspects.delete_filelist:
		os.remove(df)
	return suspects.blacklist_additions #todo: write blacklist_additions to outfile progressively. also optionally write all other evaluations to logfile

