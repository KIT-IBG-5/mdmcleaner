#!/usr/bin/env python
from Bio import SeqIO
import misc
from misc import openfile
import sys
import os
import blasthandler
import re
'''
much of this will move to the other modules when finished
'''

class comparison_hit(object):
	def __init__(self,*, taxid, seqid, domain, phylum, db, markerlevel):
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
		
	def set_extractdb(self):
		if self.seqsource == "silva":
			if self.markerlevel == "ssu_rRNA_tax":
				self.extractdb = self.db.nucdbs_ssu_rRNA[-1]
			elif self.markerlevel == "lsu_rRNA_tax":
				self.extractdb = self.db.nucdbs_lsu_rRNA[-1]
		elif self.seqsource == "gtdb":
			self.extractdb = self.db.nucdbs_genome[0]
	
	def get_seqrecord(self, blastdbcmd = "blastdbcmd"):
		import blasthandler
		self.seqrecord = blasthandler.get_contig_from_blastdb(self.seqid, self.extractdb, blastdbcmd)

	def blast_contigs(self, threads = 1, blacklist=None, outfileprefix = ""):
		import time
		if self.blast != None:
			basic_blastarglist = [ (self.seqrecord, blastdb, self.blast) for blastdb in self.blastdbs ]
			print(basic_blastarglist)
			outbasename = "{}refblast_{}_".format(outfileprefix, self.seqid)
			print(outbasename)
			start = time.time()
			resultfiles = blasthandler.run_multiple_blasts_parallel(basic_blastarglist, outfmt= "6 std qlen slen", outbasename=outbasename, total_threads=threads)
			end = time.time()
			print("blastfiles:\n\t{}\n\t".format("\n\t".join(resultfiles)))
			print("\nthis blast took {:.4f} seconds\n\n".format(end-start))
			self.blastdata = blasthandler.blastdata(*resultfiles, max_evalue = 1e-5, min_ident = 90, score_cutoff_fraction = 0, keep_max_hit_fraction = 1, keep_min_hit_count = 2, continue_from_json = False, auxilliary = False, seqtype=None, blacklist=blacklist) #todo: ignlorelistfile should be changed to ignore list. Can be a list of a filepath. If filepath, read that file as list. if list, use that
			self.blastdata.add_info_to_blastlines(taxdb_obj=self.db)
			# ~ delmeprint(self.blastdata.blastlinelist, self.db)
			print("="*50)
			if self.blast == "blastx":
				self.blastdata.filter_blasthits_by_cov_and_ident(mincov=90, filterbylen=subject)
			else:
				self.blastdata.filter_blasthits_by_cov_and_ident()		
			# ~ print("afterfilter:")
			# ~ delmeprint(self.blastdata.blastlinelist, self.db)
			# ~ print("="*50)
			return_category = count_contradictions(self.blastdata, self.db, self.domain, self.phylum, self.seqid)
			print("*"*100)
			return return_category

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
		self.blast = "blastn" #except, if the best contradiction is a "refseqprot" = eukaryotic or viral seqeunce --> then change to "diamond blastx"
		self.blastdbs = db.nucdbs_genome
		self.get_seqrecord()
		# ~ self.blast_contigs(threads = configs["threads"], blacklist=None, outfileprefix = "")
	
	def changeblastmethod(self, blasttype):
		assert blasttype in ["blastn", "blastx", "diamond blastx"], "\nERROR: blasttype={}\n".format(blasttype)
		if blasttype in ["blastx", "diamond blastx"]:
			self.blast = "diamond blastx"
			self.blastdbs = self.db.protdbs_all
		elif blasttype in ["blastn"]:
			self.blast = "blastn"
			self.blastdbs = dself.b.nucdbs_genome
				
class comp_ssurrna(comparison_hit):
	def __init__(self,*, taxid, seqid, domain, phylum, db, markerlevel, configs):
		super().__init__(taxid=taxid,seqid=seqid, domain=domain, phylum=phylum, db=db, markerlevel=markerlevel)
		self.extractdb = db.nucdbs_ssu_rRNA[-1] #todo: may need to terate through genomic db also in some cases?
		self.blastdbs = db.nucdbs_ssu_rRNA
		self.seqtype = "nucl"
		self.blast = "blastn"
		self.get_seqrecord()
		# ~ self.blast_contigs(threads = configs["threads"], blacklist=None, outfileprefix = "")
	
class comp_lsurrna(comparison_hit):
	def __init__(self,*, taxid, seqid, domain, phylum, db, markerlevel, configs):
		super().__init__(taxid=taxid,seqid=seqid, domain=domain, phylum=phylum, db=db, markerlevel=markerlevel)
		self.extractdb = db.nucdbs_lsu_rRNA[-1] #todo: may need to terate through genomic db also in some cases?
		self.blastdbs = db.nucdbs_lsu_rRNA
		self.seqtype = "nucl"
		self.blast = "blastn"
		self.get_seqrecord()
		# ~ self.blast_contigs(threads = configs["threads"], blacklist=None, outfileprefix = "")

class comparison_pair(object):
	prot_domain_types = ["d__Eukaryota", "None"] #Eukaryotes and Viruses (the latter yield Domain = "None") currently only have protein DBs
	nuc_marker_types = ["ssu_rRNA_tax", "lsu_rRNA_tax"]
	sm_lca_bh_pattern = re.compile("sm_best hit=\'([^\(\);\']+)\'\(([^;,\(\)]+),([^;,\(\)]+);\s*acc=\'([^;,\(\)]+)\'")
	sm_lca_bc_pattern = re.compile("sm_best contradiction=\'([^\(\);\']+)\'\(([^;,\(\)]+),([^;,\(\)]+);\s*acc=\'([^;,\(\)]+)\'")
	blacklist = None
	added2blacklistcount = 0
	# \1 = taxid, \2 = domain; \3= phylum, \4= accession, 
	# virus_indicator_pattern = re.compile("=\'eukcat__viral'") #in case it is needed to differentiate viral besthits/contradictions from those of other sources that also happen to yield "domain = None" (I don't think there are any such cases...)

	def __init__(self, amb_evidence,markerlevel, db, configs, blacklist = None, outfileprefix = ""):
		self.db = db
		self.markerlevel = markerlevel
		self.amb_evidence = amb_evidence
		self.best_hit, self.best_contradiction = self.parse_evidence(configs)
		if type(self).blacklist == None:
			if blacklist == None:
				type(self).blacklist = set()
			else:
				type(self).blacklist = blacklist
		for i in self.best_hit, self.best_contradiction: #todo. do this in calling function. if diamond blastx is planned, collect instances and blast together
			return_category = i.blast_contigs(threads = configs["threads"], blacklist=blacklist, outfileprefix = outfileprefix)
			if return_category == "contamination":
				type(self).blacklist.add(i.seqid)
				type(self).added2blacklistcount += 1
		
	def parse_evidence(self, configs): #todo: this should move to a compare_group_object
		return_list = []
		for p in [self.sm_lca_bh_pattern, self.sm_lca_bc_pattern]:
			patternhit = re.search(p, self.amb_evidence)
			seqid = patternhit.group(4)
			domain = patternhit.group(2)
			phylum = patternhit.group(3)
			taxid = patternhit.group(1)
			comp = None
			if domain in self.prot_domain_types and self.markerlevel not in self.nuc_marker_types:
				comp = comp_refseqprot(taxid = taxid, seqid=seqid, domain=domain, phylum=phylum, db=self.db, markerlevel=self.markerlevel, configs=configs)
			elif self.markerlevel == "ssu_rRNA_tax":
				comp = comp_ssurrna(taxid = taxid,seqid=seqid, domain=domain, phylum=phylum, db=self.db, markerlevel=self.markerlevel, configs=configs)
			elif self.markerlevel == "lsu_rRNA_tax":
				comp = comp_lsurrna(taxid = taxid,seqid=seqid, domain=domain, phylum=phylum, db=self.db, markerlevel=self.markerlevel, configs=configs)
			else:
				comp = comp_prokprotcontig(taxid = taxid,seqid=seqid, domain=domain, phylum=phylum, db=self.db, markerlevel=self.markerlevel, configs=configs)
			return_list.append(comp)
		if "prot" in [ x.seqtype for x in return_list]:
			for x in range(len(return_list)):
				if not isinstance(return_list[x], comp_refseqprot): #TODO: NOTE: extraction protein sequences from diamond DBs is rather inefficient. so this will be skipped for now, only noted here in case it should be implemented in a later version (if it turns out to be needed, which is unlikely considering only get the accession of one protein per referencecontig this way)
					return_list[x].changeblastmethod("diamond blastx") #TODO: collect blastx-jobs and run all together, instead of running one after another (is faster)
		# ~ import pdb; pdb.set_trace()
		assert len(return_list) == 2, "\nERROR: something went wrong when parsing amb_evidence: \n'{}'\n\n".format(self.amb_evidence)
		return return_list[0], return_list[1]
		


def delmeprint(inlist, db):
	wantedkeys = ["query", "subject", "ident", "alignlen", "qstart", "qend", "sstart", "send", "score", "qlen", "slen"]
	for line in inlist:
		outline="\t".join(str(line[x]) for x in wantedkeys)
		outline+="\t{:.3f}".format(line["alignlen"]/min(line["qlen"],line["slen"])*100)
		dp_tuple = db.get_domain_phylum(line["taxid"])
		outline+="\t{}\n".format(dp_tuple)
		print(outline)
		print("---")

def read_ambiguity_report(ambiguity_report, configs, blacklist = None):
	import getdb
	'''
		todo: ambiguity/contradiction at phylum level (--> e.g.: both hits are bacteria, but different phyla) are easy. But what about contradicitons at domain level? How do i doublecheck if eukaryotic, when all I have is PROTEINS for eukaryotes?)  
		todo: ANSWER: use complete amgiguity info as input:
			if sm-LCA level: check if best hit and besz contradiction come from bacteria and/or archaea --> extract contigs and blastn nucleotide level
			if however one of them come from eukaryotes and/or viruses --> just extract prokaryotic contigs and blastx against refprots (check if eukaryotes yield more/better hits than prokaryotes of same domain/phylum).
				in BOTH cases: the domain/phylum that yields the more than twice as many hits (or twice as high summed score?) within a specified identity cutoff (use total_prots species cutoff = 85%? ) than the other "wins"
				todo: in future the db should keep track which genomes are fom isolates, and wich are from MAGs/SAGs. Isolates should count double. Ideally, fully close genome should count triple
			IGNORE weighted -LCA instances for now. create specific chimera check for those later!
	'''
	import re
	sm_contam_pattern = re.compile("potential refDB-contamination \[\w+ indication sm-LCA level\]")

	print("\nREADING: {}\n\n".format(ambiguity_report))
	db = getdb.taxdb(configs)
	collect_diamond_querys = {}	
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
			cp = comparison_pair(amb_evidence, markerlevel, db, configs=configs)
			if cp.best_hit == cp.best_contradiction == None:
				continue
			print("\nadded {} new entries to blacklist!\n".format(cp.added2blacklistcount))
			print("current blacklist : \n\t{}".format("\n\t".join([x for x in cp.blacklist])))
			#TODO: NOTE: extraction protein sequences from diamond DBs is rather inefficient. so this will be skipped for now (mayble iplemented in a later version if it turns out it is needed). Will only analyse proteins for now
			#todo: blast comparison seqs against appropriate blastdb (e.g. blastx vs combined_refprots if one is eukaryotic, otherwise concat_refgenomes (with appropriate program))
			#todo: choose contig that yields higher summed-hitscores to domain/phylum other than annotated as contamination. YIeld warnfing if BOTH show such an result
			# ~ blastdatalist = blast_refdb_contigs(comparison_seqdicts, dbobject=db, outfile_basename = "refdb_ambig_{}".format(str(counter).zfill(5)), blacklist=blacklist)
			#todo: 1.) parse & filter blast results, 2.) assign taxa per blastline, 3.) count domains and phyla
			# ~ print("blastfiles:\n\t{}\n\t".format("\n\t".join(blastfiles))
			# ~ blastdata = blasthandler.blastdata(blastfiles, max_evalue = 1e-5, min_ident = 90, score_cutoff_fraction = 0, keep_max_hit_fraction = 1, keep_min_hit_count = 2, continue_from_json = False, auxilliary = False, seqtype=None, ignorelistfile=None) #todo: ignlorelistfile should be changed to ignore list. Can be a list of a filepath. If filepath, read that file as list. if list, use that
			# ~ import pdb; pdb.set_trace()

def count_contradictions(blastdata, db, comparison_domain, comparison_phylum, query_acc): #todo. add query-id to required arguments, so that hits to original query can be more safely recognized and ignored! 
	import re
	comparison_domain = re.sub("^d__","",comparison_domain)
	comparison_phylum = re.sub("^p__","",comparison_phylum)
	counts = {}
	weights = {}
	best_selfscore = 0
	firstline = True
	for line in blastdata.blastlinelist:
		# ~ print(line)
		# ~ print(firstline)
		dp_tuple = db.get_domain_phylum(line["taxid"])
		if firstline and line["subject"] == query_acc:
			best_selfscore = line["score"]
			firstline = False
		if dp_tuple in counts:
			counts[dp_tuple] += 1
			weights[dp_tuple] += line["score"]
		else:
			counts[dp_tuple] = 1
			if line["subject"] == query_acc:
				weights[dp_tuple] = 0 #ignoring scores of hits to self
			else: 
				weights[dp_tuple] = line["score"]
	print("expected: {}, {}".format(comparison_domain, comparison_phylum))
	print("best_selfscore: {}".format(best_selfscore))
	print("actual counts:\n\t{}".format("\n\t".join(["{} : {}".format(key, counts[key]) for key in counts.keys() ])))
	print("actual weights:\n\t{}".format("\n\t".join(["{} : {}".format(key, weights[key]) for key in weights.keys() ])))
	print("----------------------------")
	domain_counts = { domain: sum([counts[x] for x in counts.keys() if x[0] == domain]) for domain in set([y[0] for y in counts.keys()]) }
	domain_weights = { domain: sum([weights[x] for x in weights.keys() if x[0] == domain]) for domain in set([y[0] for y in weights.keys()]) }
	print("domain_counts:\n\t{}".format("\n\t".join(["{} : {}".format(key, domain_counts[key]) for key in domain_counts.keys() ])))
	domain_counts_expected = domain_counts[comparison_domain]
	domain_counts_contradicting = sum([domain_counts[x] for x in domain_counts.keys() if x != comparison_domain])
	print("domain_counts expected == {}".format(domain_counts_expected))
	print("domain_counts contradicting == {}".format(domain_counts_contradicting))
	domain_weights_expected = domain_weights[comparison_domain]
	domain_weights_contradicting = sum([domain_weights[x] for x in domain_weights.keys() if x != comparison_domain])
	print("domain_weights:\n\t{}".format("\n\t".join(["{} : {}".format(key, domain_weights[key]) for key in domain_weights.keys() ])))
	print("domain_weights expected == {}".format(domain_weights_expected))
	print("domain_weights contradicting == {}".format(domain_weights_contradicting))
	print("----------------------------")
	phylum_counts = { phylum: sum([counts[x] for x in counts.keys() if x[1] == phylum]) for phylum in set([y[1] for y in counts.keys()]) }
	phylum_weights= { phylum: sum([weights[x] for x in weights.keys() if x[1] == phylum]) for phylum in set([y[1] for y in weights.keys()]) }
	print("phylum_counts:\n\t{}".format("\n\t".join(["{} : {}".format(key, phylum_counts[key]) for key in phylum_counts.keys() ])))
	phylum_counts_expected = phylum_counts[comparison_phylum]
	phylum_counts_contradicting = sum([phylum_counts[x] for x in phylum_counts.keys() if x != comparison_phylum])
	print("phylum_counts expected == {}".format(phylum_counts_expected))
	print("phylum_counts contradicting == {}".format(phylum_counts_contradicting))
	phylum_weights_expected = phylum_weights[comparison_phylum]
	phylum_weights_contradicting = sum([phylum_weights[x] for x in phylum_weights.keys() if x != comparison_phylum])
	print("phylum_weights:\n\t{}".format("\n\t".join(["{} : {}".format(key, phylum_weights[key]) for key in phylum_weights.keys() ])))
	print("phylum_weights expected == {}".format(phylum_weights_expected))
	print("phylum_weights contradicting == {}".format(phylum_weights_contradicting))
	# if contradicting counts larger than expectedhitcounts AND contradicting weights >= 2x expectedhitcounts: --> contamination
	if domain_counts_contradicting > domain_counts_expected and domain_weights_contradicting >= (domain_weights_expected * 2) and domain_weights_contradicting >= (best_selfscore * 0.2): #the last requirement ensures that random matchtes of very small contigs representing e.g. partial transposases etc do not cause misclassification of very large contigs that match only in those short regions (in other words: very large contigs are only seen as contaminants if other large contigs match)
		print(" ---> '{}' is a contamination on domain level! --> ADD TO BLACKLIST!".format(query_acc))
		return "contamination"
	# if contradcting counts >=  expectedhitcounts AND contradicting weights >= expectedhitwights: suspicious ambiguity! one may be contamination, but not sure which
	if domain_counts_contradicting >= domain_counts_expected and domain_weights_contradicting >= domain_weights_expected and domain_weights_contradicting >= (best_selfscore * 0.1):
		print(" '{}' has hits to genomes of other domains! either it is a contamination, or the matching contigs in other domains are! --> evaluate independently (e.g. blast against refseq or uniprot)".format(query_acc))
		return "ambiguity"
	if phylum_counts_contradicting > phylum_counts_expected and phylum_weights_contradicting >= (phylum_weights_expected * 2) and phylum_weights_contradicting >= (best_selfscore * 0.2): #the last requirement ensures that random matchtes of very small contigs representing e.g. partial transposases etc do not cause misclassification of very large contigs that match only in those short regions (in other words: very large contigs are only seen as contaminants if other large contigs match)
		print(" ---> '{}' is a contamination on phylum level! --> ADD TO BLACKLIST!".format(query_acc))
		return "contamination"
	# if contradcting counts >=  expectedhitcounts AND contradicting weights >= expectedhitwights: suspicious ambiguity! one may be contamination, but not sure which
	if phylum_counts_contradicting >= phylum_counts_expected and phylum_weights_contradicting >= phylum_weights_expected and phylum_weights_contradicting >= (best_selfscore * 0.1):
		print(" '{}' has hits to genomes of other phyla! either it is a contamination, or the matching contigs in other phyla are! --> evaluate independently (e.g. blast against refseq or uniprot)".format(query_acc))
		return "ambiguity"
	# if contradicting counts < expectedhitcounts OR contrdicting weights < expectedhitweights: actually ok! --> drop ambiguity
	if phylum_counts_contradicting < phylum_counts_expected or phylum_counts_contradicting < phylum_weights_expected:
		print("NO sufficient proof for '{}' being a contminant. assuming it is OK!".format(query_acc))
		return "OK"
	print("there is a case i have not considered, and THIS is it!")
	#todo: add detailed evidence strings to 'contamination' and 'ambiguity' classifications...
	import pdb; pdb.set_trace()


# ~ def reevaluate_ambigeous_contigs():
	# ~ '''
	# ~ todo: this should be a function of bindata-objects
		# ~ todo: ANSWER: use complete amgiguity info as input:
			# ~ if sm-LCA level: check if best hit and besz contradiction come from bacteria and/or archaea --> extract contigs and blastn nucleotide level
			# ~ if however one of them come from eukaryotes and/or viruses --> just extract prokaryotic contigs and blastx against refprots (check if eukaryotes yield more/better hits than prokaryotes of same domain/phylum).
				# ~ in BOTH cases: the domain/phylum that yields the more than twice as many hits (or twice as high summed score?) within a specified identity cutoff (use total_prots species cutoff = 85%? ) than the other "wins"
				# ~ todo: in future the db should keep track which genomes are fom isolates, and wich are from MAGs/SAGs. Isolates should count double. Ideally, fully close genome should count triple
			# ~ IGNORE weighted -LCA instances for now. create specific chimera check for those later!	'''
	# ~ pass
