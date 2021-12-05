#!/usr/bin/env python
from Bio import SeqIO
import misc
from misc import openfile
'''
much of this will move to the other modules when finished
'''

def get_contig_from_blastdb(seqid, blastdb, blastdbcmd = "blastdbcmd"):
	'''
	this function should go directly to blasthandler.py
	a different functions should also be added to the databaseobject created at getdb, that then calls this blasthandler funciton for the correct blastdb
	seqid may be a single seqid or a list of seqids
	returns a list of Seqrecords
	'''
	import subprocess
	from io import StringIO
	if type(seqid) == list:
		seqid = ",".join(seqid)
	cmd_blastdbcmd = ["blastdbcmd", "-entry", seqid, "-db", blastdb] # todo: add 'blastdbcmd' to dependencies and to configs
	p_blastdbcmd = subprocess.run(cmd_blastdbcmd, stdout = subprocess.PIPE, stderr = subprocess.PIPE, text = True)
	try:
		p_blastdbcmd.check_returncode()
	except Exception:
		sys.stderr.write("\nAn error occured when trying to extract seqids {} from blast database {}\n".format(seqid, blastdb))
		sys.stderr.write("{}\n".format(p_blastdbcmd.stderr))
		raise RuntimeError
	fasta_io = StringIO(p_blastdbcmd.stdout)
	return(list(SeqIO.parse(fasta_io, "fasta")))



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
	sm_contam_pattern = re.compile("potential refDB-contamination [\w+ indication sm-LCA level]")
	sm_lca_bh_pattern = re.compile(".*?sm_best hit=\'[\w\d\. -_]+\'\(d__([\w\d\. _-]+),p__([\w\d\. _-]+);acc=\'([\w\d\. _-]+)\'.*")
	sm_lca_bc_pattern = re.compile(".*?sm_best contradiction=\'[\w\d\. -_]+\'\(d__([\w\d\. _-]+),p__([\w\d\. _-]+);acc=\'([\w\d\. _-]+).*?")
	# \1 = besthit_domain; \2= besthit_phylum, \3= besthit_accession, \4 = best_contradiction_domain, \5 = best_contradiciton_phylum, \6 = best_contradiction_accesion
	# virus_indicator_pattern = re.compile("=\'eukcat__viral'") #in case it is needed to differentiate viral besthits/contradictions from those of other sources that also happen to yield "domain = None" (I don't think there are any such cases...)

	def parse_evidence(amb_evidence, markertype, extractdb):
		return_dictlist = []
		prot_domain_types = ["d__Eukaryota", "None"] #Eukaryotes and Viruses (the latter yield Domain = "None") currently only have protein DBs
		nuc_marker_types = ["ssu_rRNA_tax", "lsu_rRNA_tax"]
		for p in [sm_lca_bh_pattern, sm_lca_bc_pattern]:
			seqid_dict = {}
			patternhit = re.search(sm_lca_pattern, amb_evidence)
			seqid_dict["seqid"] = patternhit.group(3)
			seqid_dict["domain"] = patternhit.group(1)
			seqid_dict["phylum"] = patternhit.group(2)
			seqid_dict["markertype"] = markertype
			if domain in prot_domain_types and markertype not in nuc_marker_types:
				seqid_dict["seqtype"] = "prot" #TODO: NOTE: extraction protein sequences from diamond DBs is rather inefficient. so this will be skipped for now, only noted here in case it should be implemented in a later version (if it turns out to be needed, which is unlikely considering only get the accession of one protein per referencecontig this way)
				seqid_dict["seqrecord"] = None
			else:
				seqid_dict["seqtype"] = "nucl"
				seqid_dict["seqrecord"] = get_contig_from_blastdb(seqid_dict["seqid"], extractdb, blastdbcmd = configs["blastdbcmd"][0])
			return_dictlist.append(seqid_dict)
		return return_dictlist # should return list with two entries: 1.) the best_hit, 2.) the best_contradiction

	def get_contig_from

	db = getdb.taxdb(configs)			
	with openfile(ambiguity_report) as infile:
		for line in infile:
			tokens = line.strip().split("\t")
			if line.startswith("magsag\t") or len(tokens) == 0:
				continue
			magsag = tokens[0]
			contig = tokens[1]
			markerlevel = tokens[2]
			### TODO: the selection of the correct blastdb should be part of "parse_evidence!" becaus
			if markerlevel == "ssu_rRNA_tax":
				extractdb = db.nucdb_ssu_rRNA #todo add blastdb paths to db-object #TODO: this actually depends if it is a silva entry or not
				blastdbs = [db.nucdb_ssu_rRNA, db.nucdb_genome]
			elif markerlevel == "lsu_rRNA_tax":
				extractdb = db.nucdb_lsu_rRNA #todo add blastdb paths to db-object
				blastdbs = [db.nucdb_lsu_rRNA, db.nucdb_genome]
			else:
				extractdb = db.nucdb_genome #todo add blastdb paths to db-object	
				blastdbs = [db.nucdb_genome]				
			ambtype = tokens[6]
			if re.search(sm_contam_pattern, ambtype) == None: #ignore everything other than potential contaminations detected on singlemarker level for now (those exclusively found on weighted LCA level are more indicative for chimeras than refDB contaminations...)
				continue
			amb_evidence = tokens[8]
			comparison_seqdicts = parse_evidence(amb_evidence, markerlevel)
			if comparison_seqsdicts == None:
				continue
			# todo: comparison_seqs is a bit convoluted! each entry (should be max 2 though) is a tuple, consisting of two entries: a.) the seqrecord, b.) the corresponding seqid_tuple
			#TODO: NOTE: extraction protein sequences from diamond DBs is rather inefficient. so this will be skipped for now (mayble iplemented in a later version if it turns out it is needed). Will only analyse proteins for now
			#todo: blast comparison seqs against appropriate blastdb (e.g. blastx vs combined_refprots if one is eukaryotic, otherwise concat_refgenomes (with appropriate program))
			#todo: choose contig that yields higher summed-hitscores to domain/phylum other than annotated as contamination. YIeld warnfing if BOTH show such an result
			blast_refdb_contigs(comparison_seqdicts, blastdbs=blastdbs, dbobject=db, blacklist=blacklist)
			

def blast_refdb_contigs(comparison_seqtuples, blastdbs, dbobject, blacklist=None):
	for cs in comparison_seqtuples:
		#blast sequence against database listed in blastdbs (withqlen and seqlen fieds, though) 
		#filter hits based on 90% identity and at least 50% coverage of the shorter contig
		#count domains/phyla of hits and compare counts to asumed contig domain/phylum in comparison_seqtuple cd[1][-1]
		
	


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


def main():
	print(get_contig_from_blastdb(["GCF_000812565.1_NZ_JTAM01000055.1","GCF_001266905.1_NZ_LHYY01000050.1"], "./deltestdb"))

main()
