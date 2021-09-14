#!/usr/bin/env python
""" 
creating carious overview- and report-files for mdm_cleaner
"""
import misc
from misc import openfile

def write_refdb_ambiguity_report(magsag, ambiguities, outfile):
	import io
	output = ""
	if len(ambiguities) > 0:
		header = dict2tsvline(indictionary, lineprefix = "magsag\t", key_header = "contig", onlyheader=True)
		if isinstance(outfile, str):
			output += header
			outfile = openfile(outfile, "wt")
			# ~ outfile.write(header)
		if isinstance(outfile, io.IOBase):
			for i in ambiguities:
				output += dict2tsvline(indictionary, lineprefix = "{}\t".format(magsag), key_header = "contig", onlyheader=True)
				line = "{}\t{}\t{}\n".format(magsag, i, "\t".join([str(v) for v in ambiguities[i].values()]))
				outfile.write(line)
		# ~ import pdb; pdb.set_trace()
	return outfile

def gather_extended_bin_metrics(bindata, outfile, cutoff=5): #todo: make a simple_binmetrics function
	def write_dictlines(indict, outfile):
		import io
		# ~ print("{} = {}".format(outfile, type(outfile)))
		# ~ sys.stdout.flush()
		# ~ sys.stderr.flush()
		if isinstance(outfile, str):
			# ~ print("IS A STRING --> CREATING NEW FILE")
			outfile = openfile(outfile, "wt")
			outfile.write("#{}\n".format("\t".join(list(indict.keys()))))
		if isinstance(outfile, io.IOBase):
			# ~ print("IS A FILE --> APPENDING")
			line = "{}\n".format("\t".join([";".join([str(y) for y in indict[x]]) if type(indict[x]) == list else str(indict[x]) for x in indict ]))
			# ~ print(outfile.name)
			# ~ print(line)
			outfile.write(line)
		return outfile
	
	
	print("WRITING TO OVERVIEWFILE")
	#bc_ = before cleanup
	#ac_ = after cleanup
	import statistics
	#todo: this is awfull and convoluted! Improve & simplify or get rid of this before final release!
	binname = bindata.bin_tempname
	totalbincontigs = len(bindata.contigdict)
	if totalbincontigs > 0:
		totalbinbp = bindata.get_total_size()
		majortaxpath = bindata.get_consensus_taxstringlist()
		fraction_trustedbp = sum([ bindata.contigdict[contig]["contiglen"] for contig in bindata.get_trusted_contignames() ])/totalbinbp
		fraction_unknownbp = sum([ bindata.contigdict[contig]["contiglen"] for contig in bindata.get_contignames_with_trustscore(5) ])/totalbinbp
		fraction_untrustedbp = sum([ bindata.contigdict[contig]["contiglen"] for contig in bindata.get_untrusted_contignames() ])/totalbinbp
		bin_trust = bindata.trust_index_from_tax_score(statistics.mean( [ bindata.contigdict[contig]["tax_score"] for contig in bindata.contigdict]))
		bin_trust_ignoring_viral = bindata.trust_index_from_tax_score(statistics.mean( [ bindata.contigdict[contig]["tax_score"] for contig in bindata.contigdict if not bindata.contigdict[contig]["viral"]] ))
		fraction_different_species = sum([ bindata.contigdict[contig]["contiglen"] for contig in bindata.contigdict if bindata.contigdict[contig]["contradict_consensus"] == "species"])/totalbinbp
		fraction_different_genus = sum([ bindata.contigdict[contig]["contiglen"] for contig in bindata.contigdict if bindata.contigdict[contig]["contradict_consensus"] == "genus"])/totalbinbp
		fraction_different_family = sum([ bindata.contigdict[contig]["contiglen"] for contig in bindata.contigdict if bindata.contigdict[contig]["contradict_consensus"] == "family"])/totalbinbp
		fraction_different_order = sum([ bindata.contigdict[contig]["contiglen"] for contig in bindata.contigdict if bindata.contigdict[contig]["contradict_consensus"] == "order"])/totalbinbp
		fraction_different_class = sum([ bindata.contigdict[contig]["contiglen"] for contig in bindata.contigdict if bindata.contigdict[contig]["contradict_consensus"] == "class"])/totalbinbp
		fraction_different_phylum = sum([ bindata.contigdict[contig]["contiglen"] for contig in bindata.contigdict if bindata.contigdict[contig]["contradict_consensus"] == "phylum"])/totalbinbp
		fraction_different_domain = sum([ bindata.contigdict[contig]["contiglen"] for contig in bindata.contigdict if bindata.contigdict[contig]["contradict_consensus"] == "domain"])/totalbinbp
		fraction_refdb_contamination = bindata.get_fraction_refdbcontamination()
		fraction_nocoding = bindata.get_fraction_nocoding()
		fraction_viral = sum([ bindata.contigdict[contig]["contiglen"] for contig in bindata.contigdict if bindata.contigdict[contig]["viral"]])/totalbinbp
		fraction_trust0 = sum([ bindata.contigdict[contig]["contiglen"] for contig in bindata.get_contignames_with_trustscore(0) ])/totalbinbp
		fraction_trust1 = sum([ bindata.contigdict[contig]["contiglen"] for contig in bindata.get_contignames_with_trustscore(1) ])/totalbinbp
		fraction_trust2 = sum([ bindata.contigdict[contig]["contiglen"] for contig in bindata.get_contignames_with_trustscore(2) ])/totalbinbp
		fraction_trust3 = sum([ bindata.contigdict[contig]["contiglen"] for contig in bindata.get_contignames_with_trustscore(3) ])/totalbinbp
		fraction_trust4 = sum([ bindata.contigdict[contig]["contiglen"] for contig in bindata.get_contignames_with_trustscore(4) ])/totalbinbp
		fraction_trust5 = sum([ bindata.contigdict[contig]["contiglen"] for contig in bindata.get_contignames_with_trustscore(5) ])/totalbinbp
		fraction_trust6 = sum([ bindata.contigdict[contig]["contiglen"] for contig in bindata.get_contignames_with_trustscore(6) ])/totalbinbp
		fraction_trust7 = sum([ bindata.contigdict[contig]["contiglen"] for contig in bindata.get_contignames_with_trustscore(7) ])/totalbinbp
		fraction_trust8 = sum([ bindata.contigdict[contig]["contiglen"] for contig in bindata.get_contignames_with_trustscore(8) ])/totalbinbp
		fraction_trust9 = sum([ bindata.contigdict[contig]["contiglen"] for contig in bindata.get_contignames_with_trustscore(9) ])/totalbinbp
		fraction_trust10 = sum([ bindata.contigdict[contig]["contiglen"] for contig in bindata.get_contignames_with_trustscore(10) ])/totalbinbp
		total_16SrRNA = len([gene for gene in bindata.markerdict if  getmarkers.seqid2contig(gene) in bindata.contigdict and bindata.markerdict[gene]["stype"] == "ssu_rRNA"])
		total_23SrRNA = len([gene for gene in bindata.markerdict if  getmarkers.seqid2contig(gene) in bindata.contigdict and bindata.markerdict[gene]["stype"] == "lsu_rRNA"])
		total_5SrRNA = len([gene for gene in bindata.markerdict if  getmarkers.seqid2contig(gene) in bindata.contigdict and bindata.markerdict[gene]["stype"] == "tsu_rRNA"])
		total_trna = len([gene for gene in bindata.markerdict if  getmarkers.seqid2contig(gene) in bindata.contigdict and bindata.markerdict[gene]["stype"] == "trna"])
		total_marker_prots = len([gene for gene in bindata.markerdict if getmarkers.seqid2contig(gene) in bindata.contigdict and "_marker " in bindata.markerdict[gene]["stype"]])
		total_proteins = len([gene for gene in bindata.markerdict if getmarkers.seqid2contig(gene) in bindata.contigdict and bindata.markerdict[gene]["stype"] == "total"]) + total_marker_prots
	else:
		# ~ print("hererherherherherhehre")
		totalbinbp = 0
		majortaxpath = 0
		fraction_trustedbp = 0
		fraction_unknownbp = 0
		fraction_untrustedbp = 0
		bin_trust = None
		bin_trust_ignoring_viral = None
		fraction_different_species = 0
		fraction_different_genus = 0
		fraction_different_family = 0
		fraction_different_order = 0
		fraction_different_class = 0
		fraction_different_phylum = 0
		fraction_different_domain = 0
		fraction_refdb_contamination = 0
		fraction_nocoding = 0
		fraction_viral = 0
		fraction_trust0 = 0
		fraction_trust1 = 0
		fraction_trust2 = 0
		fraction_trust3 = 0
		fraction_trust4 = 0
		fraction_trust5 = 0
		fraction_trust6 = 0
		fraction_trust7 = 0
		fraction_trust8 = 0
		fraction_trust9 = 0
		fraction_trust10 = 0
		total_16SrRNA = 0
		total_23SrRNA = 0
		total_5SrRNA = 0
		total_trna = 0
		total_marker_prots = 0
		total_proteins = 0

	print_dict = {  "binname" : binname, \
							"totalbincontigs" : totalbincontigs,\
							"totalbinbp" : totalbinbp,\
							"majortaxpath" : majortaxpath,\
							"fraction_trustedbp" : fraction_trustedbp, \
							"fraction_unknownbp" : fraction_unknownbp, \
							"fraction_untrustedbp" : fraction_untrustedbp, \
							"bin_trust" : bin_trust, \
							"bin_trust_ignoring_viral":bin_trust_ignoring_viral, \
							"fraction_different_species":fraction_different_species, \
							"fraction_different_genus":fraction_different_genus, \
							"fraction_different_family":fraction_different_family, \
							"fraction_different_order":fraction_different_order, \
							"fraction_different_class":fraction_different_class, \
							"fraction_different_phylum":fraction_different_phylum, \
							"fraction_different_domain":fraction_different_domain, \
							"fraction_viral":fraction_viral, \
							"fraction_refdb_contamination" : fraction_refdb_contamination, \
							"fraction_nocoding" : fraction_nocoding, \
							"fraction_trust0":fraction_trust0, \
							"fraction_trust1":fraction_trust1, \
							"fraction_trust2":fraction_trust2, \
							"fraction_trust3":fraction_trust3, \
							"fraction_trust4":fraction_trust4, \
							"fraction_trust5":fraction_trust5, \
							"fraction_trust6":fraction_trust6, \
							"fraction_trust7":fraction_trust7, \
							"fraction_trust8":fraction_trust8, \
							"fraction_trust9":fraction_trust9, \
							"fraction_trust10":fraction_trust10, \
							"total_16SrRNA":total_16SrRNA, \
							"total_23SrRNA":total_23SrRNA, \
							"total_5SrRNA":total_5SrRNA, \
							"total_trna":total_trna, \
							"total_marker_prots":total_marker_prots, \
							"total_proteins":total_proteins }
							
	return write_dictlines(print_dict, outfile)
