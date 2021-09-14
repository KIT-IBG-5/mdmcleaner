#!/usr/bin/env python
import sys
import os
import argparse
import getdb
import misc
from misc import openfile
import getdb
import getmarkers
import blasthandler
import pprint
import lca

from _version import __version__

progressdump_filename = "mdmprogress.json"

#todo: split gtdb protblastdbs into several subdbs that can be blasted in parallel (check out if faster first)
#todo: find commond names for gtdb and ncbi dbs to simplify things
#todo: find actual names for ncbi dbs
dbfiles = { "gtdb" : {	"protblastdbs" : ["gtdbplus_protdb.dmnd"], \
						"nucblastdbs" : ["concat_refgenomes", "SILVA_138.1_SSURef_NR99_tax_silva", "SILVA_138.1_LSURef_NR99_tax_silva"] ,\
						"ssu_nucblastdbs" : ["concat_refgenomes", "SILVA_138.1_SSURef_NR99_tax_silva"], \
						"lsu_nucblastdbs" : ["concat_refgenomes", "SILVA_138.1_LSURef_NR99_tax_silva"], \
						"trna_nucblastdbs" : ["concat_refgenomes"], \
						"mdmdbs" : ["gtdb_all.accession2taxid.sorted", "gtdb_taxonomy_br.json.gz", "gtdb_lcawalkdb_br.db"] },\
			"ncbi" : { "protblastdbs" : ["nr"], \
						"nucblastdbs" : ["nt"], \
						"ssu_nucblastdbs" : ["nt"], \
						"lsu_nucblastdbs" : ["nt"], \
						"trna_nucblastdbs" : ["nt"], \
						"mdmdbs" : [ "ncbi_accession2taxid", "ncbi_taxonomy_br.json.gz", "ncbi_lcawalkdb_br.db"] } }

def find_global_configfile():
	moduledir = os.path.dirname(os.path.abspath(__file__))
	if os.path.exists(os.path.join(moduledir, "mdmcleaner.config")) and os.path.isfile(os.path.join(moduledir, "mdmcleaner.config")):
		return os.path.join(moduledir, "mdmcleaner.config")
	raise Exception("\nError: a \"mdmcleaner.config\" file should exist under {}, but doesnt!\n".format(moduledir))

def find_local_configfile():
	cwd = os.getcwd()
	moduledir = os.path.dirname(os.path.abspath(__file__))
	if os.path.exists(os.path.join(cwd, "mdmcleaner.config")) and os.path.isfile(os.path.join(cwd, "mdmcleaner.config")):
		return os.path.join(cwd, "mdmcleaner.config")	

def read_configs(configfilelist, args): #todo: switch to config opbject instead of dictionary #todo: make it dynamically set blastdb_diamond (and others) to list of default-dbs, if no custom db is set for these
	"""
	Reads the config files in hierarchical order (first global, then local), with the later configs always overriding the previous in case of conflicts
	Config files must be tab-seperated text files (may be compressed though), with setting names in the first column, and the corresponding setting value(s) in subsequent columns.
	Setting values are always read as lists
	Unknown setting names will just be ignored. However, comments should optimally be marked with "#"
	"""
	setting_keys = ["blastp", "blastn", "diamond", "barrnap", "rnammer", "hmmsearch", "db_basedir","db_type", "blastdb_diamond", "blastdb_blastp", "blastdb_blastn", "threads"] # todo: complete this list 
	settings = {}
	for config in configfilelist:
		sys.stderr.write("\nreading settings from configfile: \"{}\"\n".format(config))
		with openfile(config) as configfile:
			for line in configfile:
				tokens = line.strip().split("#",1)[0].strip().split("\t")
				if len(tokens) <= 1:
					continue
				if tokens[0] in setting_keys:
					settings[tokens[0]] = tokens[1:]
				else:
					sys.stderr.write("\nWARNING: unknown setting key \"{}\" --> ignoring it!\n")
	if args.threads:
		settings["threads"] = args.threads
	else:
		settings["threads"] = int(settings["threads"][0])
	return settings

def getdatabase(acc2taxid_lookupfile, taxdictjson_file, lcawalk_file):#todo: find a system. either one directory per database (with the same filename for each of these) or specyfy each of these sperately
	return getdb.taxdb(acc2taxid_lookupfile, taxdictjson_file, lcawalk_file)

def check_progressdump(outfolder, infastas):
	if os.path.exists(outfolder):
		assert os.path.isdir(outfolder), "\nError: specified output-folder name cannot be used, because there already is a file with the same name!\n"
		progressfile = os.path.join(outfolder, progressdump_filename)
		return getdb.jsonfile2dict(progressfile)
	return { os.path.basename(i) : None for i in infastas} 

def write_refdb_inconsistency_report(magsag, inconsistencies, outfile):
	import io
	if len(inconsistencies) > 0:
		header = "magsag\tcontig\t{}\n".format("\t".join([key for key in inconsistencies[list(inconsistencies.keys())[0]]]))
		if isinstance(outfile, str):
			outfile = openfile(outfile, "wt")
			outfile.write(header)
	if isinstance(outfile, io.IOBase):
		for i in inconsistencies:
			line = "{}\t{}\t{}\n".format(magsag, i, "\t".join([str(v) for v in inconsistencies[i].values()]))
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

def main():
	#todo: barrnap should be skipped if rRNA.fastas already exist
	#todo: protein classifications should be saved and reloaded if existing
	#todo: protein classification should be the fast version (AND multithreaded)
	import traceback
	import time
	#todo: consider using "fromfile_prefix_chars" from ArgumentParser to optionally read arguments from file
	myparser = argparse.ArgumentParser(prog=os.path.basename(sys.argv[0]), description= "identifies and removes potential contamination in draft genomes and metagenomic bins, based on a hierarchically ranked contig classification pipeline")
	myparser.add_argument("-i", "--input_fastas", action = "store", dest = "input_fastas", nargs = "+", required = True, help = "input fastas of genomes and/or bins")  #todo: also allow genbanks
	myparser.add_argument("-o", "--output_folder", action = "store", dest = "output_folder", default = "mdmcleaner_output", help = "output-folder for MDMcleaner results. Default = 'mdmcleaner_output'")
	myparser.add_argument("-v", "--version", action = "version", version='%(prog)s {version}'.format(version=__version__))
	myparser.add_argument("--config", action = "store", dest = "configfile", default = find_local_configfile(), help = "provide a local config file with basic settings (such as the location of database-files). default: looks for config files named 'mdmcleaner.config' in current working directory. settings in the local config file will override settings in the global config file '{}'".format(os.path.join(os.path.dirname(os.path.abspath(__file__)), "mdmcleaner.config")))
	myparser.add_argument("-t", "--threads", action = "store", dest = "threads", type = int, default = None, help = "Number of threads to use. Can also be set in the  mdmcleaner.config file")
	myparser.add_argument("-f", "--force", action = "store_true", dest = "force", default = False, help = "Force reclassification of pre-existing blast-results")
	myparser.add_argument("--blast2pass", action = "store_true", dest = "blast2pass", default = "False", help = "add a second-pass blastx blast-run for all contigs without any classification on blastp level (default: False)")
	myparser.add_argument("--overview_files_basename", action = "store", dest = "overview_basename", default = "overview", help = "basename for overviewfiles (default=\"overview\"")
	args = myparser.parse_args()
	#print(args.configfile)
	
	configfile_hierarchy = [ cf for cf in [find_global_configfile(), args.configfile] if cf != None ]
	print(configfile_hierarchy)
	configs = read_configs(configfile_hierarchy, args) #todo finish this
	#initialize blastdbs
	ssu_nucblastdblist = [os.path.join(configs["db_basedir"][0], configs["db_type"][0], nbdb) for nbdb in dbfiles[configs["db_type"][0]]["ssu_nucblastdbs"]]
	lsu_nucblastdblist = [os.path.join(configs["db_basedir"][0], configs["db_type"][0], nbdb) for nbdb in dbfiles[configs["db_type"][0]]["lsu_nucblastdbs"]] #used for lsu and "tsu" rRNAs
	trna_nucblastdblist = [os.path.join(configs["db_basedir"][0], configs["db_type"][0], nbdb) for nbdb in dbfiles[configs["db_type"][0]]["trna_nucblastdbs"]]		#todo: implement tRNA_blasts		
	protblastdblist = [os.path.join(configs["db_basedir"][0], configs["db_type"][0], pbdb) for pbdb in dbfiles[configs["db_type"][0]]["protblastdbs"]]
	
	sys.stderr.write("\n\nSTETTINGS:\n" + pprint.pformat(configs)+ "\n\n")
	progressdump = check_progressdump(args.output_folder, args.input_fastas) #todo: this is meant to implement a "major-progressdump", consisting of multiple "mini-progressdumps" (one for each input-fasta). for each input-fasta, it should list the current progress-state [None = not started yet, stepxx = currently unfinished, "Finished" = finished]
	
	db = getdatabase(*[os.path.join(configs["db_basedir"][0], configs["db_type"][0], dbfile) for dbfile in dbfiles[configs["db_type"][0]]["mdmdbs"]]) #todo: instead have getdb() accept the congifs.dict as input?
	errorlistfile = openfile(args.overview_basename + "_errorlist.txt", "wt")
	# ~ hoho= getmarkers.get_trnas(args.input_fastas[0], aragorn = "aragorn", threads = 1)
	# ~ import pdb; pdb.set_trace()
	overview_before = args.overview_basename + "_all_before_cleanup.tsv"
	overview_after = args.overview_basename + "_all_after_cleanup.tsv"
	refdb_inconsistency_report = args.overview_basename + "_refdb_inconsistencies.tsv"
	for infasta in args.input_fastas:
		try:
			sys.stdout.flush()
			sys.stderr.flush()
			# ~ print("="*80)
			# ~ print(infasta)
			sys.stdout.flush()
			sys.stderr.flush()
			
			############### getting markers
			bindata = getmarkers.bindata(contigfile=infasta, threads=configs["threads"])
			nucblastjsonfilename = os.path.join(bindata.bin_resultfolder, "nucblasts.json.gz")
			protblastjsonfilename = os.path.join(bindata.bin_resultfolder, "protblasts.json")
			#import pdb; pdb.set_trace()

			############### getting blast data for markers
			##### first protein blasts

			sys.stdout.flush()
			sys.stderr.flush()			
			if os.path.isfile(protblastjsonfilename) and args.force != True: #for debugging. allows picking up AFTER blastlines were already classified when re-running
				protblasts = blasthandler.blastdata(protblastjsonfilename, score_cutoff_fraction = 0.75, continue_from_json = True, seqtype = "prot")
			else:
				print("blasting protein data")
				protblastfiles = []
				for blastdb in protblastdblist: #for protmarkers, every db is blasted one after another with full threads #todo: set list of pbdb names during initialization!
					pbdb = os.path.basename(blastdb)
					# ~ print("-"*20)
					# ~ print(configs["db_basedir"])
					# ~ print(configs["db_type"])
					# ~ print(pbdb)
					# ~ print("-"*20)
					sys.stdout.flush()
					sys.stderr.flush()	
					starttime = time.time()
					protblastfiles.append(blasthandler._run_any_blast(bindata.totalprotsfile, blastdb, "diamond", os.path.join(bindata.bin_resultfolder, "{}_totalprots_vs_{}.blast.tsv".format(bindata.bin_tempname, pbdb)), configs["threads"]))  #todo make choce of blast tool flexible. perhaps dependent on db (add tool/db tuple pairs to configs-dict)
					endtime = time.time()
					print("\nthis blast took {} seconds\n".format(endtime - starttime))
				sys.stderr.write("\nreading in protblast files...\n")
				protblasts = blasthandler.blastdata(*protblastfiles, score_cutoff_fraction = 0.75, seqtype = "prot")
				sys.stderr.write("looking up taxids of protein blast hits...\n")
				protblasts.add_info_to_blastlines(bindata, db)
				print("saving protblasts for reuse") #todo: make this optional. is only for debugging now!
				protblasts.to_json(protblastjsonfilename)
			##### then rnablasts
			if os.path.isfile(nucblastjsonfilename) and args.force != True: #for debugging. allows picking up AFTER blastlines were already classified when re-running
				nucblasts = blasthandler.blastdata(nucblastjsonfilename, score_cutoff_fraction = 0.8, continue_from_json = True, seqtype = "nuc")
			else:
				rnablastfiles = []
				starttime = time.time()
				#todo: not all rRNAs should be blasted against all databases.  so creating different blast_combiantions for lsu & ssu (and tsu and trna?) here. create a get_blast_combinations function for this in blasthander!
				# ~ nucblastdblist = [os.path.join(configs["db_basedir"][0], configs["db_type"][0], nbdb) for nbdb in dbfiles[configs["db_type"][0]]["nucblastdbs"]] #todo: set nucblastdblist during initialization!

				# ~ trna_nucblastdblist = [os.path.join(configs["db_basedir"][0], configs["db_type"][0], nbdb) for nbdb in dbfiles[configs["db_type"][0]]["trna_nucblastdbs"]]		#todo: implement tRNA_blasts		
				lsublastquerylist = [bindata.rRNA_fasta_dict["lsu_rRNA"]] #,  bindata.rRNA_fasta_dict["tsu_rRNA"]] was thinking about running the 5S blasts now, but actually they should only be run if needed at the end!
				ssublastquerylist = [bindata.rRNA_fasta_dict["ssu_rRNA"]]
				# ~ trnablastquerylist = [bindata.trnafile]
				import itertools #todo move up	
				print("blasting rRNA data") 
				#todo: the following blasts all against all (including 16S vs 23S database). But blasting 16S only makes sense against a 16S dabatase... --> ensure blasts are only against appropriate dbs![
				lsu_blast_combinations = blasthandler.get_blast_combinations(lsu_nucblastdblist, lsublastquerylist, blast = "blastn")
				ssu_blast_combinations = blasthandler.get_blast_combinations(ssu_nucblastdblist, ssublastquerylist, blast = "blastn")
			# ~ trna_blast_combinations = blasthandler.get_blast_combinations(trna_nucblastdblist, trnablastquerylist, blast = "blastn")
				all_blast_combinations = lsu_blast_combinations + ssu_blast_combinations
				# ~ print("queries")
				# ~ print(lsublastquerylist)
				# ~ print(ssublastquerylist)
				# ~ print("dbs")
				# ~ print(lsu_nucblastdblist)
				# ~ print(ssu_nucblastdblist)
				# ~ print("combinations")
				# ~ print(all_blast_combinations)
				# ~ print("---------___")
				rnablastfiles = blasthandler.run_multiple_blasts_parallel(all_blast_combinations, os.path.join(bindata.bin_resultfolder, "blastn"), configs["threads"])
				endtime = time.time()
				print("\nthis blast took {} seconds\n".format(endtime - starttime))
				nucblasts = blasthandler.blastdata(*rnablastfiles, score_cutoff_fraction = 0.8, seqtype = "nuc") #stricter cutoff for nucleotide blasts
				sys.stderr.write("looking up taxids of nucleotide blast hits...\n")
				nucblasts.add_info_to_blastlines(bindata, db)
				print("saving nuclasts for reuse") #todo: make this optional. is only for debugging now!
				nucblasts.to_json(nucblastjsonfilename)		
			# ~ import pdb; pdb.set_trace()

			############## getting LCA classifications
			if os.path.exists(bindata.pickle_progressfile):
				print("skipping LCA classification, because already present in pickle file") #todo: add check for toplevel_taxlevel set in contigdict and setting it if not
			else:
				sys.stderr.write("classifying protein sequences...\n")
				bindata.add_lca2markerdict(protblasts, db)
				sys.stderr.write("classifying rRNA sequences...\n")
				bindata.add_lca2markerdict(nucblasts, db)
				bindata.verify_arcNbac_marker(db) #todo: maybe skip that step and assume bac/arch-markers as more conserved even if assignable to the other domain? (after all, these archaeal and bacterial marker sets correspond to SINGLE-COPY markers and we don't care if they are single copy, only if they are conserved)
				#todo: combine prok with corresponding bac or arc markers for each contig
				bindata._save_current_status()

			testlca_dict_total = {}
			testlca_dict_prok = {}
			testlca_dict_23s = {}
			testlca_dict_16s = {}
			print("looping though contigs")
			# ~ import pdb; pdb.set_trace()

			################## compiling infos
			for contig in bindata.contigdict: #todo: create an own class in lca.py for this. that class should have options to filter, evaluate etc...
				# ~ print("*"*70)
				# ~ print(contig)
				# ~ print("totalprots")
				ctotalprottax = [bindata.markerdict[x]["tax"] for x in bindata.contigdict[contig]["totalprots"] if bindata.markerdict[x]["tax"] != None]
				testlca_dict_total[contig] = lca.weighted_lca(db, contig, ctotalprottax, taxlevel="total_prots_tax")
				if len(testlca_dict_total[contig]) != 0:
					bindata.contigdict[contig]["total_prots_tax"] = testlca_dict_total[contig]
				# ~ print("------\nprokmarkers")
				cprokprottax = [bindata.markerdict[x]["tax"] for x in bindata.contigdict[contig]["prok_marker"] + bindata.contigdict[contig]["bac_marker"] + bindata.contigdict[contig]["arc_marker"] if bindata.markerdict[x]["tax"] != None]
				testlca_dict_prok[contig] = lca.weighted_lca(db, contig, cprokprottax, taxlevel = "prok_marker_tax")
				if len(testlca_dict_prok[contig]) != 0:
					bindata.contigdict[contig]["prok_marker_tax"] = testlca_dict_prok[contig]
				# ~ print("------\n23S")
				c23srrnatax = [bindata.markerdict[x]["tax"] for x in bindata.contigdict[contig]["lsu_rRNA"] if bindata.markerdict[x]["tax"] != None]
				testlca_dict_23s[contig] = lca.weighted_lca(db, contig, c23srrnatax, taxlevel = "lsu_rRNA_tax")
				if len(testlca_dict_23s[contig]) != 0:
					bindata.contigdict[contig]["lsu_rRNA_tax"] = testlca_dict_23s[contig]
				# ~ print("------\n16S")
				c16srrnatax =[bindata.markerdict[x]["tax"] for x in bindata.contigdict[contig]["ssu_rRNA"] if bindata.markerdict[x]["tax"] != None]
				testlca_dict_16s[contig] = lca.weighted_lca(db, contig, c16srrnatax, taxlevel = "ssu_rRNA_tax") 
				if len(testlca_dict_16s[contig]) != 0:
					bindata.contigdict[contig]["ssu_rRNA_tax"] = testlca_dict_16s[contig]
				#todo: add another optional level, where all remaining contigs without protein hits are blasted via blastx 
				#todo: add 5S rRNA and trna LCAs from nucblasts as lowest level (only to be used if no annotation based on proteins to be found, and after second-pass blastx)
				#todo: everything without proteins should then be assumed noncoding eukaryotal (perhaps add 5S rRNA and tRNAscans first)
			# ~ import pdb; pdb.set_trace()
			bindata.get_topleveltax(db)
			# ~ import pdb; pdb.set_trace()
			#next steps:
			# list all contigs that have NO toplevel_tax but have auxiliary RNAs ("tsu"-rRNAs or trnas)
			#unclass_auxrna_genes = bindata.get_auxrna_from_unclass_contigs()
			# --> blast against lsu_blastdbs
			#auxrna_blast_combinations = blasthandler.get_blast_combinations(lsu_nucblastdblist, [unclass_auxrna_contigs], blast = "blastn")
			#auxrnablastfiles = blasthandler.run_multiple_blasts_parallel(auxrna_blast_combinations, os.path.join(bindata.bin_resultfolder, "blastn"), configs["threads"])
			#auxblasts = blasthandler.blastdata(*auxrnablastfiles, score_cutoff_fraction = 0.8)
			# use that info to get additional classifications
			#auxblasts.add_info_to_blastlines(bindata, db)
			
			#    					, those that still have no hit or have no RNAs --> blastx against protein-db
			#						, those that are still not assignable (on domain-level): mark as potential Eukaryote contamination (based on relatively high coding density of prokaryotic genomes) 
			# ~ import pdb; pdb.set_trace()
			bindata.calc_contig_scores()
			# ~ import pdb; pdb.set_trace()
			refdb_inconsistencies = bindata.doublecheck_refdb_ambig(db=db, nucblasts = nucblasts, protblasts = protblasts)
			refdb_inconsistency_report = write_refdb_inconsistency_report(bindata.bin_tempname, refdb_inconsistencies, refdb_inconsistency_report)
			
			bindata.print_contigdict(os.path.join(bindata.bin_resultfolder, "contigdict.tsv"))
			sys.stdout.flush()
			sys.stderr.flush()
			test_1(bindata, os.path.join(bindata.bin_resultfolder, "testcontigmarkersnew_beforecleanup.tsv"))
			overview_before = gather_extended_bin_metrics(bindata, outfile=overview_before, cutoff=5)
			# ~ import pdb; pdb.set_trace()
			bindata.clean_yourself() #todo: mdmcleaner runs into error if all contigs are removed fix this
			# ~ import pdb; pdb.set_trace()
			test_1(bindata, os.path.join(bindata.bin_resultfolder, "testcontigmarkersnew_aftercleanup.tsv"))
			misc.write_fasta(bindata.get_contig_records(), os.path.join(bindata.bin_resultfolder, bindata.bin_tempname + "after_cleanup.fasta.gz"))
			# ~ import pdb; pdb.set_trace()
			overview_after = gather_extended_bin_metrics(bindata, outfile=overview_after, cutoff=5)
			# ~ import pdb; pdb.set_trace()
			
		except Exception as e:
			sys.stderr.write("\nTHERE WAS A EXCEPTION WHILE HANDLING {}\n".format(infasta))
			print(e)
			traceback.print_exc()
			errorlistfile.write(infasta + "\n")
	import io
	for f in [overview_before, overview_after, refdb_inconsistency_report]:
		if isinstance(f, io.IOBase):
			f.close() #make sure whole buffer is written to files before script terminates!
	print("finished")
	print("==="*100)	

			
def test_1(bindata, outfilename):
	outfile = openfile(outfilename, "wt")
	sys.stderr.write("\n" + outfile.name + "\n")
	if len(bindata.contigdict) >= 1:
		sys.stderr.write("\nwriting results\n")
		outfile.write("contig\t{}\n".format("\t".join([x for x in bindata.contigdict[list(bindata.contigdict.keys())[0]]])))
		for contig in bindata.contigdict:
			#print("")
			#print(bindata.contigdict[contig]["totalprots"])
			#print("="*20)
			line = "{}\t{}\n".format(contig, "\t".join([";".join([str(y) for y in bindata.contigdict[contig][x]]) if type(bindata.contigdict[contig][x]) == list else str(bindata.contigdict[contig][x]) for x in bindata.contigdict[contig] ])) #todo: in protmarkerdicts change "protid" to "seqid". Add "seqid" and "marker" keys to ssu and lsu entries
			outfile.write(line)	
	else:
		sys.stderr.write("\nNo contigs left! Nothing left to write here!\n")
main()
	
