#!/usr/bin/env python
import sys
import os
import argparse
import getdb
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
						"mdmdbs" : ["gtdb_all.accession2taxid.sorted", "gtdb_taxonomy_br.json.gz", "gtdb_lcawalkdb_br.db"] },\
			"ncbi" : { "protblastdbs" : ["nr"], \
						"nucblastdbs" : ["nt"], \
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

def main():
	import time
	#todo: consider using "fromfile_prefix_chars" from ArgumentParser to optionally read arguments from file
	myparser = argparse.ArgumentParser(prog=os.path.basename(sys.argv[0]), description= "identifies and removes potential contamination in draft genomes and metagenomic bins, based on a hierarchically ranked contig classification pipeline")
	myparser.add_argument("-i", "--input_fastas", action = "store", dest = "input_fastas", nargs = "+", required = True, help = "input fastas of genomes and/or bins")  #todo: also allow genbanks
	myparser.add_argument("-o", "--output_folder", action = "store", dest = "output_folder", default = "mdmcleaner_output", help = "output-folder for MDMcleaner results. Default = 'mdmcleaner_output'")
	myparser.add_argument("-v", "--version", action = "version", version='%(prog)s {version}'.format(version=__version__))
	myparser.add_argument("--config", action = "store", dest = "configfile", default = find_local_configfile(), help = "provide a local config file with basic settings (such as the location of database-files). default: looks for config files named 'mdmcleaner.config' in current working directory. settings in the local config file will override settings in the global config file '{}'".format(os.path.join(os.path.dirname(os.path.abspath(__file__)), "mdmcleaner.config")))
	myparser.add_argument("-t", "--threads", action = "store", dest = "threads", type = int, default = None, help = "Number of threads to use. Can also be set in the  mdmcleaner.config file")
	args = myparser.parse_args()
	#print(args.configfile)
	
	
	configfile_hierarchy = [ cf for cf in [find_global_configfile(), args.configfile] if cf != None ]
	print(configfile_hierarchy)
	configs = read_configs(configfile_hierarchy, args) #todo finish this
	sys.stderr.write("\n\nSTETTINGS:\n" + pprint.pformat(configs)+ "\n\n")
	progressdump = check_progressdump(args.output_folder, args.input_fastas) #todo: this is meant to implement a "major-progressdump", consisting of multiple "mini-progressdumps" (one for each input-fasta). for each input-fasta, it should list the current progress-state [None = not started yet, stepxx = currently unfinished, "Finished" = finished]
	
	db = getdatabase(*[os.path.join(configs["db_basedir"][0], configs["db_type"][0], dbfile) for dbfile in dbfiles[configs["db_type"][0]]["mdmdbs"]]) #todo: instead have getdb() accept the congifs.dict as input?
	for infasta in args.input_fastas:
		#get markers
		bindata = getmarkers.bindata(contigfile=infasta, threads=configs["threads"])
		#import pdb; pdb.set_trace()
		test_1(bindata)
		#blast markers
		protblastfiles = []
		for pbdb in dbfiles[configs["db_type"][0]]["protblastdbs"]: #for protmarkers, every db is blasted one after another with full threads #todo: set list of pbdb names during initialization!
			print("-"*20)
			print(configs["db_basedir"])
			print(configs["db_type"])
			print(pbdb)
			print("-"*20)
			blastdb = os.path.join(configs["db_basedir"][0], configs["db_type"][0], pbdb)
			starttime = time.time()
			protblastfiles.append(blasthandler._run_any_blast(bindata.totalprotfile, blastdb, "diamond", os.path.join(bindata.bin_resultfolder, "{}_totalprots_vs_{}.blast.tsv".format(bindata.bin_tempname, pbdb)), configs["threads"]))  #todo make choce of blast tool flexible. perhaps dependent on db (add tool/db tuple pairs to configs-dict)
			endtime = time.time()
			print("\nthis blast took {} seconds\n".format(endtime - starttime))
		rnablastfiles = []
		starttime = time.time()
		nucblastdblist = [os.path.join(configs["db_basedir"][0], configs["db_type"][0], nbdb) for nbdb in dbfiles[configs["db_type"][0]]["nucblastdbs"]] #todo: set nucblastdblist during initialization!
		nucblastquerylist = list(bindata.rRNA_fasta_dict.values())
		import itertools #todo move up	
		print("blasting rRNA data") 
		#todo: the following blasts all against all (including 16S vs 23S database). But blasting 16S only makes sense against a 16S dabatase... --> ensure blasts are only against appropriate dbs![
		all_blast_combinations = [ blasttuple + ("blastn",) for blasttuple in list(itertools.chain(*list(zip(nucblastquerylist, permu) for permu in itertools.permutations(nucblastdblist, len(nucblastquerylist)))))] #Todo see if this works correctly. Only works as long as nucblastdbist is longer or equal to nucblastquerylist...
		print(all_blast_combinations)
		rnablastfiles = blasthandler.run_multiple_blasts_parallel(all_blast_combinations, os.path.join(bindata.bin_resultfolder, "blastn"), configs["threads"])
		endtime = time.time()
		print("\nthis blast took {} seconds\n".format(endtime - starttime))
		#import pdb; pdb.set_trace()
		sys.stderr.write("\nreading in blast files...\n")
		protblasts = blasthandler.blastdata(*protblastfiles, score_cutoff_fraction = 0.75)
		# ~ print("="*50)
		# ~ print("NUCBLASTS")
		nucblasts = blasthandler.blastdata(*rnablastfiles, score_cutoff_fraction = 0.8) #stricter cutoff for nucleotide blasts
		# ~ import pdb; pdb.set_trace()
		sys.stderr.write("looking up taxids of protein blast hits...\n")
		protblasts.add_info_to_blastlines(bindata, db)
		sys.stderr.write("looking up taxids of nucleotide blast hits...\n")
		nucblasts.add_info_to_blastlines(bindata, db)
		sys.stderr.write("classifying protein sequences...\n")
		bindata.add_lca2markerdict(protblasts, db)
		sys.stderr.write("classifying rRNA sequences...\n")
		bindata.add_lca2markerdict(nucblasts, db)
		bindata.verify_arcNbac_marker(db) #todo: maybe skip that step and assume bac/arch-markers as more conserved even if assignable to the other domain? (after all, these archaeal and bacterial marker sets correspond to SINGLE-COPY markers and we don't care if they are single copy, only if they are conserved)
		#todo: combine prok with corresponding bac or arc markers for each contig
		testlca_dict_total = {}
		testlca_dict_prok = {}
		testlca_dict_23s = {}
		testlca_dict_16s = {}
		print("looping though contigs")
		for contig in bindata.contigdict: #todo: create an own class in lca.py for this. that class should have options to filter, evaluate etc...
			# ~ print("*"*70)
			print(contig)
			# ~ print("totalprots")
			ctotalprottax = [bindata.markerdict[x]["tax"] for x in bindata.contigdict[contig]["totalprots"] if bindata.markerdict[x]["tax"] != None]
			testlca_dict_total[contig] = lca.weighted_lca(db, contig, ctotalprottax)
			if len(testlca_dict_total[contig]) != 0:
				bindata.contigdict[contig]["total_prots_tax"] = testlca_dict_total[contig]
			# ~ print("------\nprokmarkers")
			cprokprottax = [bindata.markerdict[x]["tax"] for x in bindata.contigdict[contig]["prok_marker"] + bindata.contigdict[contig]["bac_marker"] + bindata.contigdict[contig]["arc_marker"] if bindata.markerdict[x]["tax"] != None]
			testlca_dict_prok[contig] = lca.weighted_lca(db, contig, cprokprottax)
			if len(testlca_dict_prok[contig]) != 0:
				bindata.contigdict[contig]["prok_marker_tax"] = testlca_dict_prok[contig]
			# ~ print("------\n23S")
			c23srrnatax = [bindata.markerdict[x]["tax"] for x in bindata.contigdict[contig]["lsu_rRNA"] if bindata.markerdict[x]["tax"] != None]
			testlca_dict_23s[contig] = lca.weighted_lca(db, contig, c23srrnatax)
			if len(testlca_dict_23s[contig]) != 0:
				bindata.contigdict[contig]["lsu_rRNA_tax"] = testlca_dict_23s[contig]
			# ~ print("------\n16S")
			c16srrnatax =[bindata.markerdict[x]["tax"] for x in bindata.contigdict[contig]["ssu_rRNA"] if bindata.markerdict[x]["tax"] != None]
			testlca_dict_16s[contig] = lca.weighted_lca(db, contig, c16srrnatax) 
			if len(testlca_dict_16s[contig]) != 0:
				bindata.contigdict[contig]["ssu_rRNA_tax"] = testlca_dict_16s[contig]
			#todo: add another optional level, where all remaining contigs without protein hits are blasted via blastx 
		# ~ import pdb; pdb.set_trace()
		bindata.get_major_taxon(db)
		import pdb; pdb.set_trace()
		bindata.calc_contig_scores()
		import pdb; pdb.set_trace()
		bindata.print_contigdict("contigdict.tsv")
	print("finished")

				
			

def test_1(bindata):
	outfile = openfile(os.path.join(bindata.bin_resultfolder, "testcontigmarkers.tsv"), "wt")
	sys.stderr.write("\n" + outfile.name + "\n")
	sys.stderr.write("\nwriting results\n")
	outfile.write("contig\t{}\n".format("\t".join([x for x in bindata.contigdict[list(bindata.contigdict.keys())[0]]])))
	for contig in bindata.contigdict:
		#print("")
		#print(bindata.contigdict[contig]["totalprots"])
		#print("="*20)
		line = "{}\t{}\n".format(contig, "\t".join([";".join([str(y) for y in bindata.contigdict[contig][x]]) if type(bindata.contigdict[contig][x]) == list else str(bindata.contigdict[contig][x]) for x in bindata.contigdict[contig] ])) #todo: in protmarkerdicts change "protid" to "seqid". Add "seqid" and "marker" keys to ssu and lsu entries
		outfile.write(line)	

main()
	
