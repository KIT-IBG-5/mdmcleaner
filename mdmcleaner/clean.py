#!/usr/bin/env python
import sys
import os
import argparse
from mdmcleaner import getdb
from mdmcleaner import misc
from mdmcleaner.misc import openfile
from mdmcleaner import getdb
from mdmcleaner import getmarkers
from mdmcleaner import blasthandler

from mdmcleaner import lca
from mdmcleaner import reporting
import itertools

from mdmcleaner._version import __version__

progressdump_filename = "mdmprogress.json"

def check_progressdump(outfolder, infastas):
	if os.path.exists(outfolder):
		assert os.path.isdir(outfolder), "\nError: specified output-folder name cannot be used, because there already is a folder with the same name!\n"
		progressfile = os.path.join(outfolder, progressdump_filename)
		return getdb.jsonfile2dict(progressfile)
	return { os.path.basename(i) : None for i in infastas} 

class FastaFileNotFoundError(Exception):
    pass

def main(args, configs):
	#todo: barrnap should be skipped if rRNA.fastas already exist
	#todo: protein classification should be the fast version (AND multithreaded)
	import traceback
	import time
	#todo: consider using "fromfile_prefix_chars" from ArgumentParser to optionally read arguments from file
	# ~ myparser = argparse.ArgumentParser(prog=os.path.basename(sys.argv[0]), description= "identifies and removes potential contamination in draft genomes and metagenomic bins, based on a hierarchically ranked contig classification pipeline")
	# ~ myparser.add_argument("-i", "--input_fastas", action = "store", dest = "input_fastas", nargs = "+", required = True, help = "input fastas of genomes and/or bins")  #todo: also allow genbanks
	# ~ myparser.add_argument("-o", "--output_folder", action = "store", dest = "output_folder", default = "mdmcleaner_output", help = "output-folder for MDMcleaner results. Default = 'mdmcleaner_output'")
	# ~ myparser.add_argument("-v", "--version", action = "version", version='%(prog)s {version}'.format(version=__version__))
	# ~ myparser.add_argument("--config", action = "store", dest = "configfile", default = find_local_configfile(), help = "provide a local config file with basic settings (such as the location of database-files). default: looks for config files named 'mdmcleaner.config' in current working directory. settings in the local config file will override settings in the global config file '{}'".format(os.path.join(os.path.dirname(os.path.abspath(__file__)), "mdmcleaner.config")))
	# ~ myparser.add_argument("-t", "--threads", action = "store", dest = "threads", type = int, default = None, help = "Number of threads to use. Can also be set in the  mdmcleaner.config file")
	# ~ myparser.add_argument("-f", "--force", action = "store_true", dest = "force", default = False, help = "Force reclassification of pre-existing blast-results")
	# ~ #myparser.add_argument("--blast2pass", action = "store_true", dest = "blast2pass", default = "False", help = "add a second-pass blastx blast-run for all contigs without any classification on blastp level (default: False)")
	# ~ myparser.add_argument("--overview_files_basename", action = "store", dest = "overview_basename", default = "overview", help = "basename for overviewfiles (default=\"overview\"")
	# ~ myparser.add_argument("--ignorelistfile", action = "store", dest = "ignorelistfile", default = None, help = "File listing reference-DB sequence-names that should be ignored during blast-analyses (e.g. known refDB-contaminations...")
	# ~ myparser.add_argument("--no_filterfasta", action = "store_true", dest = "no_filterfasta", default = False, help = "Do not write filtered contigs to final output fastas (Default = False)")
	# ~ args = myparser.parse_args()
	#print(args.configfile)
	
	# ~ configfile_hierarchy = [ cf for cf in [find_global_configfile(), args.configfile] if cf != None ]
	# ~ print(configfile_hierarchy)
	# ~ configs = read_configs(configfile_hierarchy, args) #todo finish this
	#initialize blastdbs
	if not "db_basedir" in configs.settings or len(configs.settings["db_basedir"]) == 0:
		sys.exit("\n\nERROR: No Path to reference database given! 'db_basedir' needs to be specified in the configs file!\n")
	if not "db_type" in configs.settings or len(configs.settings["db_type"]) == 0:
		sys.exit("\n\nERROR: 'db_type' needs to be specified in the configs file (even if currently 'gtdb' is still the only option...!\n")
	ssu_nucblastdblist = [os.path.join(configs.settings["db_basedir"][0], configs.settings["db_type"][0], nbdb) for nbdb in getdb.dbfiles[configs.settings["db_type"][0]]["ssu_nucblastdbs"]]
	lsu_nucblastdblist = [os.path.join(configs.settings["db_basedir"][0], configs.settings["db_type"][0], nbdb) for nbdb in getdb.dbfiles[configs.settings["db_type"][0]]["lsu_nucblastdbs"]] #used for lsu and "tsu" rRNAs
	# ~ trna_nucblastdblist = [os.path.join(configs.settings["db_basedir"][0], configs.settings["db_type"][0], nbdb) for nbdb in getdb.dbfiles[configs.settings["db_type"][0]]["genome_nucblastdbs"]]		#todo: implement tRNA_blasts		
	protblastdblist = [os.path.join(configs.settings["db_basedir"][0], configs.settings["db_type"][0], pbdb) for pbdb in getdb.dbfiles[configs.settings["db_type"][0]]["protblastdbs"]]
	
	
	#progressdump = check_progressdump(args.output_folder, args.input_fastas) #todo: this is meant to implement a "major-progressdump", consisting of multiple "mini-progressdumps" (one for each input-fasta). for each input-fasta, it should list the current progress-state [None = not started yet, stepxx = currently unfinished, "Finished" = finished]
	db = getdb.taxdb(configs)
	errorlistfile = openfile(args.overview_basename + "_errorlist.txt", "wt")

	overview_before = args.overview_basename + "_all_before_cleanup.tsv"
	overview_after = args.overview_basename + "_all_after_cleanup.tsv"
	refdb_ambiguity_report = args.overview_basename + "_refdb_ambiguities.tsv"
	db_suspects = None
	for infasta in args.input_fastas:
		try:
			if not os.path.exists(infasta) or not os.path.isfile(infasta):
				raise FastaFileNotFoundError
				
			############### getting markers
			bindata = getmarkers.bindata(contigfile=infasta, outbasedir=args.output_folder, configs = configs)
			nucblastjsonfilename = os.path.join(bindata.bin_resultfolder, "nucblasts.json.gz")
			protblastjsonfilename = os.path.join(bindata.bin_resultfolder, "protblasts.json.gz")

			############### getting blast data for markers
			##### first protein blasts

			sys.stdout.flush()
			sys.stderr.flush()	
			if os.path.isfile(protblastjsonfilename) and args.force != True: #for debugging. allows picking up AFTER blastlines were already classified when re-running
				sys.stderr.write("\n-->using preexisting protein blast results in resultfolder!\n")
				protblasts = blasthandler.blastdata(protblastjsonfilename, score_cutoff_fraction = 0.75, continue_from_json = True, seqtype = "prot", blacklist = configs.blacklist)
			else:
				sys.stderr.write("\n-->blasting protein data\n")
				protblastfiles = []
				for blastdb in protblastdblist: #for protmarkers, every db is blasted one after another with full threads #todo: set list of pbdb names during initialization!
					pbdb = os.path.basename(blastdb)
					sys.stdout.flush()
					sys.stderr.flush()	
					# ~ starttime = time.time()
					protblastfiles.append(blasthandler._run_any_blast(bindata.totalprotsfile, blastdb, configs.settings["diamond"], outname=os.path.join(bindata.bin_resultfolder, "{}_totalprots_vs_{}.blast.tsv".format(bindata.bin_tempname, pbdb)), threads=configs.settings["threads"]))  #todo make choce of blast tool flexible. perhaps dependent on db (add tool/db tuple pairs to configs-dict)
					# ~ endtime = time.time()
					# ~ print("\nthis blast took {} seconds\n".format(endtime - starttime))
				sys.stderr.write("\treading in protblast files...\n")
				protblasts = blasthandler.blastdata(*protblastfiles, score_cutoff_fraction = 0.75, seqtype = "prot", blacklist = configs.blacklist)
				sys.stderr.write("\tlooking up taxids of protein blast hits...\n")
				protblasts.add_info_to_blastlines(bindata, db)
				sys.stderr.write("\tsaving protblasts for reuse\n") #todo: make this optional. is only for debugging now!
				protblasts.to_json(protblastjsonfilename)
			##### then rnablasts
			if os.path.isfile(nucblastjsonfilename) and args.force != True: #for debugging. allows picking up AFTER blastlines were already classified when re-running
				sys.stderr.write("\n-->using preexisting ribosomal rRNA blast results in resultfolder!\n")
				nucblasts = blasthandler.blastdata(nucblastjsonfilename, score_cutoff_fraction = 0.8, continue_from_json = True, seqtype = "nuc", blacklist = configs.blacklist)
			else:
				rnablastfiles = []
				# ~ starttime = time.time()
				#todo: not all rRNAs should be blasted against all databases.  so creating different blast_combiantions for lsu & ssu (and tsu and trna?) here. create a get_blast_combinations function for this in blasthander!
				# ~ nucblastdblist = [os.path.join(configs.settings["db_basedir"][0], configs.settings["db_type"][0], nbdb) for nbdb in getdb.dbfiles[configs.settings["db_type"][0]]["nucblastdbs"]] #todo: set nucblastdblist during initialization!

				# ~ trna_nucblastdblist = [os.path.join(configs.settings["db_basedir"][0], configs.settings["db_type"][0], nbdb) for nbdb in getdb.dbfiles[configs.settings["db_type"][0]]["trna_nucblastdbs"]]		#todo: implement tRNA_blasts		
				lsublastquerylist = [bindata.rRNA_fasta_dict["lsu_rRNA"]] #,  bindata.rRNA_fasta_dict["tsu_rRNA"]] was thinking about running the 5S blasts now, but actually they should only be run if needed at the end!
				ssublastquerylist = [bindata.rRNA_fasta_dict["ssu_rRNA"]]
				# ~ trnablastquerylist = [bindata.trnafile]
	
				sys.stderr.write("\n-->blasting rRNA data\n") 
				#todo: the following blasts all against all (including 16S vs 23S database). But blasting 16S only makes sense against a 16S dabatase... --> ensure blasts are only against appropriate dbs![
				lsu_blast_combinations = blasthandler.get_blast_combinations(lsu_nucblastdblist, lsublastquerylist, blast = configs.settings["blastn"])
				ssu_blast_combinations = blasthandler.get_blast_combinations(ssu_nucblastdblist, ssublastquerylist, blast = configs.settings["blastn"])
			# ~ trna_blast_combinations = blasthandler.get_blast_combinations(trna_nucblastdblist, trnablastquerylist, blast = "blastn")
				all_blast_combinations = lsu_blast_combinations + ssu_blast_combinations
				rnablastfiles = blasthandler.run_multiple_blasts_parallel(all_blast_combinations, outbasename=os.path.join(bindata.bin_resultfolder, "blastn"), total_threads=configs.settings["threads"])
				# ~ endtime = time.time()
				# ~ print("\nthis blast took {} seconds\n".format(endtime - starttime))
				nucblasts = blasthandler.blastdata(*rnablastfiles, score_cutoff_fraction = 0.8, seqtype = "nuc", blacklist = configs.blacklist) #stricter cutoff for nucleotide blasts
				sys.stderr.write("\tlooking up taxids of nucleotide blast hits...\n")
				nucblasts.add_info_to_blastlines(bindata, db)
				sys.stderr.write("\tsaving nuclasts for reuse\n") #todo: make this optional. is only for debugging now!
				nucblasts.to_json(nucblastjsonfilename)		


			############## getting LCA classifications
			sys.stderr.write("\n-->LCA classifications\n")
			if os.path.exists(bindata.pickle_progressfile):
				sys.stderr.write("\tskipping LCA classification, because already present in pickle file\n") #todo: add check for toplevel_taxlevel set in contigdict and setting it if not
			else:
				sys.stderr.write("\tclassifying protein sequences...\n")
				bindata.add_lca2markerdict(protblasts, db)
				sys.stderr.write("\tclassifying rRNA sequences...\n")
				bindata.add_lca2markerdict(nucblasts, db)
				bindata.verify_arcNbac_marker(db) #todo: maybe skip that step and assume bac/arch-markers as more conserved even if assignable to the other domain? (after all, these archaeal and bacterial marker sets correspond to SINGLE-COPY markers and we don't care if they are single copy, only if they are conserved)
				#todo: combine prok with corresponding bac or arc markers for each contig
				bindata._save_current_status()

			testlca_dict_total = {}
			testlca_dict_prok = {}
			testlca_dict_23s = {}
			testlca_dict_16s = {}

			################## compiling infos #todo: condense this!
			for contig in bindata.contigdict: #todo: create an own class in lca.py for this. that class should have options to filter, evaluate etc...
				ctotalprottax = [bindata.markerdict[x]["tax"] for x in bindata.contigdict[contig]["totalprots"] if bindata.markerdict[x]["tax"] != None]
				testlca_dict_total[contig] = lca.weighted_lca(db, contig, ctotalprottax, taxlevel="totalprots_tax")
				if len(testlca_dict_total[contig]) != 0:
					bindata.contigdict[contig]["totalprots_tax"] = testlca_dict_total[contig]

				cprokprottax = [bindata.markerdict[x]["tax"] for x in bindata.contigdict[contig]["prok_marker"] + bindata.contigdict[contig]["bac_marker"] + bindata.contigdict[contig]["arc_marker"] if bindata.markerdict[x]["tax"] != None]
				testlca_dict_prok[contig] = lca.weighted_lca(db, contig, cprokprottax, taxlevel = "prok_marker_tax")
				if len(testlca_dict_prok[contig]) != 0:
					bindata.contigdict[contig]["prok_marker_tax"] = testlca_dict_prok[contig]

				c23srrnatax = [bindata.markerdict[x]["tax"] for x in bindata.contigdict[contig]["lsu_rRNA"] if bindata.markerdict[x]["tax"] != None]
				testlca_dict_23s[contig] = lca.weighted_lca(db, contig, c23srrnatax, taxlevel = "lsu_rRNA_tax")
				if len(testlca_dict_23s[contig]) != 0:
					bindata.contigdict[contig]["lsu_rRNA_tax"] = testlca_dict_23s[contig]

				c16srrnatax =[bindata.markerdict[x]["tax"] for x in bindata.contigdict[contig]["ssu_rRNA"] if bindata.markerdict[x]["tax"] != None]
				testlca_dict_16s[contig] = lca.weighted_lca(db, contig, c16srrnatax, taxlevel = "ssu_rRNA_tax") 
				if len(testlca_dict_16s[contig]) != 0:
					bindata.contigdict[contig]["ssu_rRNA_tax"] = testlca_dict_16s[contig]
				#todo: add another optional level, where all remaining contigs without protein hits are blasted via blastx 
				#todo: add 5S rRNA and trna LCAs from nucblasts as lowest level (only to be used if no annotation based on proteins to be found, and after second-pass blastx)
				#todo: everything without proteins should then be assumed noncoding eukaryotal (perhaps add 5S rRNA and tRNAscans first)
				#todo: instead of using eukaryotic proteins, do a bacterial stile ORF prediction on eukaryotic genome fragments, use THOSE as protein referencces and all fragments WITHOUT a CDS as nucleotide references!

			bindata.get_topleveltax(db)
			# ~ import pdb; pdb.set_trace()
			#next steps:
			# list all contigs that have NO toplevel_tax but have auxiliary RNAs ("tsu"-rRNAs or trnas)
			#unclass_auxrna_genes = bindata.get_auxrna_from_unclass_contigs()
			# --> blast against lsu_blastdbs
			#auxrna_blast_combinations = blasthandler.get_blast_combinations(lsu_nucblastdblist, [unclass_auxrna_contigs], blast = "blastn")
			#auxrnablastfiles = blasthandler.run_multiple_blasts_parallel(auxrna_blast_combinations, os.path.join(bindata.bin_resultfolder, "blastn"), configs.settings["threads"])
			#auxblasts = blasthandler.blastdata(*auxrnablastfiles, score_cutoff_fraction = 0.8)
			# use that info to get additional classifications
			#auxblasts.add_info_to_blastlines(bindata, db)
			
			#    					, those that still have no hit or have no RNAs --> blastx against protein-db
			#						, those that are still not assignable (on domain-level): mark as potential Eukaryote contamination (based on relatively high coding density of prokaryotic genomes) 
			# ~ print("db_suspects = {}".format(db_suspects))
			db_suspects = bindata.evaluate_and_flag_all_contigs(db=db, protblasts=protblasts, nucblasts=nucblasts, db_suspects=db_suspects, fast_run=args.fast_run)
			# ~ print("db_suspects_tempdir = {}".format(db_suspects.tempdir))
			sys.stderr.write("\n--> writing to output files\n")
			
			refdb_ambiguity_report = reporting.write_refdb_ambiguity_report(bindata.bin_tempname, bindata.ref_db_ambiguity_overview, refdb_ambiguity_report)
			# ~ bindata.print_contigdict(os.path.join(bindata.bin_resultfolder, "contigdict.tsv"))
			
			sys.stderr.flush()

			reporting.write_full_bindata(bindata, os.path.join(bindata.bin_resultfolder, "fullcontiginfos_beforecleanup.tsv"))
			overview_before = reporting.gather_extended_bin_metrics(bindata, outfile=overview_before, cutoff=5)
			if not args.no_filterfasta:
				bindata.sort_and_write_contigs()
			bindata.write_krona_inputtable(db)
			
			# ~ bindata.clean_yourself() #todo: mdmcleaner runs into error if all contigs are removed fix this
			# ~ test_1(bindata, os.path.join(bindata.bin_resultfolder, "testcontigmarkersnew_aftercleanup.tsv"))
			# ~ misc.write_fasta(bindata.get_contig_records(), os.path.join(bindata.bin_resultfolder, bindata.bin_tempname + "after_cleanup.fasta.gz"))
			# ~ overview_after = gather_extended_bin_metrics(bindata, outfile=overview_after, cutoff=5)
			if not args.fast_run:
				sys.stderr.write("reference-database contaminations detected during this run: {}".format(len(db_suspects.blacklist_additions)))
			
		except FastaFileNotFoundError:
			sys.stderr.write("\nERROR: Input File {} does not exist! --> skipping it!\n".format(infasta))
			errorlistfile.write(infasta + "\n")
		
		except Exception as e:
			sys.stderr.write("\nTHERE WAS A EXCEPTION WHILE HANDLING {}\n".format(infasta))
			sys.stderr.write("\n{}\n".format(e))
			traceback.print_exc()
			errorlistfile.write(infasta + "\n")

	if not args.fast_run and db_suspects != None:
		if "contamination" in db_suspects.collective_diamondblast():
			sys.stderr.write("\nWARNING: potential eukaryotic contaminants were determined in reference genomes. some classifications may need to be ajdjusted.\nIt is recommended to run the pipeline again (as is), with the updated blacklist to correct that (most intermediate results can be reused, so this will be faster than the original run)\n\n")
		if len(db_suspects.blacklist_additions) > 0:
			sys.stderr.write("\nA total of {} reference-database entries were newly detected as contaminants! Please note the updated blacklist!\n".format(len(db_suspects.blacklist_additions)))
		
	import io
	for f in [overview_before, overview_after, refdb_ambiguity_report]:
		if isinstance(f, io.IOBase):
			f.close() #make sure whole buffer is written to files before script terminates!
	sys.stderr.write("\n{line}finished{line}\n".format(line="-"*30))

	if db_suspects != None: #removing tempfiles from refdb-ambiguity assessments
		for f in db_suspects.delete_filelist:
			os.remove(f)
		os.rmdir(db_suspects.tempdir)
		return db_suspects.blacklist_additions

if __name__ == '__main__':
	main()		

#todo: add more argparse options and multiple parsers
	
