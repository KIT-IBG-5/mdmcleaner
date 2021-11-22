#!/usr/bin/env python
import sys
import os
import argparse
import misc
from misc import openfile


from _version import __version__

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

def read_configs(configfilelist, args): #todo: switch to config object instead of dictionary
	"""
	Reads the config files in hierarchical order (first global, then local), with the later configs always overriding the previous in case of conflicts
	Config files must be tab-seperated text files (may be compressed though), with setting names in the first column, and the corresponding setting value(s) in subsequent columns.
	Setting values are always read as lists
	Unknown setting names will just be ignored. However, comments should optimally be marked with "#"
	"""
	setting_keys = ["blastp", "blastn", "diamond", "barrnap", "rnammer", "hmmsearch", "aragorn", "db_basedir","db_type", "blastdb_diamond", "blastdb_blastp", "blastdb_blastn", "threads"] # todo: complete this list 
	settings = {}
	for config in configfilelist:
		sys.stderr.write("\nreading settings from configfile: \"{}\"\n".format(config))
		with openfile(config) as configfile:
			for line in configfile:
				tokens = line.strip().split("#",1)[0].strip().split()
				if len(tokens) <= 1:
					continue
				if tokens[0] in setting_keys:
					settings[tokens[0]] = tokens[1:]
				else:
					sys.stderr.write("\nWARNING: unknown setting key \"{}\" --> ignoring it!\n")
	if "threads" in vars(args):
		settings["threads"] = args.threads
	else:
		settings["threads"] = int(settings["threads"][0])
	return settings

def main():
	myparser = argparse.ArgumentParser(prog=os.path.basename(sys.argv[0]), description= "MDMcleaner pipeline v{} for decontaminating and classifying microbial dark matter MAGs and SAGs".format(__version__))
	subparsers = myparser.add_subparsers(dest = "command")
	
	cleansagmag_args = subparsers.add_parser("clean", help = "classify and filter contigs from microbial dark matter MAGs and SAGs")
	cleansagmag_args.add_argument("-i", "--input_fastas", action = "store", dest = "input_fastas", nargs = "+", required = True, help = "input fastas of genomes and/or bins")  #todo: also allow genbanks
	cleansagmag_args.add_argument("-o", "--output_folder", action = "store", dest = "output_folder", default = "mdmcleaner_output", help = "output-folder for MDMcleaner results. Default = 'mdmcleaner_output'")
	cleansagmag_args.add_argument("-v", "--version", action = "version", version='%(prog)s {version}'.format(version=__version__))
	cleansagmag_args.add_argument("-c", "--config", action = "store", dest = "configfile", default = find_local_configfile(), help = "provide a local config file with basic settings (such as the location of database-files). default: looks for config files named 'mdmcleaner.config' in current working directory. settings in the local config file will override settings in the global config file '{}'".format(os.path.join(os.path.dirname(os.path.abspath(__file__)), "mdmcleaner.config")))
	cleansagmag_args.add_argument("-t", "--threads", action = "store", dest = "threads", type = int, default = None, help = "Number of threads to use. Can also be set in the  mdmcleaner.config file")
	cleansagmag_args.add_argument("-f", "--force", action = "store_true", dest = "force", default = False, help = "Force reclassification of pre-existing blast-results")
	# ~ cleansagmag_args.add_argument("--blast2pass", action = "store_true", dest = "blast2pass", default = "False", help = "add a second-pass blastx blast-run for all contigs without any classification on blastp level (default: False)")
	cleansagmag_args.add_argument("--overview_files_basename", action = "store", dest = "overview_basename", default = "overview", help = "basename for overviewfiles (default=\"overview\"")
	cleansagmag_args.add_argument("-I", "--ignorelistfile", action = "store", dest = "ignorelistfile", default = None, help = "File listing reference-DB sequence-names that should be ignored during blast-analyses (e.g. known refDB-contaminations...")
	cleansagmag_args.add_argument("--no_filterfasta", action = "store_true", dest = "no_filterfasta", default = False, help = "Do not write filtered contigs to final output fastas (Default = False)")

	makedb_args = subparsers.add_parser("makedb", help = "Download and create MDMcleaner database")
	makedb_args.add_argument("-o", "--outdir", action = "store", dest = "outdir", default = None, help = "target base directory for reference-data. may not be the current working directory. Needs >100GB space! Default = './db/gtdb'")
	makedb_args.add_argument("-c", "--config", action = "store", dest = "configfile", default = find_local_configfile(), help = "provide a local config file with the target location to store database-files. default: looks for config files named 'mdmcleaner.config' in current working directory. settings in the local config file will override settings in the global config file '{}'".format(os.path.join(os.path.dirname(os.path.abspath(__file__)), "mdmcleaner.config")))
	makedb_args.add_argument("--verbose", action = "store_true", dest = "verbose", default = False, help = "verbose output (download progress etc)") #todo: finish implementing
	makedb_args.add_argument("--quiet", action = "store_true", dest = "quiet", default = False, help = "quiet mode (suppress any status messages except Errors and Warnings)") #todo: implement
	
	get_trna_args = subparsers.add_parser("get_trnas", help = "extracts trna sequences using Aragorn")
	get_trna_args.add_argument("infastas", nargs = "+", help = "input fasta(s). May be gzip-compressed")
	get_trna_args.add_argument("--outdir", dest = "outdir", default = ".", help = "Output directory for temporary files, etc. Default = '.'")
	get_trna_args.add_argument("-b", "--binary", default = "binary", help = "aragorn executable (with path if not in PATH). Default= assume aragorn is in PATH")
	get_trna_args.add_argument("-t", "--threads", default = 1, type = int, help = "number of parallel threads. Is only used when multiple input files are passed")
	
	get_rrnas_args = subparsers.add_parser("get_rrnas", help = "extract_rRNA sequences using barrnap")
	get_rrnas_args.add_argument("infasta", help = "input fasta. May be gzip-compressed")
	get_rrnas_args.add_argument("-o", "--outbasename", action="store", dest="outfilebasename", default = "rRNA_barrnap", help = "basename of output files (default = 'rRNA_barrnap')")
	get_rrnas_args.add_argument("-t", "--threads", action="store", dest="threads", type = int, default = 1, help = "number of threads to use (default = 1)")
	get_rrnas_args.add_argument("-b", "--binary", action="store", dest="barrnap", default = "barrnap", help = "path to barrnap binaries (if not in PATH)")
	get_rrnas_args.add_argument("--outdir", dest = "outdir", default = ".", help = "Output directory for temporary files, etc. Default = '.'")
	args = myparser.parse_args()

	if args.command in ["clean", "makedb"]:
		configfile_hierarchy = [ cf for cf in [find_global_configfile(), args.configfile] if cf != None ]
		configs = read_configs(configfile_hierarchy, args)
	
	# ~ import pdb; pdb.set_trace()
	if args.command == "clean":
		import mdmcleaner
		mdmcleaner.main(args, configs)
	
	if args.command == "makedb":
		import read_gtdb_taxonomy
		if args.outdir == None:
			assert "db_basedir" in configs, ("\n\nERROR: either 'outdir' must be specified as argument or 'db_basedir' needs to be specified in config file!\n\n")
			args.outdir = os.path.join(configs["db_basedir"][0], configs["db_type"][0])
		else:
			args.outdir = os.path.join(args.outdir, configs["db_type"][0])
		read_gtdb_taxonomy.main(args, configs)
main()
