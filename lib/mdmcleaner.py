#!/usr/bin/env python
import sys
import os
import argparse
import misc
from misc import openfile
import pprint

setting_keys = ["blacklistfile","blastp", "blastn", "blastdbcmd", "diamond", "barrnap", "rnammer", "hmmsearch", "aragorn", "db_basedir","db_type", "blastdb_diamond", "blastdb_blastp", "blastdb_blastn", "threads"] # todo: complete this list 
from _version import __version__

def write_blacklist(blacklist, outfilename):
	with misc.openfile(outfilename, "at") as outfile:
		for b in blacklist:
			outfile.write("{}\n".format(b))

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
	import os
	settings, settings_source = {}, {}
	settings["blacklistfile"] = []
	if "ignore_default_blacklist" in vars(args) and args.ignore_default_blacklist == True:
		settings_source["blacklistfile"] = None
	else:
		mdmcleaner_lib_path = os.path.dirname(os.path.realpath(__file__))
		default_blacklist = os.path.join(mdmcleaner_lib_path, "blacklist.list")
		settings["blacklistfile"].append(default_blacklist)
		settings_source["blacklistfile"] = "default"
	for config in configfilelist: #todo: maybe instead of replacing blacklist with that in config-file, just append it to the list of blacklists?
		sys.stderr.write("\nreading settings from configfile: \"{}\"\n".format(config))
		with openfile(config) as configfile:
			for line in configfile:
				tokens = line.strip().split("#",1)[0].strip().split()
				if len(tokens) <= 1:
					continue
				if tokens[0] in setting_keys:
					settings[tokens[0]] = tokens[1:]
					settings_source[tokens[0]] = config
				else:
					sys.stderr.write("\nWARNING: unknown setting key \"{}\" --> ignoring it!\n")
	if "threads" in vars(args):
		settings["threads"] = args.threads
	else:
		settings["threads"] = int(settings["threads"][0])
	if "blacklistfile" in vars(args) and args.blacklistfile != None:
		settings["blacklistfile"].append(args.blacklistfile)
	return settings, settings_source

def _read_blacklistfiles(blacklistfilelist):
	blacklist = []
	for blacklistfile in blacklistfilelist:
		with openfile(blacklistfile) as blf:
			for line in blf:
				tokens = line.strip().split("#")
				if len(tokens) >0 and tokens[0] != "":
					blacklist.append(tokens[0])
	return set(blacklist)

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
	cleansagmag_args.add_argument("-b", "--blacklistfile", action = "store", dest = "blacklistfile", default = None, help = "File listing reference-DB sequence-names that should be ignored during blast-analyses (e.g. known refDB-contaminations...")
	cleansagmag_args.add_argument("--no_filterfasta", action = "store_true", dest = "no_filterfasta", default = False, help = "Do not write filtered contigs to final output fastas (Default = False)")
	cleansagmag_args.add_argument("--ignore_default_blacklist", action = "store_true", dest = "ignore_default_blacklist", default = False, help = "Ignore the default blacklist (Default = False)")


	makedb_args = subparsers.add_parser("makedb", help = "Download and create MDMcleaner database")
	makedb_args.add_argument("-o", "--outdir", action = "store", dest = "outdir", default = None, help = "target base directory for reference-data. may not be the current working directory. Needs >100GB space! Default = './db/gtdb'")
	makedb_args.add_argument("-c", "--config", action = "store", dest = "configfile", default = find_local_configfile(), help = "provide a local config file with the target location to store database-files. default: looks for config files named 'mdmcleaner.config' in current working directory. settings in the local config file will override settings in the global config file '{}'".format(os.path.join(os.path.dirname(os.path.abspath(__file__)), "mdmcleaner.config")))
	makedb_args.add_argument("--verbose", action = "store_true", dest = "verbose", default = False, help = "verbose output (download progress etc)") #todo: finish implementing
	makedb_args.add_argument("--quiet", action = "store_true", dest = "quiet", default = False, help = "quiet mode (suppress any status messages except Errors and Warnings)") #todo: implement
	
	get_marker_args = subparsers.add_parser("get_markers", help = "extracts protein coding and/or rRNA gene sequences from input genome(s)")
	get_marker_args.add_argument("-i", "--input_fastas", action = "store", dest = "input_fastas", nargs = "+", help = "input fasta(s). May be gzip-compressed")
	get_marker_args.add_argument("-m", "--markertype", action = "store", dest = "markertype", default = "all", choices = ["rrna", "trna", "totalprots", "markerprots", "all"], help = "type of marker gene that should be extracted (default = 'all')")
	get_marker_args.add_argument("-o", "--outdir", dest = "outdir", default = ".", help = "Output directory (will be created if it does not exist). Default = '.'")
	get_marker_args.add_argument("-c", "--config", action = "store", dest = "configfile", default = find_local_configfile(), help = "provide a local config file with the target location to store database-files. default: looks for config files named 'mdmcleaner.config' in current working directory. settings in the local config file will override settings in the global config file '{}'".format(os.path.join(os.path.dirname(os.path.abspath(__file__)), "mdmcleaner.config")))
	get_marker_args.add_argument("-t", "--threads", action="store", dest="threads", type = int, default = 1, help = "number of threads to use (default = 1)")
	get_marker_args.add_argument("-M", "--mincontiglngth", action="store", dest="mincontiglength", type = int, default = 0, help = "minimum contig length (contigs shorter than this will be ignored)")	

	acc2taxpath = subparsers.add_parser("acc2taxpath", help = "Get full taxonomic path assorciated with a specific acession number")
	acc2taxpath.add_argument("accessions", action = "store", nargs = "+", help = "(space seperated list of) input accessions. Or just pass \"interactive\" for interactive mode")
	acc2taxpath.add_argument("-c", "--config", action = "store", dest = "configfile", default = find_local_configfile(), help = "provide a local config file with the target location to store database-files. default: looks for config files named 'mdmcleaner.config' in current working directory. settings in the local config file will override settings in the global config file '{}'".format(os.path.join(os.path.dirname(os.path.abspath(__file__)), "mdmcleaner.config")))

	evaluate_refdbcontam_args = subparsers.add_parser("refdb_contams", help = "EXPERIMENTAL: evaluate potentiel refDB-contaminations")
	evaluate_refdbcontam_args.add_argument("ambiguity_report", nargs = "?", default = "./overview_refdb_ambiguities.tsv")
	evaluate_refdbcontam_args.add_argument("-c", "--config", action = "store", dest = "configfile", default = find_local_configfile(), help = "provide a local config file with the target location to store database-files. default: looks for config files named 'mdmcleaner.config' in current working directory. settings in the local config file will override settings in the global config file '{}'".format(os.path.join(os.path.dirname(os.path.abspath(__file__)), "mdmcleaner.config")))
	evaluate_refdbcontam_args.add_argument("-t", "--threads", action="store", dest="threads", type = int, default = 1, help = "number of threads to use (default = 1)")
	evaluate_refdbcontam_args.add_argument("-o", "--outblacklist", action="store", dest="outblacklist", default = "new_blacklist_additions.tsv", help = "Outputfile for new blacklist additions. If a preexisting file is selected, additions will be appended to end of that file")
	evaluate_refdbcontam_args.add_argument("--tempbasename", action="store", dest="tempbasename", default = "mdmtempfile", help = "basename for temporary files (default = 'mdmtempfile'")
	
	
	set_configs_args = subparsers.add_parser("set_configs", help = "setting or changing settings in config files")
	set_configs_args.add_argument("-s", "--scope", action = "store", dest = "scope", choices = ["local", "global"], default = "local", help = "change settings in local or global config file. 'global' likely require admin privileges. 'local' will modify or create a mdmcleaner.config file in the current working directory. default = 'local'")
	set_configs_args.add_argument("--blastp", action = "store", dest = "blastp", help = "path to blastp binaries (if not in PATH)")
	set_configs_args.add_argument("--blastn", action = "store", dest = "blastn", help = "path to blastn binaries (if not in PATH)")
	set_configs_args.add_argument("--diamond", action = "store", dest = "diamond", help = "path to diamond binaries (if not in PATH)")
	set_configs_args.add_argument("--barrnap", action = "store", dest = "barrnap", help = "path to barrnap binaries (if not in PATH)")
	set_configs_args.add_argument("--hmmsearch", action = "store", dest = "hmmsearch", help = "path to hmmsearch binaries (if not in PATH)")
	set_configs_args.add_argument("--aragorn", action = "store", dest = "aragorn", help = "path to aragorn binaries (if not in PATH)")
	set_configs_args.add_argument("--db_basedir", action = "store", dest = "db_basedir", help = "path to basedirectory for reference database")
	# ~ set_configs_args.add_argument("--db_type", action = "store", dest = "db_type", choices = ["gtdb"], help = "taxonomic system used (currently only gtdb available)") #ncbi nor re-implemented yet. no actual choice right now
	set_configs_args.add_argument("--threads", action = "store", dest = "threads", help = "threads to use by default")

	show_configs_args = subparsers.add_parser("show_configs", help = "check settings in config files (with highest ranking source for each setting)")
	show_configs_args.add_argument("-c", "--config", action = "store", dest = "configfile", default = find_local_configfile(), help = "local config file (default: search for config file in current working directory, use None if not present)")

	version_args = subparsers.add_parser("version", help = "show version info and exit")
	if len(sys.argv)==1:
		myparser.print_help(sys.stderr)
		sys.exit(1)
	args = myparser.parse_args()
	
	if args.command == "version":
		sys.stderr.write("MDMcleaner v{}\n".format(__version__))

	if args.command in ["clean", "makedb", "show_configs", "get_markers", "acc2taxpath", "refdb_contams"]:
		configfile_hierarchy = [ cf for cf in [find_global_configfile(), args.configfile] if cf != None ]
		configs, settings_source = read_configs(configfile_hierarchy, args)
		sys.stderr.write("\n\nSETTINGS:\n" + pprint.pformat(configs)+ "\n\n")
		if args.command in ["clean", "get_markers", "refdb_contams"]:
			configs["blacklist"] =  _read_blacklistfiles(configs["blacklistfile"])
	
	# ~ import pdb; pdb.set_trace()
	if args.command == "clean":
		import clean
		clean.main(args, configs)
	
	if args.command == "makedb":
		import read_gtdb_taxonomy
		if args.outdir == None:
			assert "db_basedir" in configs, ("\n\nERROR: either 'outdir' must be specified as argument or 'db_basedir' needs to be specified in config file!\n\n")
			args.outdir = os.path.join(configs["db_basedir"][0], configs["db_type"][0])
		else:
			args.outdir = os.path.join(args.outdir, configs["db_type"][0])
		read_gtdb_taxonomy.main(args, configs)
		
	# ~ if args.command == "get_rrnas":
		# ~ import getmarkers
		# ~ getmarkers.get_rrnas_only(args)
		
	# ~ if args.command == "get_rrnas":
		# ~ pass
		
	if args.command == "get_markers":
		import getmarkers
		getmarkers.get_only_marker_seqs(args, configs)

	if args.command == "set_configs":
		# ~ import pdb; pdb.set_trace()
		if args.scope == "global":
			configfile_hierarchy = [find_global_configfile()]
			outfile = misc.openfile(configfile_hierarchy[0], "wt")
		else:
			configfile_hierarchy = [ cf for cf in [find_global_configfile(), find_local_configfile()] if cf != None ]
		configs, settings_source = read_configs(configfile_hierarchy, args)
		if args.scope == "local":
			outfile = misc.openfile("mdmcleaner.config", "wt")
		arg_settings = vars(args)
		for a in arg_settings:
			if a != "scope" and arg_settings[a] != None:
				configs[a] = [arg_settings[a]]
		for c in configs:
			if c in setting_keys and configs[c]:
				outfile.write("{}\t{}\n".format(c, configs[c][0]))
		outfile.close()
		sys.stderr.write("wrote settings to 'mdmcleaner.config'\n")
	
	if args.command == "show_configs":
		sys.stderr.write("setting\tvalue\tsource\n")
		# ~ print(configs)
		# ~ print("*"*100)
		# ~ print(settings_source)
		for c in configs:
			sys.stderr.write("{}\t{}\t{}\n".format(c,configs[c],settings_source[c]))

	if args.command == "acc2taxpath":
		import getdb
		db = getdb.taxdb(configs)
		if args.accessions[0] == "interactive":
			# ~ import pdb; pdb.set_trace()
			sys.stderr.write("\n\nrunning {} acc2taxpath {} in interactive mode. type \"exit\" or press \"ctrl+D\" to exit\n\n".format(os.path.basename(sys.argv[0]), __version__))
			userinput = ""
			while True:
				try:
					userinput = input("enter accession here:")
					if userinput == "":
						continue
					elif userinput == "exit":
						break
					sys.stderr.write("\n--{}--\n".format(userinput))
					taxid, _ = db.acc2taxid(userinput)
					taxpath = db.taxid2pathstring(taxid)
					sys.stderr.write("{}\n".format(taxpath))
				except EOFError:
					break
		else:
			for acc in args.accessions:
				taxid, _ = db.acc2taxid(acc)
				# ~ sys.stderr.write("=={}==".format(taxid))
				taxpath = db.taxid2pathstring(taxid)
				sys.stderr.write("{}\t{}".format(taxid, taxpath))
		sys.stderr.write("\n")
	
	if args.command == "refdb_contams":
		import review_refdbcontams
		blacklist_additions = review_refdbcontams.read_ambiguity_report(args.ambiguity_report, configs, args.tempbasename)
		print("\n------> writing {} blacklist additions to {}".format(len(blacklist_additions), args.outblacklist))
		write_blacklist(blacklist_additions, args.outblacklist)
main()

			
