#!/usr/bin/env python
import sys
if sys.version_info[0] < 3:
	sys.exit("\nERROR: python 3 or higher required!\n")
from mdmcleaner import check_dependencies
import os
import argparse
from mdmcleaner import misc
from mdmcleaner.misc import openfile
from mdmcleaner._version import __version__

setting_keys = ["blacklistfile","blastp", "blastn", "blastdbcmd", "diamond", "barrnap", "rnammer", "hmmsearch", "aragorn", "prodigal", "db_basedir","db_type", "blastdb_diamond", "blastdb_blastp", "blastdb_blastn", "threads"] # todo: complete this list 
MDMCLEANER_LIBPATH = os.path.dirname(os.path.realpath(misc.__file__))

def write_blacklist(blacklist, outfilename):
	with misc.openfile(outfilename, "at") as outfile:
		for b in blacklist:
			outfile.write("{}\n".format(b))

def find_global_configfile():
	if os.path.exists(os.path.join(MDMCLEANER_LIBPATH, "mdmcleaner.config")) and os.path.isfile(os.path.join(MDMCLEANER_LIBPATH, "mdmcleaner.config")):
		return os.path.join(MDMCLEANER_LIBPATH, "mdmcleaner.config")
	sys.exit("\nError: a \"mdmcleaner.config\" file should exist under {}, but doesnt!\n".format(MDMCLEANER_LIBPATH))

def find_local_configfile():
	cwd = os.getcwd()
	if os.path.exists(os.path.join(cwd, "mdmcleaner.config")) and os.path.isfile(os.path.join(cwd, "mdmcleaner.config")):
		return os.path.join(cwd, "mdmcleaner.config")


class config_object(object):
	execs = {	"blastpath" : ["blastn", "blastp", "makeblastdb", "blastdbcmd"],\
				"diamondpath" : ["diamond"],\
				"barrnappath" : ["barrnap"],\
				"hmmerpath" : ["hmmsearch"],\
				"aragornpath" : ["aragorn"],\
				"prodigalpath" : ["prodigal"] } #lists which excecutables should be found in which path. todo: integrate this with check_dependencies...
	config_file_exepathkeys = list(execs.keys())
	config_file_misckeys = ["threads"]
	config_file_dbkeys = ["db_basedir", "db_type", "blacklistfile"]
	config_file_setting_keys = config_file_exepathkeys + config_file_misckeys + config_file_dbkeys
					
	def __init__(self, args, read_blacklist = True):
		self.configfile_hierarchy = [ cf for cf in [find_global_configfile(), args.configfile] if cf != None ]
		self.settings = {key : [] for key in self.config_file_setting_keys}
		self.settings_source = {key : "default" for key in self.config_file_setting_keys}
		self.blacklist = set()
		self.read_configs(args)
		for execpath in self.execs:
			for e in self.execs[execpath]:
				if len(self.settings[execpath]) == 0:
					self.settings[e] = e
				else:
					self.settings[e] = os.path.join(self.settings[execpath][0], e)
		if read_blacklist:
			self.blacklist = self.read_blacklistfiles()	
	
	def print_settings(self):
		sys.stderr.write("\n\tsettings:\n")
		for key in self.settings:
			if self.settings[key] != []:
				sys.stderr.write("\t\t{} = '{}'\n".format(key, self.settings[key]))
		sys.stderr.write("\n")	
	
	def read_configs(self, args):
		"""
		Reads the config files in hierarchical order (first global, then local), with the later configs always overriding the previous in case of conflicts
		Config files must be tab-seperated text files (may be compressed though), with setting names in the first column, and the corresponding setting value(s) in subsequent columns.
		Setting values are always read as lists (except "threads" which is read as integer)
		Unknown setting names will just be ignored. However, comments should optimally be marked with "#"
		"""
		import os
		self.settings["blacklistfile"] = []
		if "ignore_default_blacklist" in vars(args) and args.ignore_default_blacklist == True:
			self.settings_source["blacklistfile"] = None
		else:
			self.mdmcleaner_lib_path = MDMCLEANER_LIBPATH
			default_blacklist = os.path.join(self.mdmcleaner_lib_path, "blacklist.list")
			self.settings["blacklistfile"].append(default_blacklist)
			self.settings_source["blacklistfile"] = "default"
		for config in self.configfile_hierarchy: #todo: maybe instead of replacing blacklist with that in config-file, just append it to the list of blacklists?
			sys.stderr.write("\nreading settings from configfile: \"{}\"\n".format(config))
			with openfile(config) as configfile:
				for line in configfile:
					tokens = line.strip().split("#",1)[0].strip().split()
					if len(tokens) <= 1:
						continue
					if tokens[0] in self.config_file_setting_keys:
						self.settings[tokens[0]] = tokens[1:]
						self.settings_source[tokens[0]] = config
					else:
						sys.stderr.write("\tWARNING: unknown setting key \"{}\" in config file '{}'--> ignoring it!\n".format(tokens[0], config))
		if "threads" in vars(args) and args.threads:
			self.settings["threads"] = args.threads
		else:
			self.settings["threads"] = int(self.settings["threads"][0])
		if "blacklistfile" in vars(args) and args.blacklistfile != None:
			self.settings["blacklistfile"].append(args.blacklistfile)	

	def read_blacklistfiles(self):
		blacklist = []
		for blacklistfile in self.settings["blacklistfile"]:
			if not (os.path.exists(blacklistfile) and os.path.isfile(blacklistfile)):
				sys.stderr.write("\nWARNING: could not file blacklistfile : \"{}\"! --> skipping it!\n".format(blacklistfile))
				continue
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
	cleansagmag_args.add_argument("--outblacklist", action="store", dest="outblacklist", default = "new_blacklist_additions.tsv", help = "Outputfile for new blacklist additions. If a preexisting file is selected, additions will be appended to end of that file")
	cleansagmag_args.add_argument("-v", "--version", action = "version", version='%(prog)s {version}'.format(version=__version__))
	cleansagmag_args.add_argument("-c", "--config", action = "store", dest = "configfile", default = find_local_configfile(), help = "provide a local config file with basic settings (such as the location of database-files). default: looks for config files named 'mdmcleaner.config' in current working directory. settings in the local config file will override settings in the global config file '{}'".format(os.path.join(os.path.dirname(os.path.abspath(__file__)), "mdmcleaner.config")))
	cleansagmag_args.add_argument("-t", "--threads", action = "store", dest = "threads", type = int, help = "Number of threads to use. (default = use setting from config file)")
	cleansagmag_args.add_argument("-f", "--force", action = "store_true", dest = "force", default = False, help = "Force reclassification of pre-existing blast-results")
	# ~ cleansagmag_args.add_argument("--blast2pass", action = "store_true", dest = "blast2pass", default = "False", help = "add a second-pass blastx blast-run for all contigs without any classification on blastp level (default: False)")
	cleansagmag_args.add_argument("--overview_files_basename", action = "store", dest = "overview_basename", default = "overview", help = "basename for overviewfiles (default=\"overview\"")
	cleansagmag_args.add_argument("-b", "--blacklistfile", action = "store", dest = "blacklistfile", default = None, help = "File listing reference-DB sequence-names that should be ignored during blast-analyses (e.g. known refDB-contaminations...")
	cleansagmag_args.add_argument("--no_filterfasta", action = "store_true", dest = "no_filterfasta", default = False, help = "Do not write filtered contigs to final output fastas (Default = False)")
	cleansagmag_args.add_argument("--ignore_default_blacklist", action = "store_true", dest = "ignore_default_blacklist", default = False, help = "Ignore the default blacklist (Default = False)")
	cleansagmag_args.add_argument("--fast_run", action = "store_true", dest = "fast_run", default = False, help = "skip detailed analyses of potential reference ambiguities (runs may be faster but also classification may be less exact, and potential reference database contaminations will not be verified)")


	makedb_args = subparsers.add_parser("makedb", help = "Download and create MDMcleaner database")
	makedb_args.add_argument("-o", "--outdir", action = "store", dest = "outdir", default = None, help = "target base directory for reference-data. may not be the current working directory. Needs >100GB space! Default = './db' (if downloading publication-dataset or if not already specified in config file)")
	makedb_args.add_argument("-c", "--config", action = "store", dest = "configfile", default = find_local_configfile(), help = "provide a local config file with the target location to store database-files. default: looks for config files named 'mdmcleaner.config' in current working directory. settings in the local config file will override settings in the global config file '{}'".format(os.path.join(os.path.dirname(os.path.abspath(__file__)), "mdmcleaner.config")))
	makedb_args.add_argument("--get_pub_data", action = "store_true", dest = "get_pub_data", default = False, help = "simply download the reference dataset from the MDMcleaner publication from zenodo (WARNING: this is outdated and therefore NOT optimal! Only meant for testing/verification purposes)")
	makedb_args.add_argument("--verbose", action = "store_true", dest = "verbose", default = False, help = "verbose output (download progress etc)") #todo: finish implementing
	makedb_args.add_argument("--quiet", action = "store_true", dest = "quiet", default = False, help = "quiet mode (suppress any status messages except Errors and Warnings)") #todo: implement
	
	get_marker_args = subparsers.add_parser("get_markers", help = "extracts protein coding and/or rRNA gene sequences from input genome(s)")
	get_marker_args.add_argument("-i", "--input_fastas", action = "store", dest = "input_fastas", nargs = "+", required = True, help = "input fasta(s). May be gzip-compressed")
	get_marker_args.add_argument("-m", "--markertype", action = "store", dest = "markertype", default = "all", choices = ["rrna", "trna", "totalprots", "markerprots", "all"], help = "type of marker gene that should be extracted (default = 'all')")
	get_marker_args.add_argument("-o", "--outdir", dest = "outdir", default = ".", help = "Output directory (will be created if it does not exist). Default = '.'")
	get_marker_args.add_argument("-c", "--config", action = "store", dest = "configfile", default = find_local_configfile(), help = "provide a local config file with the target location to store database-files. default: looks for config files named 'mdmcleaner.config' in current working directory. settings in the local config file will override settings in the global config file '{}'".format(os.path.join(os.path.dirname(os.path.abspath(__file__)), "mdmcleaner.config")))
	get_marker_args.add_argument("-t", "--threads", action="store", dest="threads", type = int, help = "number of threads to use (default = use setting from config file)")
	get_marker_args.add_argument("-M", "--mincontiglength", action="store", dest="mincontiglength", type = int, default = 0, help = "minimum contig length (contigs shorter than this will be ignored)")	

	completeness_args = subparsers.add_parser("completeness", help = "estimate completeness (roughly based on presence of universally required tRNA types). Results are printed directly to stdout")
	completeness_args.add_argument("-i", "--input_fastas", action = "store", dest = "input_fastas", nargs = "+", required = True, help = "input fasta(s). May be gzip-compressed")
	completeness_args.add_argument("-c", "--config", action = "store", dest = "configfile", default = find_local_configfile(), help = "provide a local config file with the target location to store database-files. default: looks for config files named 'mdmcleaner.config' in current working directory. settings in the local config file will override settings in the global config file '{}'".format(os.path.join(os.path.dirname(os.path.abspath(__file__)), "mdmcleaner.config")))
	completeness_args.add_argument("-t", "--threads", action="store", dest="threads", type = int, help = "number of threads to use (default = use setting from config file)")
	completeness_args.add_argument("-M", "--mincontiglength", action="store", dest="mincontiglength", type = int, default = 0, help = "minimum contig length (contigs shorter than this will be ignored)")	

	acc2taxpath = subparsers.add_parser("acc2taxpath", help = "Get full taxonomic path assorciated with a specific acession number")
	acc2taxpath.add_argument("accessions", action = "store", nargs = "+", help = "(space seperated list of) input accessions. Or just pass \"interactive\" for interactive mode")
	acc2taxpath.add_argument("-c", "--config", action = "store", dest = "configfile", default = find_local_configfile(), help = "provide a local config file with the target location to store database-files. default: looks for config files named 'mdmcleaner.config' in current working directory. settings in the local config file will override settings in the global config file '{}'".format(os.path.join(os.path.dirname(os.path.abspath(__file__)), "mdmcleaner.config")))

	evaluate_refdbcontam_args = subparsers.add_parser("refdb_contams", help = "EXPERIMENTAL: evaluate potentiel refDB-contaminations")
	evaluate_refdbcontam_args.add_argument("ambiguity_report", help =  "The ambiguity report file to analyse (Usually './overview_refdb_ambiguities.tsv')")
	evaluate_refdbcontam_args.add_argument("-c", "--config", action = "store", dest = "configfile", default = find_local_configfile(), help = "provide a local config file with the target location to store database-files. default: looks for config files named 'mdmcleaner.config' in current working directory. settings in the local config file will override settings in the global config file '{}'".format(os.path.join(os.path.dirname(os.path.abspath(__file__)), "mdmcleaner.config")))
	evaluate_refdbcontam_args.add_argument("-t", "--threads", action="store", dest="threads", type = int, help = "number of threads to use (default = use setting from config file)")
	evaluate_refdbcontam_args.add_argument("-o", "--outblacklist", action="store", dest="outblacklist", default = "new_blacklist_additions.tsv", help = "Outputfile for new blacklist additions. If a preexisting file is selected, additions will be appended to end of that file")
	# ~ evaluate_refdbcontam_args.add_argument("--tempbasename", action="store", dest="tempbasename", default = "mdmtempfile", help = "basename for temporary files (default = 'mdmtempfile'")
	
	
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

	check_dependencies_args = subparsers.add_parser("check_dependencies", help = "checks if all dependencies for MDMcleaner are being met")
	check_dependencies_args.add_argument("-c", "--config", action = "store", dest = "configfile", default = find_local_configfile(), help = "provide a local config file with basic settings (such as the location of database-files). default: looks for config files named 'mdmcleaner.config' in current working directory. settings in the local config file will override settings in the global config file '{}'".format(os.path.join(os.path.dirname(os.path.abspath(__file__)), "mdmcleaner.config")))

	version_args = subparsers.add_parser("version", help = "show version info and exit")
	if len(sys.argv)==1:
		myparser.print_help(sys.stderr)
		sys.exit(1)
	args = myparser.parse_args()
	
	if args.command == "version":
		sys.stderr.write("MDMcleaner v{}\n".format(__version__))

	# ~ if args.command in ["clean", "makedb", "show_configs", "get_markers", "acc2taxpath", "refdb_contams", "check_dependencies"]:
		# ~ configfile_hierarchy = [ cf for cf in [find_global_configfile(), args.configfile] if cf != None ]
		# ~ configs, settings_source = read_configs(configfile_hierarchy, args)
		# ~ sys.stderr.write("\n\nSETTINGS:\n" + pprint.pformat(configs)+ "\n\n")
		# ~ if args.command in ["clean", "get_markers", "refdb_contams"]:
			# ~ configs["blacklist"] =  _read_blacklistfiles(configs["blacklistfile"])

	sys.stderr.write("You are running the following MDMcleaner command:\n\t '{}'".format(" ".join(sys.argv)))

	if args.command in ["clean", "makedb", "show_configs", "get_markers", "completeness", "acc2taxpath", "refdb_contams", "check_dependencies"]:
		if args.command in ["clean", "get_markers", "refdb_contams"]:
			configs = config_object(args, read_blacklist = True)
		else:
			configs = config_object(args)
		configs.print_settings()


	# ~ import pdb; pdb.set_trace()
	if args.command == "clean":
		from mdmcleaner import clean
		check_dependencies.check_dependencies("blastp", "blastn", "diamond", "aragorn", "barrnap", "hmmsearch", configs=configs)
		blacklist_additions = clean.main(args, configs)
		if not args.fast_run and blacklist_additions != None:
			print("\n------> writing {} blacklist additions to {}".format(len(blacklist_additions), args.outblacklist))
			write_blacklist(blacklist_additions, args.outblacklist)
	
	if args.command == "makedb":
		from mdmcleaner  import read_gtdb_taxonomy
		if args.outdir == None:
			if args.get_pub_data:
				args.outdir = "./db"
			else:
				if not ("db_basedir" in configs.settings and len(configs.settings["db_basedir"]) != 0): 
					makedb_args.error("either 'outdir' must be specified as argument or 'db_basedir' needs to be specified in config file")
				args.outdir = os.path.join(configs.settings["db_basedir"][0], configs.settings["db_type"][0]) # todo: read_gtdb_taxonomy should only get basedir as target-dir and simply assume the gtdb part!
		elif not args.get_pub_data:
			args.outdir = os.path.join(args.outdir, configs.settings["db_type"][0]) # todo: read_gtdb_taxonomy should only get basedir as target-dir and simply assume the gtdb part!
		check_dependencies.check_dependencies("makeblastdb", "diamond", "wget", configs=configs)
		read_gtdb_taxonomy.main(args, configs)
		
	if args.command == "get_markers":
		from mdmcleaner import getmarkers
		check_dependencies.check_dependencies("aragorn", "barrnap", "hmmsearch", configs=configs)
		getmarkers.get_only_marker_seqs(args, configs)
	
	if args.command == "completeness":
		check_dependencies.check_dependencies("aragorn", "barrnap", configs=configs)
		from mdmcleaner import getmarkers
		args.outdir = None
		getmarkers.get_only_trna_completeness(args, configs)

	if args.command == "set_configs":
		args.configfile = None
		conf = config_object(args)
		# ~ import pdb; pdb.set_trace()
		if args.scope == "global":
			outfile = misc.openfile(conf.configfile_hierarchy[0], "wt")
		elif args.scope == "local":
			outfile = misc.openfile("mdmcleaner.config", "wt")
		arg_settings = vars(args)
		for a in arg_settings:
			if a != "scope" and arg_settings[a] != None:
				conf.settings[a] = [arg_settings[a]]
		for c in conf.settings:
			if c in conf.config_file_setting_keys and conf.settings[c]:
				if isinstance(conf.settings[c], list): #todo: this if statement is only necessary because of converting thread-entries to ints instead of lists. not sure why i did that. --> change it back to consistently sava all settings as lists, if possible
					outfile.write("{}\t{}\n".format(c, conf.settings[c][0]))
				else:
					outfile.write("{}\t{}\n".format(c, conf.settings[c]))
		outfile.close()
		sys.stderr.write("wrote settings to 'mdmcleaner.config'\n")
	
	if args.command == "show_configs":
		sys.stderr.write("setting\tvalue\tsource\n")
		# ~ print(configs)
		# ~ print("*"*100)
		# ~ print(settings_source)
		for c in configs.config_file_setting_keys:
			sys.stderr.write("{}\t{}\t{}\n".format(c,configs.settings[c],configs.settings_source[c]))

	if args.command == "acc2taxpath":
		from mdmcleaner import getdb
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

	if args.command == "check_dependencies":
		check_dependencies.check_dependencies(configs=configs)
		sys.stderr.write("\nSUCCESS: all dependencies are being met\n")

	if args.command == "refdb_contams":
		from mdmcleaner import review_refdbcontams
		check_dependencies.check_dependencies("blastp", "blastn", "diamond", "wget", configs=configs)
		blacklist_additions = review_refdbcontams.read_ambiguity_report(args.ambiguity_report, configs)
		print("\n------> writing {} blacklist additions to {}".format(len(blacklist_additions), args.outblacklist))
		write_blacklist(blacklist_additions, args.outblacklist)

if __name__ == '__main__':
	main()

			
