#!/usr/bin/env python
""" 
extracting protein-coding as rRNA genes from assemblies, and identifying markergenes.
"""

# note to self:
# cutoff values were determined in different ways:
# for tigrfam and pfam: strict and sensitive values were parsed from the "GA" and "NC" fields, respectively. moderate values were calculated as the average betrween the respective strict and sensitive cutoffs
# for cogs all marker models were aligned seperately against the component merkergene-alignments and against all nonmarker-cog-alignments
# strict: the higher value of the cutoff that yielded 95% of the true positives (component markergene-alignents) and the cutoff that yielded less than 5% false positives (all nonmarker alignments)
# sensitive: the LOWER value of the cutoff that yielded 95% of the true positives (component markergene-alignents) and the cutoff that yielded less than 5% false positives (all nonmarker alignments)
# moderate: the average between strict and sensitive

import subprocess
import sys
import os
from misc import openfile
#import Bio.SearchIO.HmmerIO.hmmer3_domtab.Hmmer3DomtabHmmhitParser #probably better to parse it my self

#currently the marker-hmms only encompass universal SINGLE-COPY genes. It would be interesting to include the multicopy-universal genes as well! --> parse the COG-database for this...?

libpath = os.path.dirname(os.path.realpath(__file__))
hmmpath = os.path.realpath(os.path.join(libpath, "../hmms/"))
hmmpath_prok = os.path.realpath(os.path.join(hmmpath, "prok/"))
hmmpath_bact = os.path.realpath(os.path.join(hmmpath, "bact/"))
hmmpath_arch = os.path.realpath(os.path.join(hmmpath, "arch/"))
hmmpathdict={	"prok" : [hmmpath_prok], \
				"bact" : [hmmpath_prok, hmmpath_bact], \
				"arch" : [hmmpath_prok, hmmpath_arch], \
				"all" : [hmmpath_prok, hmmpath_bact, hmmpath_arch] }
cutofftablefile = os.path.join(hmmpath, "cutofftable_combined.tsv")
#each path in hmmpathdict should contain a number of hmm files named e.g. COG.hmm, PFAM.hmm or TIGR.hmm, containing concatenated hmm models for each level/db-type



def split_fasta_for_parallelruns(infasta, outputbasename, threads = 2):
	"""
	splits large multifastas into several portions, to enable parallel runs of single threaded processes, such as rnammer or prodigal
	requires subdivide_multifas.py
	"""
	from subdivide_multifas import subdivide
	return subdivide(fastaname, threads, outputbasename)

# def run_parallel(input_list, 

def runprodigal(infasta, outfilename, prodigal="prodigal"): #todo: allow piping via stdin
	"""
	creates a protein fasta file of all CDS identified in the inputfasta via prodigal,
	all other prodigal output is ignored
	prodigal is called using the "-p meta" argument for metagenomes, in the assumption that the input fasta MAY consist of mutliple organisms
	The return value is simply the value of 'outfilename'
	"""
	prodigal_cmd = subprocess.run([prodigal, "-a", outfilename, "-p", "meta", "-q", "-o", "/dev/null", "-i", infasta], stdout = subprocess.PIPE, stderr = subprocess.PIPE, text = True) # TODO: maybe add option to change translation table ("-g")? Although table 11 should be general enough? 
	try:
		prodigal_cmd.check_returncode()
	except Exception:
		sys.stderr.write(hmmsearch_cmd.stderr)
		raise Exception("\nERROR: Something went wrong while trying to call prodigal...\n")
	return outfilename
	
def runbarrnap(infasta, outfilename, barrnap="barrnap"): #todo: allow piping via stdin
	barrnap_cmd = subprocess.run([barrnap, "-a", outfilename, "-p", "meta", "-q", "-o", "/dev/null"], stdout = subprocess.PIPE, stderr = subprocess.PIPE, text = True) 
	try:
		barrnap_cmd.check_returncode()
	except Exception:
		sys.stderr.write(hmmsearch_cmd.stderr)
		raise Exception("\nERROR: Something went wrong while trying to call barrnap...\n")
	return outfilename

def runrnammer(infasta, outfilename, threads = 1): #todo: allow piping via stdin
	pass

def hmmersearch(hmmsearch, model, query, outfilename, score_cutoff = None, eval_cutoff = None, threads = 1):# todo: strict parameters = gathering threshold (GA), sensitive parameters = noise cutoff (NC)
	"""
	runs hmmsearch
	score and/or evalue cutoffs can be specified seperately.
	if neither score-, nor eval_cutoff are supplied, it will try to obtain the cutoff values from the "GA" field of the model ("Gathering Threshold"; if available).
	alternatively the score_cutoff can be non-explicetly set either as "strict" or "sensitive". In this case the evalue_cutoff will ignored and the following cutoffs will be used from the hmm file:
		- strict: GA (=Gathering threshold)
		- sensitive: NC (Noise Cutoff)
	note that for self-built hmms without "GA" and "NC" keys, cutoffs will need to be specified explicitely.
	"""
	eval_cutoff_arg, score_cutoff_arg = [], []
	if (eval_cutoff == None and score_cutoff == None) or score_cutoff == "strict":
		score_cutoff = ["--cut_tc"] # use trusted cutoff of hmm model (if available). consider only using gathering threshold (GA) uinstead
	elif score_cutoff == "sensitive":
		score_cutoff = ["--cut_nc"] # use noise cutoff of hmm model (if available).	
	elif score_cutoff == "moderate":
		score_cutoff = ["--cut_ga"] # use gathering cutoff of hmm model (if available)
	else:	
		if eval_cutoff != None:
			eval_cutoff_arg = ["-E", eval_cutoff]
		if score_cutoff != None:
			score_cutoff_arg = ["-T", score_cutoff]
	if type(query) == str: #TODO: This assumes "if query is a string, it must be a filename." That is obviously BS! implement a check that tests if string is a fasta-record!
		hmmsearch_cmd = subprocess.run([hmmsearch, "--noali", "--cpu", str(threads), "--domtblout", outfilename] \
		 + eval_cutoff_arg + score_cutoff_arg + [model, query], stdout = subprocess.PIPE, stderr = subprocess.PIPE, text = True)
	else:	#otherwise it must be something to pipe via stdin
		pass #TODO: add this!
	try:
		hmmsearch_cmd.check_returncode()
	except Exception:# Todo: define/choose more detailed exception categories
		sys.stderr.write(hmmsearch_cmd.stderr)
		raise Exception("\nERROR: something went wrong while trying to call hmmsearch...\n")
	return outfilename

### NOTE TO SELF: perform hmmsearch always mit "sensitive" cutoff, and only PARSE hits with higher cutoffs --> enables reanalyses with different cutoffs without redoing hmmsearch!

def get_cutoff_dict(cutofffilename = cutofftablefile): #todo lookupfile with cutoffs for ALL used models. TODO: better: parse this from model.hmm files (require GA, TC and NC fields in all used models!)
	"""
	reads cutoff values from cutoff_file into a dictonary
	each model is represented as a seperate line with 4 columns:
		- first column = model name
		- second column = strict cutoff
		- third column = moderate cutoff
		- fourth column = sensitive cutoff"
	"""
	cutofffile = openfile(cutofffilename)
	cutoff_dict = {}
	for line in cutofffile:
		if line.startswith("#"):
			continue
		tokens = line.split()
		model = tokens[0].split(".")[0]
		strict = float(tokens[1])
		moderate = float(tokens[2])
		sensitive = float(tokens[3])
		cutoff_dict[model] = {"strict" : strict, "moderate" : moderate, "sensitive" : sensitive}
	return cutoff_dict
	
def parse_hmmer(hmmerfile, cutoff_dict = cutofftablefile, cmode = "moderate", prev_results = None):
	"""
	parses hmmer result file, using cutoff-thresholds passed as a dictionary "cutoff_dict", as returned by "get_cutoff_dict()"
	returns a dictionary containing protein-identifiers as keys and respective marker-designations as values for each hmm hit that passed cutoff criteria
	cutoff_dict should be a dictinary with the "strict", "moderate" and "sensitive" cutoff-values for each marker-model, but CAN also be a filename from which to parse that dict (default = parse from default file)
	prev_results may be a previous hit_dictionary that shlud be updated with hits form the present one
	"""
	assert cmode in  ["strict", "moderate", "sensitive"], "\nError: dont recognize mode \"{}\"! mode must be one of [\"strict\", \"moderate\", \"sensitive\"]\n"
	if type(cutoff_dict) != dict: #alternative for parsing cutoff_dict will be read from a file (better to pass it as dict, though)
		#TODO: add logger message that cutoff dict is being read from file
		cutoff_dict = get_cutoff_dict(cutoff_dict)
	infile = openfile(hmmerfile)
	if prev_results == None:
		markerdict = {}
	elif type(prev_results) == dict:
		markerdict = prev_results
	else:
		raise RuntimeError("\nArgument 'prev_results' should be either None or of type dict\n")
	for line in infile:
		if line.startswith("#"):
			continue
		tokens = line.split()
		prot = tokens[0]
		marker = tokens[4].split(".")[0]
		fscore = float(tokens[7])
		#dscore = float(tokens[13]) #not sure if i will use this
		#print(" found '{}' (which is marker '{}') with score = {}.  cutoff is {}".format(prot, marker, fscore, cutoff_dict[marker][cmode]))
		if fscore < cutoff_dict[marker][cmode]:
			#print("    --> score not goud enough")
			continue 
		if prot not in markerdict:
			#print("            --> {} is being stored".format(marker))
			markerdict[prot] = { "marker" : marker, "fscore" : fscore } #may need to add dscore here
			#print(markerdict)
	return markerdict 

def get_markernames(proteinfastafile, cutoff_dict = cutofftablefile, hmmsearch = "hmmsearch", outdir = ".", cmode = "moderate", level = "prok", threads = "1"):
	"""
	runs hmmersearch and and parse_hmmer on designated proteinfasta using models for designated level. Requires a cutoff_dict as returned by "get_gutoff_dict()"
	cutoff_dict should be a dictinary with the "strict", "moderate" and "sensitive" cutoff-values for each marker-model, but CAN also be a filename from which to parse that dict (default = parse from default file)
	returns a nested dictionary containing protein-identifiers as main keys and subdictionaries with respective marker-designations (key = "marker") and score values (key = "fscore")  as values for each hmm hit that passed cutoff criteria
	"""
	assert level in ["prok", "bact", "arch", "all"], "\nError: dont recognize level \"{}\"! mode must be one of [\"prok\", \"bact\", \"arch\", \"all\"]\n"
	assert cmode in  ["strict", "moderate", "sensitive"], "\nError: dont recognize mode \"{}\"! mode must be one of [\"strict\", \"moderate\", \"sensitive\"]\n" 
	if type(cutoff_dict) != dict: #alternative for parsing cutoff_dict will be read from a file (better to pass it as dict, though)
		#TODO: add logger message that cutoff dict is being read from file
		cutoff_dict = get_cutoff_dict(cutoff_dict)
	list_of_markerdicts = []
	print("getting markerdicts")
	for hmmpath in hmmpathdict[level]:
		hmmfiles = [ os.path.join(hmmpath, hmmfile) for hmmfile in os.listdir(hmmpath) if hmmfile.endswith(".hmm") ]
		markerdict = {}
		for hmmfile in hmmfiles:
			sys.stderr.write("\nsearching {} ...".format(hmmfile))
			outfile = hmmersearch(hmmsearch, hmmfile, proteinfastafile, os.path.join(outdir, os.path.basename(hmmfile) + ".domtblout"), "sensitive", None, threads)
			markerdict = parse_hmmer(outfile, cutoff_dict, cmode, markerdict)
			print(hmmfile)
			print(len(markerdict))
			#print(markerdict)
			print("--------------------------")
		list_of_markerdicts.append(markerdict)
	return deduplicate_markers(list_of_markerdicts) #list_of_markerdicts will be in tis order: [prok[, bact[, arch]]]

def deduplicate_markers(list_of_markerdicts): # For proteins with hits to different models, just keep the hit with the highest score. This function is a highly convoluted way to do this, but it is late and my brain is tired
	print("deduplicating")
	print("{}".format(", ".join([str(len(x)) for x in list_of_markerdicts])))
	keys = set([ key for md in list_of_markerdicts for key in md ])
	for key in keys:
		a, b = 0, 1
		while a < len(list_of_markerdicts) and b < len(list_of_markerdicts):
			if key in list_of_markerdicts[a]:		
				while b < len(list_of_markerdicts[a:]):
					if key in list_of_markerdicts[b]:
						if list_of_markerdicts[a][key]["fscore"] >= list_of_markerdicts[b][key]["fscore"]:
							list_of_markerdicts[b].pop(key)
						else:
							list_of_markerdicts[a].pop(key)
					b += 1
			a += 1
			b += 1
	print("whats left:")
	print("{}".format(", ".join([str(len(x)) for x in list_of_markerdicts])))
	return list_of_markerdicts
			
	

def __get_markerseqs(proteinfastafile, markerdict):
	"""
	returns a list of proteinsequences corresponding to markers found in markerdict
	marker designation and score alue are written to the description of each protein sequence
	"""
	from Bio import SeqIO
	protfastafile = openfile(proteinfastafile)
	protrecords  = SeqIO.parse(protfastafile, "fasta")
	markerlist = []
	for prot in protrecords:
		#print("checking if '{}' in markerdict\n".format(prot.id))
		if prot.id in markerdict:
			prot.description = "marker={};score={};desc={}".format(markerdict[prot.id]["marker"], markerdict[prot.id]["fscore"], prot.description)
			markerlist.append(prot)
	#print(markerdict)
	return markerlist

def get_markers(proteinfastafile, cutoff_dict = cutofftablefile, cmode = "moderate", level = "prok", outfile_basename = "markerprots", threads = 1):
	"""
	writes fasta sequences of detected markergenes in fasta format to outfile, with marker-designation and hmm score value in description
	'cmode' refers to "cutoff_mode" and can be one of ["strict", "moderate", or "sensitive"]. Sets the score cutoff_values to use for selecting hits. For each marker-designation and cutoff-mode 
	cutoff_dict should be a dictinary with the "strict", "moderate" and "sensitive" cutoff-values for each marker-model, but CAN also be a filename from which to parse that dict (default = parse from default file)
	return value is simply the name of the outfile
	"""
	from Bio import SeqIO
	levelorder = ["prok", "bact", "arch"]
	outdir = os.path.dirname(outfile_basename)
	if type(cutoff_dict) != dict: #alternative for parsing cutoff_dict will be read from a file (better to pass it as dict, though)
		#TODO: add logger message that cutoff dict is being read from file
		cutoff_dict = get_cutoff_dict(cutoff_dict)
	list_of_markerdicts = get_markernames(proteinfastafile, cutoff_dict, hmmsearch = "hmmsearch", outdir = outdir, cmode = "moderate", level = level, threads = threads)
	outfilelist = []
	for i in range(len(list_of_markerdicts)):
		markerseqs = __get_markerseqs(proteinfastafile, list_of_markerdicts[i])
		outfilename = "{}_{}.faa".format(outfile_basename, levelorder[i])
		outfile = openfile(outfilename, "wt")
		SeqIO.write(markerseqs, outfile, "fasta")
		outfile.close()
		outfilelist.append(outfilename)
	return outfilelist
	
def write_markerdict(markerdict, outfilename):
	"""
	writes the marker dictionary, obtained by get_markernames(), to an overview file in tab-seperated text-table format
	return value is simply the name of the outfile
	"""
	outfile = openfile(outfilename, "wt")
	for m in markerdict:
		outline = "{}\t{}\n".format(m, "\t".join([ str(markerdict[m][v]) for v in markerdict[m].keys() ]))
		outfile.write(outline)
	return outfilename
	
def _test_markernames():
	sys.stderr.write("\ntesting get_markernames...")
	proteinfastafile = sys.argv[1]
	cutofftable = sys.argv[2]
	sys.stderr.write("\nreading cutofftable")
	cutoffdict = get_cutoff_dict(cutofftable)
	sys.stderr.write("\nsearching markers")
	markerdict = get_markernames(proteinfastafile, cutoffdict, hmmsearch = "hmmsearch", outdir = ".", cmode = "moderate", level = "prok", threads = "4")
	sys.stderr.write("\nwriting results\n")
	write_markerdict(markerdict, "delmetestresults.tsv")
	
def _test_basicmarkers():
	infasta = sys.argv[1]
	tempdir = sys.argv[2]
	if not os.path.exists(tempdir):
		os.mkdir(tempdir) #todo: implement tempfile module if available a base module
	#else:
		#raise Exception("\n'{}' already exists\n".format(tempdir))
	cutofftable = os.path.join(hmmpath, "cutofftable_combined.tsv")
	cutoff_dict = get_cutoff_dict(cutofftable)
	sys.stderr.write("\nrunning prodigal...\n")
	protfasta = runprodigal(infasta, os.path.join(tempdir, "delme_protfasta"), prodigal="prodigal")
	#protfasta = os.path.join(tempdir, "delme_protfasta")
	#todo: create a "runparallel function in misc or here
	level = "all"
	sys.stderr.write("\nextracting markers for level {}\n".format(level))
	outfastalist = get_markers(protfasta, cutoff_dict, level = level, outfile_basename = os.path.join(tempdir, "markers".format(level)), threads = 4)
	sys.stderr.write("  --> created files: '{}'".format(", ".join(outfastalist)))
	#todo: implement automatic blasts
	#todo implement actual lca
	
def main():
	#_test_markernames()
	_test_basicmarkers()

if __name__ == '__main__':
	main()
