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
from Bio import SeqIO
import misc
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

protmarkerlevel_dict = { 0 : "prok_marker", 1 : "bac_marker", 2 : "arc_marker" }

def split_fasta_for_parallelruns(infasta, minlength = 0, number_of_fractions = 2, outfilebasename = None):
	"""
	splits large multifastas into several portions, to enable parallel runs of single threaded processes, such as rnammer or prodigal
	requires subdivide_multifas.py
	returns a list of lists of sequence records
	each list of sequence records should then be passed to the stdin of a seperate RNAmmer or prodigal call (not necessary for barrnap, because that already supports multithreading)
	"""
	import random
	from Bio import SeqIO

	sys.stderr.write("\nsubdividing contigs of {} for multiprocessing\n".format(infasta))
	fastafile = openfile(infasta)
	records = SeqIO.parse(fastafile, "fasta")
	contigdict = {}
	#print("whooooooo --{}-- whooooooo".format(number_of_fractions))
	outlist = [[] for x in range(int(number_of_fractions))]
	contigcounter = 0
	index = 0
	direction = 1

	for record in records: #distribute contigs evenly over fractions by iterating up and down the fractions again and again --> ensures even size distribution if possible...
		if len(record) < minlength:
			continue
		contigcounter += 1
		contigdict[record.id] = {"contiglen": [len(record)], "totalprotcount" : [0], "ssu_rRNA" : [], "lsu_rRNA" : [], "prok_marker" : [], "bac_marker" : [], "arc_marker" : [], "totalprots" : []}
		if index > len(outlist)-1:
			direction = -1
			index = len(outlist) -1
		if index < 0:
			direction = 1
			index = 0
		outlist[index].append(record)
		index += direction
		
	outlist = [ x for x in outlist if len(x) > 0 ] #removing any leftover fractions that did not get contigs (in case number of contigs was lower than number of fractions)
	sys.stderr.write("divided {} contigs into {} fractions\n".format(contigcounter, len(outlist)))

	if outfilebasename != None: #IF an outfilenasename is specified --> Do NOT return a list of lists of seqrecords, but instead write fractions fo tempfiles and return list of linemanes instead
		outfilenamelist = []
		for i in range(len(outlist)):
			outfilenamelist.append("{}_temp_fraction_{}.fasta".format(outfilebasename, i))
			with openfile(outfilenamelist[-1], "wt") as outfile:
				print("writing to {}".format(outfile.name))
				SeqIO.write(outlist[i], outfile, "fasta")
		return outfilenamelist, contigdict

	return outlist, contigdict 

def runprodigal(infasta, outfilename, prodigal="prodigal", threads = 1): #todo: allow piping via stdin (if input is a list: simply expect it to be a list of seqrecords)
	"""
	creates a protein fasta file of all CDS identified in the inputfasta via prodigal,
	all other prodigal output is ignored
	prodigal is called using the "-p meta" argument for metagenomes, in the assumption that the input fasta MAY consist of mutliple organisms
	The return value is simply the value of 'outfilename'
	"""
	prodigal_cmd = [prodigal, "-a", outfilename, "-p", "meta", "-q", "-o", "/dev/null"] # TODO: maybe add option to change translation table ("-g")? Although table 11 should be general enough?
	if type(infasta) == str and os.path.isfile(infasta):
		prodigal_cmd += ["-i", infasta]
		inputarg = None
	elif type(infasta) == list:
		inputarg =  "\n".join([record.format("fasta") for record in infasta])
	else:
		raise IOError("\nERROR: don't recognize query argument\n")
	try:
		prodigal_proc = subprocess.run(prodigal_cmd, input = inputarg, stdout = subprocess.PIPE, stderr = subprocess.PIPE, text = True)  
		prodigal_proc.check_returncode()
	except Exception:
		sys.stderr.write(prodigal_proc.stderr)
		raise Exception("\nERROR: Something went wrong while trying to call prodigal...\n")
	return outfilename
	
def runbarrnap_single(infasta, barrnap="barrnap", kingdom = "bac", threads=1): #todo: allow piping via stdin #todo instead of  the function doing a call for all kindoms at the same time, take a "kingdom" argument and do a seperate run for each kingdom --> allows better parallel multiprocessing!
	#tempfastalist, gffoutputs = [], []
	#for kingdom in ["bac", "arc", "euk"]:
	tempfasta = "temp_barrnap_{}.fasta".format(kingdom)
	barrnap_cmd = [barrnap, "--kingdom", kingdom, "--outseq", tempfasta, "--threads", str(threads), "-q", infasta] #todo: enable piping via stdin 
	assert os.path.isfile(infasta), "Error. can't find input file {}".format(infasta) # since barrnap can do multithreading, do not accet subdivided input_fasta for this
	try:
		barrnap_proc = subprocess.run(barrnap_cmd, stdout = subprocess.PIPE, stderr = subprocess.PIPE, text = True)
		barrnap_proc.check_returncode()
	except Exception:
		sys.stderr.write(barrnap_proc.stderr)
		raise Exception("\nERROR: Something went wrong while trying to call barrnap...\n")
	gff_output = barrnap_proc.stdout
	#todo: need to parse barrnap results from stdout (gff-output) rather than output-fasta-headers
	return (tempfasta, gff_output) #todo: make sure these results are then collected for each kingdom and run through deduplicate_barrnap_results()

def runbarrnap_all(infasta, outfilebasename, barrnap="barrnap", threads=3): #todo: split rRNA outputfiles by TYPE (16S & 23S)!. Also: return a contig2rRNA_dict in addition to the outfastaname!
	from Bio import SeqIO
	import misc
	joblist = []
	for kingdom in ["bac", "arc", "euk"]:
		joblist.append(("getmarkers", "runbarrnap_single", {"infasta" : infasta, "barrnap" : barrnap, "kingdom" : kingdom }))
	outputlist = misc.run_multiple_functions_parallel(joblist, threads)
	tempfasta_list = [op[0] for op in outputlist]
	gff_outputlist = [ op[1] for op in outputlist]
	final_fastadict, contig_rrnadict = deduplicate_barrnap_results(tempfasta_list, gff_outputlist) #todo: also get a dictionary which which markers are on which contig
	for ff in final_fastadict:
		outfile = openfile("{}_{}.fasta".format(outfilebasename, ff), "wt")
		SeqIO.write(final_fastadict[ff], outfile, "fasta")
		outfile.close()
		final_fastadict[ff] = outfile.name
	return final_fastadict, contig_rrnadict

def deduplicate_barrnap_results(tempfastas, gff_outputs):
	from Bio import SeqIO #todo: check if this is even necessary if SeqIO was already imported globally for this module
	import os
	contig_hit_dict = {}
	seqtype_dict = { "18S_rRNA": "ssu_rRNA", "16S_rRNA": "ssu_rRNA",  "23S_rRNA": "lsu_rRNA",  "23S_rRNA": "lsu_rRNA"}
	for gff in gff_outputs:
		for line in gff.rstrip().split("\n"):
			if line.startswith("#"):
				continue
			tokens = line.split()
			contig = tokens[0]
			start = int(tokens[3])
			stop = int(tokens[4])
			evalue = float(tokens[5])
			orient = tokens[6]
			rrna = tokens[8].split(";")[0][5:]
			seqid = "{}::{}:{}-{}({})".format(rrna, contig, start, stop, orient)
			altseqid = "{}::{}:{}-{}({})".format(rrna, contig, start-1, stop, orient) ##todo: remove this if barrnap issue is resolved. barrnap currently (v.0.9) gives different start position in fasta header and in gff output. Until i am sure what is the reason, or to make this work when if that is fixed in barrnap, i have to check for both variants
			if rrna == "5S_rRNA":
				continue #ignoring 5S rRNA for now
			if contig in contig_hit_dict:
				#todo find out if overlaps by more than 50%
				#if yes take only the one that has lower evalue
				#otherwise take both
				redundant = False
				index = 0
				while index < len(contig_hit_dict[contig]):
					evalueold = contig_hit_dict[contig][index]["evalue"]
					rangenew = set(range(start, stop))
					rangeold = range(contig_hit_dict[contig][index]["coords"][0], contig_hit_dict[contig][index]["coords"][1])
					intersection = rangenew.intersection(rangeold)
					if len(intersection)/min([len(rangenew), len(rangeold)]) > 0.5: #if it intersects by more than 50%, keep only the one with the better evalue
						if evalue < evalueold:
							contig_hit_dict[contig].pop(index)
							continue
						else:
							redundant = True
					index += 1
				if not redundant:
					contig_hit_dict[contig].append({"seqid" : seqid, "altseqid" : altseqid, "coords" : (start, stop, orient), "evalue" : evalue, "seqtype" : seqtype_dict[rrna] })	#todo: "altseqid" key not needed if barrnap issue is resolved				 
			else:
				contig_hit_dict[contig] = [{"seqid" : seqid, "altseqid" : altseqid, "coords" : (start, stop, orient), "evalue" : evalue, "seqtype" : seqtype_dict[rrna] }] #todo: "altseqid" key not needed if barrnap issue is resolved
	finalseqids = set()
	for contig in contig_hit_dict:
		for seq in contig_hit_dict[contig]:
			finalseqids.add(seq["seqid"])
			finalseqids.add(seq["altseqid"]) #todo: remove this if barrnap issue is resolved
	finalfastadict = {"ssu_rRNA" : [], "lsu_rRNA" : []} #16S & 18S are "ssu_rRNAs", 23S & 28S are "lsu_rRNAs". 5S is ignored
	#beforecounter = 0
	contig_rrna_dict = {}
	for fasta in tempfastas:
		#todo: add a seqcounter for before and after dedup
		infile = openfile(fasta)
		for record in SeqIO.parse(infile, "fasta"):
			recordtype, contig = parse_barrnap_headers(record.id)
			#recordtype = record.id[0:9] #todo: not the best way to get the marker type (16S_rRNA or 23S_rRNA) from the fasta-headers. but good enough for now...
			#beforecounter += 1
			if record.id in finalseqids: # todo:/note: I realize that if two models (e.g. arc & bac) detect the exact same region with the exact same coordinates, this would lead to a dupicate genesequence in the rRNA-predictions. However, currently it seems this would be without consequences for the further workflow
				finalfastadict[seqtype_dict[recordtype]].append(record)
				if contig not in contig_rrna_dict:
					contig_rrna_dict[contig] = {"ssu_rRNA" : [], "lsu_rRNA" : [] }
				contig_rrna_dict[contig][seqtype_dict[recordtype]].append({"seqid": record.id, "marker" : recordtype})
	for fasta in tempfastas: #currently doing this AFTER the previous loop, to make sure the files are only deleted when everything went well (debugging purposes)
		os.remove(fasta)

	return finalfastadict, contig_rrna_dict		#todo: also return a dictionary with contignames as keys and type of marker as values?
					
def parse_barrnap_headers(header):
	tokens = header.lstrip(">").split(":") #todo: add some kind of test to verify that this is actually barrnap-result-fasta-header
	recordtype = tokens[0]
	contig = tokens[2]
	return recordtype, contig

def runrnammer(infasta, outfilename, threads = 1): #todo: allow piping via stdin
	pass #todo: implement this (not a priority since rnammer is painful to install for most users)

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
	hmmsearch_cmd = [hmmsearch, "--noali", "--cpu", str(threads), "--domtblout", outfilename] + eval_cutoff_arg + score_cutoff_arg + [model]
	#print("\nquery = {}\n".format(query))
	if type(query) == str: #TODO: This assumes "if query is a string, it must be a filename." That is obviously BS! implement a check that tests if string is a fasta-record! #note to to self: for now i will assume fasta via stdin if query is a list of seqrecords 
		hmmsearch_cmd.append(query)
		inputarg = None
	elif type(query) == list:	#otherwise, if it is a list of seqrecords it must be something to pipe via stdin
		inputarg = "\n".join([record.format("fasta") for record in query])
	else:
		raise IOError("\nERROR: don't recognize query argument\n")
	try:
		hmmsearch_proc = subprocess.run(hmmsearch_cmd, input = inputarg, stdout = subprocess.PIPE, stderr = subprocess.PIPE, text = True)
		hmmsearch_proc.check_returncode()
	except Exception:# Todo: define/choose more detailed exception categories
		sys.stderr.write(hmmsearch_proc.stderr)
		raise Exception("\nERROR: something went wrong while trying to call hmmsearch...\n")
	return outfilename

### NOTE TO SELF: perform hmmsearch always with "sensitive" cutoff, and only PARSE hits with higher cutoffs --> enables reanalyses with different cutoffs without redoing hmmsearch!

def get_cutoff_dict(cutofffilename = cutofftablefile): #todo lookupfile with cutoffs for ALL used models. TODO: better: parse this from model.hmm files (require GA, TC and NC fields in all used models!)
	"""
	reads cutoff values from cutoff_file into a dictonary
	each model is represented as a seperate line with 4 columns:
		- first column = model name
		- second column = strict cutoff
		- third column = moderate cutoff
		- fourth column = sensitive cutoff"
	""" #also todo: make sure this is loaded only once for multiple input fastas (not reloaded again and again for each input)
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

def get_markerprotnames(proteinfastafile, cutoff_dict = cutofftablefile, hmmsearch = "hmmsearch", outdir = ".", cmode = "moderate", level = "prok", threads = "1"): #todo: turn list of markerdicts into dict of markerdits
	"""
	runs hmmersearch and and parse_hmmer on designated proteinfasta using models for designated level. Requires a cutoff_dict as returned by "get_gutoff_dict()"
	cutoff_dict should be a dictinary with the "strict", "moderate" and "sensitive" cutoff-values for each marker-model, but CAN also be a filename from which to parse that dict (default = parse from default file)
	returns a nested dictionary containing protein-identifiers as main keys and subdictionaries with respective marker-designations (key = "marker") and score values (key = "fscore")  as values for each hmm hit that passed cutoff criteria
	"""
	#print("\nget_markerprotnames()  --> proteinfastafile = {}\n".format(proteinfastafile))
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
			#print(hmmfile)
			#print(len(markerdict))
			#print(markerdict)
			#print("--------------------------")
		list_of_markerdicts.append(markerdict)
	return deduplicate_markerprots(list_of_markerdicts) #list_of_markerdicts will be in this order: [prok[, bact[, arch]]]

def deduplicate_markerprots(list_of_markerdicts): # For proteins with hits to different models, just keep the hit with the highest score. This function is a highly convoluted way to do this, but it is late and my brain is tired
	#todo: turn list of markerdicts into dict of markerdits
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
			
	

def __get_markerprotseqs(proteinfastafile, markerdict): #todo: implement piping proteinfastafile from stdin
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

def get_markerprots(proteinfastafile, cutoff_dict = cutofftablefile, cmode = "moderate", level = "prok", outfile_basename = "markerprots", threads = 1): #todo: turn list of markerdicts into dict of markerdits
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
	list_of_markerdicts = get_markerprotnames(proteinfastafile, cutoff_dict, hmmsearch = "hmmsearch", outdir = outdir, cmode = "moderate", level = level, threads = threads)
	outfilelist = []
	for i in range(len(list_of_markerdicts)):
		markerseqs = __get_markerprotseqs(proteinfastafile, list_of_markerdicts[i])
		outfilename = "{}_{}.faa".format(outfile_basename, levelorder[i])
		outfile = openfile(outfilename, "wt")
		SeqIO.write(markerseqs, outfile, "fasta")
		outfile.close()
		outfilelist.append(outfilename)
	return outfilelist
	
def write_markerdict(markerdict, outfilename):# todo: improve markerdict
	"""
	writes the marker dictionary, obtained by get_markerprotnames(), to an overview file in tab-seperated text-table format
	return value is simply the name of the outfile
	"""
	outfile = openfile(outfilename, "wt")
	for m in markerdict:
		outline = "{}\t{}\n".format(m, "\t".join([ str(markerdict[m][v]) for v in markerdict[m].keys() ]))
		outfile.write(outline)
	return outfilename

def combine_multiple_fastas(infastalist, outfilename = None, delete_original = True, contigdict = None): #pass contigdict in order to ba able to capture totalproteincounts per contig. currently only works for prodigal_output# todo: find a more flexible solution!
	"""
	different steps in getmarkers may subdivide input into fractions for better multiprocessing, and subsequently produce multiple output files
	This function is meant to combine such fastas to either a single output file (outfilename) or a list of seqrecords (if outfilename==None)
	Will delete the original fraction-fastas unless delete_original is set to False
	"""
	#todo: create an alternative version that writes to the outfile on the fly, for parsing huge assemblies
	import re
	from Bio import SeqIO
	recordcount = 0
	pattern = re.compile("_\d+$")
	outrecordlist=[]
	if outfilename != None:
		outfile = openfile(outfilename, "wt")
	for f in infastalist:
		infile=openfile(f)
		if outfilename!= None:
			for record in SeqIO.parse(f, "fasta"):
				recordcount += 1
				SeqIO.write([record], outfile, "fasta")
				if contigdict:
					contigname = re.sub(pattern, "", record.id)
					#print(contigdict[contigname].keys())
					contigdict[contigname]["totalprots"].append(record.id) #todo: if i understand python scopes correctly, te dictionary should be changed globally, even if not explicitely returned... check this!				
					contigdict[contigname]["totalprotcount"][0] += 1
					#print(contigdict[contigname]["totalprots"])
		else:
			outrecordlist.extend(list(SeqIO.parse(f, "fasta")))
			for record in outrecordlist:
				contigname = re.sub(pattern, "", record.id)
				contigdict[contigname]["totalprots"] += record.id #todo: if i understand python scopes correctly, te dictionary should be changed globally, even if not explicitely returned... check this!
		infile.close()
	if outfilename != None:
		outfile.close()
		output = outfilename
	else:
		output = outrecordlist
	if delete_original:
		for f in infastalist:
			os.remove(f)
	print("protein_recordcount = {}".format(recordcount))
	return output

def parse_protmarkerdict(protmarkerdict, contigdict, protmarkerlevel):
	import re
	pattern = re.compile("_\d+$")
	marker = protmarkerlevel_dict[protmarkerlevel]
	for protid in protmarkerdict:
		contigname = re.sub(pattern, "", protid)
		contigdict[contigname][marker].append({"seqid" : protid, "marker" : protmarkerdict[protid]["marker"]})
	return contigdict

def add_rrnamarker_to_contigdict(rrnamarkerdict, contigdict):
	for contig in rrnamarkerdict:
		#print(contigdict[contig])
		#print("-"*50)
		#print(rrnamarkerdict[contig])
		contigdict[contig].update(rrnamarkerdict[contig])
	return contigdict

def get_all_markers(infasta, outfilebasename, threads, cutoffdict = cutofftablefile): #todo: cutoffdict should be a dict optimally, but can be the name of a cutoffdict-file (from which to parse cutoffdict). also todo: for now running steps prodigal, barrnap etc one after another. but infuture find a way how to optimize paralellization
	#todo: this is probably going to be replaced by instantiation of a "bindata" object
	import misc
	subfastas, contigdict = split_fasta_for_parallelruns(infasta, minlength = 0, number_of_fractions = threads)
	commandlist = [("getmarkers", "runprodigal", {"infasta" : subfastas[i], "outfilename" : "test_prodigal_{}.faa".format(i) }) for i in range(len(subfastas))]
	tempprotfiles = misc.run_multiple_functions_parallel(commandlist, threads)
	totalprotfile = combine_multiple_fastas(tempprotfiles, outfilename = outfilebasename + "_totalprots.faa", delete_original = True, contigdict = contigdict) #todo: add support for gzipped protfiles in hmmsearch!
	protmarkerdictlist = get_markerprotnames(totalprotfile, cutoffdict, hmmsearch = "hmmsearch", outdir = ".", cmode = "moderate", level = "all", threads = "4")
	for pml in range(len(protmarkerdictlist)):
		contigdict = parse_protmarkerdict(protmarkerdictlist[pml], contigdict, pml)
	rRNA_fasta_dict, rrnamarkerdict = runbarrnap_all(infasta=infasta, outfilebasename=outfilebasename + "_rRNA", barrnap="barrnap", threads=threads)
	contigdict = add_rrnamarker_to_contigdict(rrnamarkerdict, contigdict)
	progressdump = { "ssu_rRNA_file" :  rRNA_fasta_dict["ssu_rRNA"], "ssu_rRNA_file" :  rRNA_fasta_dict["ssu_rRNA"], "totalprot_file" : totalprotfile, "protmarkerdictlist" : protmarkerdictlist, "contigdict" : contigdict }
	return progressdump

class bindata(object): #meant for gathering all contig/protein/marker info
	def __init__(self, contigfile, threads = 1, outbasedir = "mdmcleaner_results", mincontiglength = 0, cutofftable = cutofftablefile): #todo: enable init with additional precalculated infos
		self.binfastafile = contigfile
		bin_tempname = os.path.basename(contigfile)
		for suffix in [".gz", ".fa", ".fasta", ".fna", ".fas", ".fsa"]:
			if bin_tempname.endswith(suffix):
				bin_tempname = bin_tempname[:-len(suffix)]
		self.outbasedir = outbasedir		
		self.bin_tempname = bin_tempname
		self.bin_resultfolder = os.path.join(self.outbasedir, self.bin_tempname)
		for d in [self.outbasedir, self.bin_resultfolder]:
			if not os.path.exists(d):
				print("creating {}".format(d))
				os.mkdir(d)
		self._get_all_markers(threads, mincontiglength, cutofftable)
		#todo:finish this
	    #todo: simplify all those dicts
	def _get_all_markers(self, threads, mincontiglength, cutofftable): #todo: split into a.) get totalprots b.) get_markerprots c.) get rRNA genes!
		subfastas, self.contigdict = split_fasta_for_parallelruns(self.binfastafile, minlength = mincontiglength, number_of_fractions = threads)
		commandlist = [("getmarkers", "runprodigal", {"infasta" : subfastas[i], "outfilename" : os.path.join(self.bin_resultfolder, "tempfile_{}_prodigal_{}.faa".format(self.bin_tempname, i)) }) for i in range(len(subfastas))]
		tempprotfiles = misc.run_multiple_functions_parallel(commandlist, threads)
		self.totalprotfile = combine_multiple_fastas(tempprotfiles, outfilename = os.path.join(self.bin_resultfolder, self.bin_tempname + "_totalprots.faa"), delete_original = True, contigdict = self.contigdict)
		self.protmarkerdictlist = get_markerprotnames(self.totalprotfile, cutofftable, hmmsearch = "hmmsearch", outdir = self.bin_resultfolder, cmode = "moderate", level = "all", threads = threads) #todo: delete hmm_intermediate_results
		for pml in range(len(self.protmarkerdictlist)):
			self.contigdict = parse_protmarkerdict(self.protmarkerdictlist[pml], self.contigdict, pml)
		self.rRNA_fasta_dict, self.rrnamarkerdict = runbarrnap_all(infasta=self.binfastafile, outfilebasename=os.path.join(self.bin_resultfolder, self.bin_tempname + "_rRNA"), barrnap="barrnap", threads=threads) #todo add option for rnammer (using the subdivided fastafiles)?
		self.contigdict = add_rrnamarker_to_contigdict(self.rrnamarkerdict, self.contigdict)
		
		#todo: save progress as pickle of this the currenc instance of this class-object
	
	def _prep_contigsANDtotalprots():
		subfastas, self.contigdict = split_fasta_for_parallelruns(self.binfastafile, minlength = mincontiglength, number_of_fractions = threads)
		commandlist = [("getmarkers", "runprodigal", {"infasta" : subfastas[i], "outfilename" : os.path.join(self.bin_resultfolder, "tempfile_{}_prodigal_{}.faa".format(self.bin_tempname, i)) }) for i in range(len(subfastas))]
		tempprotfiles = misc.run_multiple_functions_parallel(commandlist, threads)
		self.totalprotfile = combine_multiple_fastas(tempprotfiles, outfilename = os.path.join(self.bin_resultfolder, self.bin_tempname + "_totalprots.faa"), delete_original = True, contigdict = self.contigdict)
	
	def _prep_protmarker():
		self.protmarkerdictlist = get_markerprotnames(self.totalprotfile, cutofftable, hmmsearch = "hmmsearch", outdir = self.bin_resultfolder, cmode = "moderate", level = "all", threads = "4") #todo: delete hmm_intermediate_results
		for pml in range(len(self.protmarkerdictlist)):
			self.contigdict = parse_protmarkerdict(self.protmarkerdictlist[pml], self.contigdict, pml)	

	def _prep_rRNAmarker():
		self.rRNA_fasta_dict, self.rrnamarkerdict = runbarrnap_all(infasta=self.binfastafile, outfilebasename=os.path.join(self.bin_resultfolder, self.bin_tempname + "_rRNA"), barrnap="barrnap", threads=threads) #todo add option for rnammer (using the subdivided fastafiles)?
		self.contigdict = add_rrnamarker_to_contigdict(self.rrnamarkerdict, self.contigdict)
				
	def _prep_onlycontigs():
		pass #todo: create a function that reads in contig fasta to intitiate contig_dct. For use in cases where ORF-calling was done on the complete metagenome 
	
	def get_contig_fastas(self):
		"""
		returns the bin contig-/scaffold-records in fasta format as a list
		"""
		with openfile(self.binfastafile) as infile:
			return list(SeqIO.parse(openfile, "fasta"))
	
	def pickleyourself(self):
		pass
		
	def unpickleyourself(self):
		pass	
		
	def get_contig2prot_dict():
		pass #todo: make this
	
	def get_prot2contig_dict():
		pass #todo: make this
	
	def get_prot2marker_dict():
		pass #todo: make this
		
	
		
######################################################
# test functions below (can be deleted)
def _test_markernames():
	sys.stderr.write("\ntesting get_markernames...")
	proteinfastafile = sys.argv[1]
	cutofftable = sys.argv[2]
	sys.stderr.write("\nreading cutofftable")
	cutoffdict = get_cutoff_dict(cutofftable)
	sys.stderr.write("\nsearching markers")
	markerdict = get_markerprotnames(proteinfastafile, cutoffdict, hmmsearch = "hmmsearch", outdir = ".", cmode = "moderate", level = "prok", threads = "4")
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
	outfastalist = get_markerprots(protfasta, cutoff_dict, level = level, outfile_basename = os.path.join(tempdir, "markers".format(level)), threads = 4)
	sys.stderr.write("  --> created files: '{}'".format(", ".join(outfastalist)))
	#todo: implement automatic blasts
	#todo implement actual lca

def _test_pipeline():
	import getdb
	infasta = sys.argv[1]
	threads = int(sys.argv[2])
	import misc
	outfilebasename = "testtesttest2"
	progressdump = get_all_markers(infasta, outfilebasename, threads, cutoffdict = cutofftablefile)
	getdb.dict2jsonfile(progressdump, "gelallmarkers.json")
	outfile = openfile("testcontigmarkers.tsv", "wt")
	contigdict = progressdump["contigdict"]
	sys.stderr.write("\nwriting results\n")
	outfile.write("contig\t{}\n".format("\t".join([x for x in contigdict[list(contigdict.keys())[0]]])))
	for contig in contigdict:
		line = "{}\t{}\n".format(contig, "\t".join([";".join([y["seqid"] if type(y) == dict else str(y) for y in contigdict[contig][x] ]) for x in contigdict[contig]])) #todo: in protmarkerdicts change "protid" to "seqid". Add "seqid" and "marker" keys to ssu and lsu entries
		outfile.write(line)

def _test_barrnap():
	from Bio import SeqIO
	infasta = sys.argv[1]
	threads = int(sys.argv[2])
	tempfilelist, gfflist = [], []
	rRNA_fasta = runbarrnap_all(infasta=infasta, outfilename="new_test_barrnap_results_dedup.fasta", barrnap="barrnap", threads=threads)

def _test_pipelineobj():
	import getdb
	infasta = sys.argv[1]
	threads = int(sys.argv[2])
	import misc
	outfilebasename = "testtesttest2"
	testbin = bindata(contigfile=infasta, threads=threads)
	outfile = openfile("testcontigmarkers.tsv", "wt")
	sys.stderr.write("\nwriting results\n")
	outfile.write("contig\t{}\n".format("\t".join([x for x in testbin.contigdict[list(testbin.contigdict.keys())[0]]])))
	for contig in testbin.contigdict:
		line = "{}\t{}\n".format(contig, "\t".join([";".join([y["seqid"] if type(y) == dict else str(y) for y in testbin.contigdict[contig][x] ]) for x in testbin.contigdict[contig]])) #todo: in protmarkerdicts change "protid" to "seqid". Add "seqid" and "marker" keys to ssu and lsu entries
		outfile.write(line)	

def _test_splitfasta2file():
	infasta = sys.argv[1]
	threads = int(sys.argv[2])
	outfilebasename = "huhudelmetest/fractiontest"
	a,b=split_fasta_for_parallelruns(infasta = infasta, number_of_fractions = threads, outfilebasename = outfilebasename)
	print(a)
def main():
	#_test_markernames()
	#_test_basicmarkers()
	#_test_multiprodigal()
	#_test_barrnap()
	#_test_pipeline()
	_test_pipelineobj()
	#_test_splitfasta2file()
if __name__ == '__main__':
	main()
