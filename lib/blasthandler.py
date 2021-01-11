#!/usr/bin/env python
""" 
creating, reading and parsing of blast files (Diamond and well as blast+).
only creates/handles tabular blast files in "-outfmt 6" format
"""

#todo: assign blasts to contigs-proteins in a nother module.
#todo:   --> in that module create "contigs" classes, that combine all blast hits (and corresponding marker-gene status) for each contig
#todo:		--> in that module combine contig-objects of an assembly dataset into a dataset-object

#todo: use differend score cutoffs for protein/unknown-contig-blasts and 16S/23S-blasts (für proteins accept hits down to 50% of best score, for 16S/23S accept hits only down to 90-95%% of best score)

import sys
import os
assert sys.version_info >= (3, 7), "This module requires python 3.7+! You, however, are running '{}'".format(sys.version)
from collections import namedtuple

from misc import openfile, run_multiple_functions_parallel

def read_lookup_table(lookup_table, table_format = None): #should be able to read tsv, csv and gff. should output a simple dict with locus_tags as keys and contignames as values
	#consider putting this under "lookup_handler.py"? Probably naaaw...
	#start of subfunctions
	def _read_anything(lookup_table):
		global table_format
		linelist = []
		infile = openfile(lookup_table)
		firstline = True
		for line in infile:
			if line.strip().startswith("#"):
				if firstline and line.startswith("##gff-version"):
					table_format = "gff"
				if line == "##FASTA":
					table_format = "gff"
					break 
				continue
			linelist.append(line)
		return linelist
	
	def _guess_format(linelist):
		"""
		format is assumed to be tsv, if splitting by tab yields exactly two columns per line
		format is assumed to be csv, if splitting by comma yields exactly two columns per line
		otherwise format is assumed to be gff
		"""
		no_tab_fields = linelist[0].split("\t")
		no_comma_fields = linelist[0].split(",")
		if no_tab_fields == 9:
			return "gff"
		if no_tab_fields == 2:
			return "tsv"
		if no_tab_comma_fields == 2:
			return "csv"
		else:
			raise IOError("Can't recognize format of lookup_table! It should be either gff OR a simple tab- or comma seperated table with TWO columns \
			(locus_tag in the first, and contig_name in the second column)")
	
	def _parse_gff(linelist): #there may be some customized gffs out there (e.g. only "CDS"-lines grepped). Does not harm to include them also...
		lookup_dict = {}
		for line in linelist:
			tokens = line.split("\t")
			contig = tokens[0]
			attributes = tokens[9].split[";"]
			for att in attributes:
				if att.startswith("ID="):
					locustag = att[3:]
					break
			lookup_dict[locustag] = contig
		return lookup_dict
		
	def _parse_simpletable(linelist, seperator="\t"):
		lookup_dict = {}
		for line in linelist:
			line.split(seperator)[0] = locustag
			line.split(seperator)[0] = contig
			lookup_dict[locustag] = contig
		return lookup_dict

	#end of subfunctions
	
	possible_table_formats = ["gff", "tsv", "csv", None]
	assert table_format in possible_table_formats, "ERROR: tableformat must be one of {}".format(",".join([str(x) for x in possible_table_formats]))
	linelist = _read_anything(lookup_table)
	if table_format == None:
		table_format = _guess_format(linelist)
	if table_format == "gff":
		return _parse_gff(linelist)
	if table_format("tsv"):
		return _parse_simpletable("\t")
	if table_format("csv"):
		return _parse_simpletable(",")

def read_blast_tsv(infilename, max_evalue = None, min_ident = None):
	""" 
	returns a a list of dictionaries, each representing the data in a blast line
	if max_evalue or min_ident are set, lines above or below these cutoffs are ignored
	each blast will later be assigned to a contig and to a "subject-type" (=stype), which can be either of ["16S", "23S", "univ_marker", "proc_marker", "bact_marker", "arch_marker", "other"]
	"""
	#note to self: wanted to use namedtuples here, but namedtuples won't work well here, because i need them to be mutable.
	columninfos = {"query" : 0, "subject" : 1, "ident" : 2, "alignlen" : 3, "qstart": 6, "qend" : 7, "sstart" : 8, "ssend" : 9, "evalue" : 10, "score" : 11, "contig" : None, "stype" : None, "taxid": None} #Since python 3.7 all dicts are now ordered by default (YAY!). So the order of keys is maintained, without having to keep a seperate "orderlist". --> THEREFORE:  Ensure/Enforce python 3.7+ !!!
	infile = openfile(infilename)
	blastlinelist = []
	for line in infile:
		tokens = line.strip().split("\t")
		bl = { x : tokens[columninfos[x]] if type(columninfos[x]) == int else None for x in columninfos}
		if max_evalue and bl["evalue"] > max_evalue: #bl.evalue > max_evalue:
			continue
		if min_ident and bl["ident"] < min_ident: #bl.ident < min_ident:
			continue
		blastlinelist.append(bl)
	infile.close()
	return blastlinelist

def contigs2blasthits(blastlinelist, parsetype = "prodigal", lookup_table = None):
	"""
	assigns contigs to blast hits
	parsetype must be one of ["prodigal", "rnammer", "barrnap", "lookup"]
	if parsetype is set to "lookup", a lookup_table must also be provided
	"""
	#start of subfunctions
	def _parse_prodigal(blastlinelist):
		#sys.stderr.write("\nPARSING PRODIGAL STYLE\n")
		import re
		pattern = re.compile("_\d+$") #this pattern only works for standard naming scheme of standard prodigal output
		for blindex in range(len(blastlinelist)):
			assert re.search(pattern, blastlinelist[blindex]["query"]), "\nERROR: protein {} does not seem to be named according to prodigal standards. You should supply a lookup-table\n".format(blastlinelist[blindex]["query"])
			contig = re.sub(pattern, "", blastlinelist[blindex]["query"])
			blastlinelist[blindex]["contig"] = contig
			#sys.stderr.write("\n{} --> {}\n".format(blastlinelist[blindex]["query"], contig))
		return blastlinelist
	
	def _parse_lookup_table(blastlinelist, lookup_table):
		for blindex in range(len(blastlinelist)):
			contig = lookup_table[blastlinelist[blindex]["query"]]
			blastlinelist[blindex]["contig"] = contig
		return blastlinelist
		
	def _parse_rnammer(blastlinelist):
		import re
		patternprefix = re.compile("^\d+S_rRNA_") #this pattern only works for standard naming scheme of standard rnammer output
		patternsuffix = re.compile("_\d+-\d+_DIR[+-]$") #this pattern only works for standard naming scheme of standard rnammer output
		assert re.search(patternprefix, blastlinelist[blindex]["query"]) and re.search(patternsuffix, blastlinelist[blindex]["query"]), "\nERROR: Query {} does not seem to be named according to rnammer standards. You should supply a lookup-table\n".format(blastlinelist[blindex]["query"])
		for blindex in range(len(blastlinelist)):
			contig = re.sub(patternprefix, "", blastlinelist[blindex]["query"])
			contig = re.sub(patternsuffix, "", contig)
			blastlinelist[blindex]["contig"] = contig
		return blastlinelist
	
	def _parse_barrnap(blastlinelist):
		import re
		patternprefix = re.compile("^\d+S_rRNA::") #this pattern only works for standard naming scheme of standard barrnap output
		patternsuffix = re.compile(":\d+-\d+_\([+-]\)$") #this pattern only works for standard naming scheme of standard barrnap output
		assert re.search(patternprefix, blastlinelist[blindex]["query"]) and re.search(patternsuffix, blastlinelist[blindex]["query"]), "\nERROR: Query {} does not seem to be named according to barrnap standards. You should supply a lookup-table\n".format(blastlinelist[blindex]["query"])
		for blindex in range(len(blastlinelist)):
			contig = re.sub(patternprefix, "", blastlinelist[blindex]["query"])
			contig = re.sub(patternsuffix, "", contig)
			blastlinelist[blindex]["contig"] = contig
		return blastlinelist
		
	#end of subfunctions
	assert parsetype in ["prodigal", "rnammer", "barrnap", "lookup"], "Do not recognize parsetype {}".format(parsetype)
	if parsetype == "lookup":
		assert lookup_table != None, "if parsetype is set to \"lookup\", a corresponding lookup_table must also be provided"
		return _parse_lookup_table(blastlinelist, lookup_table)
	if parsetype == "prodigal":
		return _parse_prodigal(blastlinelist)
	if parsetype == "rnammer":
		return _parse_rnammer(blastlinelist)
	if parsetype == "barrnap":
		return _parse_barrnap(blastlinelist)
			
	
	
def stype2blasthits(blastlinelist, markersetdict): #markersetdict should be derived from a different module "fastahandler.py" that reads the markerfastas and returns either the sequences, or just the locus_tags/fasta_headers.
	#"markersetlist" should be a list of sets. sets should be ordered in decreasing hierarchy (16S first, then 23S then markerprots then totalprots)
	#each set contains the locus_tags of the markers used in the blast
	for blindex in range(len(blastlinelist)):
		for ms in markersetdict:
			if blastlinelist[bl]["query"] in markersetdict[ms]:
				blastlinelist[bl]["stype"] = ms
	return blastlinelist

def run_single_blastp(query, db, blast, outname, threads = 1):
#	from Bio.Blast.Applications import NcbiblastpCommandline
#	blastcmd = NcbiblastpCommandline(cmd = blast, query = query, db = db, evalue = 1e-10, out = outname + ".tmp", \
#								outfmt = 6, num_threads = threads)
	import subprocess
	blastcmd = subprocess.run([blast, "-query", query, "-db", db, "-evalue", "1e-10",\
							   "-outfmt", "6", "-num_threads", str(threads), "-out", outname + ".tmp"], \
							   stdout = subprocess.PIPE, stderr = subprocess.PIPE, text = True)
	try:
		blastcmd.check_returncode()
	except Exception:
		sys.stderr.write("\nAn error occured during blastp run with query '{}'\n".format(query))
		sys.stderr.write("{}\n".format(blastcmd.stderr))
		raise RuntimeError
	return outname

def run_single_blastn(query, db, blast, outname, threads = 1):
	# ~ from Bio.Blast.Applications import NcbiblastnCommandline
	# ~ blastcmd = NcbiblastnCommandline(cmd = blast, query = query, db = db, evalue = 1e-10, out = outname + ".tmp", \
								# ~ outfmt = 6, num_threads = threads)
	import subprocess
	blastcmd = subprocess.run([blast, "-query", query, "-db", db, "-evalue", "1e-10",\
							   "-outfmt", "6", "-num_threads", str(threads), "-out", outname + ".tmp"], \
							   stdout = subprocess.PIPE, stderr = subprocess.PIPE, text = True)
	try:
		blastcmd.check_returncode()
	except Exception:
		sys.stderr.write("\nAn error occured during blastn run with query '{}'\n".format(query))
		sys.stderr.write("{}\n".format(blastcmd.stderr))
		raise RuntimeError
	return outname
	
def run_single_diamondblastp(query, db, diamond, outname, threads = 1): #TODO: currently not setting "--tmpdir" & "--parallel-tmpdir" here! figure something out if this turns out to be problematic on hpc systems
	#TODO: add a maxmem arguemt that states how much memory can be used. use this to determine optimal blocksize and chunks for more efficient blasting. BLOCKSIZE=INT(MEMORY/6) CHUNKS=4/2/1 IF MEMORY >=12/24/48
	import subprocess
	blastcmd = subprocess.run([diamond, "--query", query, "--db", db, "--evalue", "1e-10",\
							   "--outfmt", "6", "--threads", str(threads), "--out", outname + ".tmp"], \
							   stdout = subprocess.PIPE, stderr = subprocess.PIPE, text = True)
	try:
		blastcmd.check_returncode()
	except Exception:
		sys.stderr.write("\nAn error occured during diamond blastp run with query '{}'\n".format(query))
		sys.stderr.write("{}\n".format(blastcmd.stderr))
		raise RuntimeError
	return outname
	
def _run_any_blast(query, db, path_appl, outname, threads):
	#TODO: fix stupid problem that blast (unfortunately) cannot handle gzipped query files
	appl = os.path.basename(path_appl)
	assert appl in ["blastp", "blastn", "diamond"], "\nError: unknown aligner '{}'\n".format(appl)
	_command_available(command=appl)
	if appl == "blastp":
		outname = run_single_blastp(query, db, appl, outname, threads)
	elif appl == "blastn":
		outname = run_single_blastn(query, db, appl, outname, threads)
	elif appl == "diamond":
		run_single_diamondblastp(query, db, appl, outname, threads)
	os.rename(outname + ".tmp", outname)
	return outname

def _distribute_threads_over_jobs(total_threads, num_jobs): # to distribute N threads over M groups as evenly as possible, put (N/M) +1 in (N mod M) groups, and (N/M) in the rest
	##TODO: test using misc.run_multiple_functions_parallel() for this instead! DELETE THIS IF MISC VERSION WORKS!
	if numjobs < total_threads:
		more_threads_list = [ (total_threads / num_jobs) + 1 ] * (total_threads % num_jobs) #these should go to the lower priority markers, as there will be more of them
		fewer_threads_list = [ (total_threads / num_jobs) ] * (num_jobs - len(more_threads_list))
		return fewer_threads_list + more_threads_list, num_jobs
	return [1] * num_jobs, total_threads

def run_multiple_blasts_parallel(basic_blastarg_list, outbasename, total_threads): #basic_blastarg_list = list of tuples such as [(query1, db1, blast1), (query2, db2, blast2),...])
	#TODO: test using misc.run_multiple_functions_parallel() for this instead! DELETE THIS IF MISC VERSION WORKS!
	if __name__ == '__main__':
		from multiprocessing import Pool
		thread_args, no_processes = _distribute_threads_over_jobs(total_threads, len(basic_blastarg_list))
		arglist = []
		for i in range(len(basic_blastarg_list)):
			outname = "{prefix}_{counter:03d}{appl}_{query}_vs_{db}.tsv".format(prefix = outbasename, counter = i, \
						appl = os.path.basename(basic_blastarg_list[i][2]), query = os.path.basename(basic_blastarg_list[i][0]), \
						db = os.path.basename(basic_blastarg_list[i][1]))
			arglist.append(tuple(bba for bba in basic_blastarg_list[i]) + (outname, thread_args[i]))
		masterblaster = Pool(processes = no_processes)
		outfile_list = masterblaster.starmap(_run_any_blast, arglist)
		masterblaster.close()
		masterblaster.join()
		return outfile_list

def make_diamond_db(infasta, outfilename, diamond = "diamond", threads = 1):
	import subprocess
	_command_available(command="diamond")
	mkdbcmd = [ diamond, "makedb", "--in", infasta, "--db", outfilename, "--threads", str(threads) ]
	mkdbproc = subprocess.run(mkdbcmd, stdout = subprocess.PIPE, stderr = subprocess.PIPE, text = True)
	sys.stderr.write("\n creating diamond-db \"{}\"...".format(outfilename))
	sys.stderr.flush()
	try:
		mkdbproc.check_returncode()
		sys.stderr.write("...finished!\n")
		sys.stderr.flush()
	except Exception:
		sys.stderr.write("\nAn error occured while creating '{}'\n".format(outfilename))
		sys.stderr.write("{}\n".format(mkdbproc.stderr))
		raise RuntimeError
	return outfilename	

def make_blast_db(infasta, outfilename, makeblastdb="makeblastdb", db_type="nucl", threads = 1): #threads  argument is ignored. Is only there to work with multiprocessing function
	import subprocess
	#todo: combine with "make_blast_db_from_gz" below: check if infasta ends with .gz. If yes: pipe from zcat, otherwise run makeblastdb directly
	assert db_type in ["nucl", "prot"],  "dbtype has to be either 'nucl' or 'prot'"
	sys.stderr.write("\n-->creating blastdb {}...\n".format(outfilename))
	sys.stderr.flush()
	mkdbcmd = [ makeblastdb, "-in", infasta, "-title", outfilename, "-out", outfilename, "-dbtype", db_type, "-max_file_sz", "4GB"] #hash_index apparently not possible for very large datasets
	print(" ".join(mkdbcmd))
	mkdbproc = subprocess.run(mkdbcmd, stdout = subprocess.PIPE, stderr = subprocess.PIPE, text = True)
	sys.stderr.write("\n creating blastdb-db \"{}\"...".format(outfilename))
	sys.stderr.flush()
	try:
		mkdbproc.check_returncode()
		sys.stderr.write("...finished!\n")
		sys.stderr.flush()
	except Exception:
		sys.stderr.write("\nAn error occured while creating '{}'\n".format(outfilename))
		sys.stderr.write("{}\n".format(mkdbproc.stderr))
		raise RuntimeError
	return outfilename	

def make_blast_db_from_gz(infasta, outfilename, makeblastdb="makeblastdb", db_type="nucl", threads = 1): #threads  argument is ignored. Is only there to work with multiprocessing function
	import subprocess
	_command_available(command="makeblastdb")
	# todo: if this works,incorporate it with "make_blast_db" above! check if infasta ends with .gz. If yes: pipe from zcat, otherwise run makeblastdb directly
	sys.stderr.write("\ncreating blastdb {}...\n".format(outfilename))
	sys.stderr.flush()
	if not infasta.endswith(".gz"):
		sys.stderr.write("\ninput is not compressed! using standard make_blast_db function".format(outfilename))
		return make_blast_db(infasta, outfilename, makeblastdb, db_type, threads) #todo: make this the other way round. make it call this function if it finds input to be compressed
	assert db_type in ["nucl", "prot"],  "dbtype has to be either 'nucl' or 'prot'"
	zcatcmd = [ "zcat", infasta]
	mkdbcmd = [ makeblastdb, "-title", outfilename, "-out", outfilename, "-dbtype", db_type, "-max_file_sz", "4GB"] #hash_index apparently not possible for very large datasets
	import time
	start = time.time()
	zcatproc = subprocess.Popen(zcatcmd, stdout = subprocess.PIPE, stderr = subprocess.PIPE, text = True)
	print(" ".join(mkdbcmd))
	sys.stdout.flush()
	mkdbproc = subprocess.Popen(mkdbcmd, stdin = zcatproc.stdout, stdout = subprocess.PIPE, stderr = subprocess.PIPE, text = True)
	stdout, stderr = mkdbproc.communicate()
	mkdbproc.wait()
	end = time.time()
	print("piping it this way took : {}".format(end - start))
	sys.stdout.flush()
	return outfilename
	#todo: compare how long it take to pipe this via python subprocess and ow long to do the same on the bash shell

def _command_available(command="diamond"): #todo: move this to misc, maybe?
	'''
	for checking if diamond or blast are actually installed in the current environment
	'''
	import subprocess
	try:
		subprocess.call([command],stderr=subprocess.PIPE, stdout=subprocess.PIPE)
	except FileNotFoundError:
		raise FileNotFoundError("\nError: can't run command '{}'. Maybe it is not installen in PATH, or you do not have the right conda environment activated?\n")
	else:
		return True
	return False
	
###### Testing functions below:
def test_prodigal_blasts():
	blastfile = sys.argv[1]
	sys.stderr.write("\nreading blast file \"{}\n".format(blastfile))
	blastlinelist = read_blast_tsv(blastfile)
	for bl in blastlinelist[:5]:
		#sys.stderr.write("\n{}\n".format(str(bl)))
		#sys.stderr.write("\n{}\n".format(type(bl)))
		sys.stderr.write("{}\n".format("\t".join([str(bl[x]) for x in bl])))
		sys.stderr.flush()
	sys.stderr.write("\n-------------------------\n")
	sys.stderr.write("\nadding contig_info:\n")
	blastlinelist = contigs2blasthits(blastlinelist)
	for bl in blastlinelist[:5]:
		#sys.stderr.write("\n{}\n".format(str(bl)))
		#sys.stderr.write("\n{}\n".format(type(bl)))
		sys.stderr.write("{}\n".format("\t".join([str(bl[x]) for x in bl])))
		sys.stderr.flush()
	sys.stderr.write("DONE (for now...")

def test_multiblastjobs():
	print("test_multiblastjobs")
	query_list = ["query1.fasta", "query2.fasta", "query3.fasta"]
	db = "/opt/db/blast/nr"
	tool = "blastp"
	jobtuplelist = []
	for query in query_list:
		jobtuplelist.append(("blasthandler", "_run_any_blast", {"query":query, "db":db, "path_appl": tool, "outname":"{}_{}_vs_{}.tab.tmp".format(tool, os.path.basename(query), os.path.basename(db))}))
	run_multiple_functions_parallel(jobtuplelist, 8)
	
def main():
	""" for testing purposes only"""
	#test_prodigal_blasts()
	print("main")
	test_multiblastjobs()

####### End of test functions

if __name__ == '__main__':
	main()



