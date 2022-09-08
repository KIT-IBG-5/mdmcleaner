#!/usr/bin/env python
""" 
miscellaneous basic functions required by multiple modules 
"""

import os
import sys

def has_gzip_suffix(infilename):
	gzip_suffixes = [".gz", ".gzip"]
	for s in gzip_suffixes:
		if infilename.endswith(s):
			return True
	return False

def openfile(infilename, filemode = "rt"):
	""" a convenience function for creating file handles for compressed as well as uncompressed text files, for reading or writing as needed
	accepts a filename(required) and a filemode (default = "rt") argument.
	filemode may be any of ["w", "wt", "r", "rt"]
	automatically assumes gzip.compression if filename-suffix equals ".gz" 
	"""
	if has_gzip_suffix(infilename):
		import gzip
		filehandle = gzip.open(infilename, filemode)
	else:
		filehandle = open(infilename, filemode)
	return filehandle 

def read_fasta(infilename, mincontiglen=0):
	from Bio import SeqIO
	infile = openfile(infilename)
	return [record for record in SeqIO.parse(infile, "fasta") if len(record) >= mincontiglen ]

def write_fasta(outrecords, outfilename):
	from Bio import SeqIO
	# ~ import pdb; pdb.set_trace()
	with openfile(outfilename, "wt") as outfile: 
		SeqIO.write(outrecords, outfile, "fasta")
	
	
def unixzcat(*infilelist, outfilename): #my guess is, that this is probably much faster than doing it natively with python...
	import subprocess
	for i in infilelist:
		assert os.path.exists(i) and os.path.isfile(i), "\nError: \"{}\" doesn't exist or isn't a file".format(i)
		assert has_gzip_suffix(i), "\nError: trying to run unixzcat on a file without \".gz\"-suffix: {}\n".format(i)
	try:
		zcat_cmd = ["zcat"] + list(infilelist)
		# ~ print(" ".join(zcat_cmd))
		zcat_proc = subprocess.run(zcat_cmd, stdout = subprocess.PIPE, stderr = subprocess.PIPE, text = True)
		zcat_proc.check_returncode()
	except Exception:
		sys.stderr.write(zcat_proc.stderr)
		raise Exception("\nERROR: Something went wrong while trying to run unixscat with the following files: {}...\n".format(str(infilelist)))
	outfile=openfile(outfilename, "wt")
	outfile.write(zcat_proc.stdout)
	outfile.close()
	return outfilename

def untar(infilename, targetdir=".", filemode = None, removetar = False, verbose=False): #todo: make sure calling function passes verbose variable
	""" a convenience function for unpacking compressed and ancompressed tar files.
	accepts a filename(required) and an optional filemode (default = None) argument.
	filemode may be any of ["r:", "r:gz", "r:bz2", None ]. If filemode == None, it will try to determine filemode based on filename-extension
	".tar" = uncompressed tar --> filemode = "r:", ".tar.gz" = compressed tar --> filemode "r:gz"
	unpacks tar to targetdir and returns a list of unpacked filenames
	"""
	import tarfile
	import os
	import sys
	def track_progress(filelist):
		for f in range(len(filelist)):
			if (f % 50 == 0 or f == len(filelist)) and verbose:
				sys.stderr.write("\r\textracted {} of {} objects".format(f, len(filelist)))
				sys.stderr.flush()
			yield filelist[f]
		
	assert filemode in ["r:", "r:gz", "r:bz2", None], "\nERROR: filemode not allowed : '{}'\n".format(filemode)
	sys.stderr.write("\tassessing contents of '{}'\n".format(infilename))
	sys.stderr.flush() 
	if filemode == None:
		if infilename.endswith(".tar.gz"):
			filemode = "r:gz"
		elif infilename.endswith(".tar.bz2"):
			filemode = "r:bz2"
		elif infilename.endswith(".tar"):
			filemode = "r:"
		else:
			raise RunTimeError("\nError: can't guess filemode for '{}'. Please provide it!\n".format(infilename))
	infile = tarfile.open(infilename, filemode)
	contentlist = infile.getmembers()
	#print(contentlist)
	if verbose:
		sys.stderr.write("\tfound {} files in tar\n".format(len(contentlist)))
	sys.stderr.write("\textracting contents of '{}'\n".format(infilename))
	sys.stderr.flush() 
	infile.extractall(path=targetdir, members = track_progress(contentlist))
	sys.stderr.write("\n\tfinished extracting {}\n".format(infilename))
	infile.close()
	if removetar:
		sys.stderr.write("\tdeleting {}\n".format(infilename))
		os.remove(infilename)
	return [ os.path.join(targetdir, f.name) for f in contentlist ] #todo: add check if file or dir

def _run_any_function(modulename, functionname, arg_dict, threads=1):
	import importlib
	# ~ sys.stderr.write("\nimporting lib '{}'; function '{}'\n".format(modulename, functionname))
	module = importlib.import_module(modulename, "mdmcleaner")
	function = getattr(module, functionname)
	# ~ sys.stderr.write("\n" + str(function) + "\n")
	return function(**arg_dict, threads=threads) #all functions meant for multiprocessing must accept a threads argument (even if they are not really multiprocessing enabled!) and output an outputfilename

def run_multiple_functions_parallel(jobtuple_list, total_threads): #jobtuple_list should be a list of tuples with the following [(module1, function1, {kwargs.dict1}), (module2, function2, {kwargs.dict1}),...] 
	"""
	accepts a list of tuples, each containing the module and function as well as a dictionary of function-arguments for parallel jobs to run
	also accepts an integer, representing the total amount of threads to distribute ofer parallel jobs.
	the list of tuples should be set up similar to this example: [(module1, function1, {kwargs_dict1}), (module2, function2, {kwargs_dict2}), ...]
	each function MUST accept a "threads" argument, and each function should return a filename
	however, threads must NOT be set in the kwargs_dict. Instead they be distributed automatically, based on the number of jobs and the total available threads (passed as "total_threads")
	"total_threads" may be lower than the number of parallel jobs (jobs will be automatically queued)
	"""
	#TODO: allow option for UNEVEN distribution of threads over jobs. Figure aout a way how to do so... (e.g. protein-diamond vs nr should get much more threads than hmmer-searches!)
	def _distribute_threads_over_jobs(total_threads, num_jobs): # to distribute N threads over M groups as evenly as possible, put (N/M) +1 in (N mod M) groups, and (N/M) in the rest
		# ~ print("_distribute_threads_over_jobs")
		if num_jobs < total_threads:
			more_threads_list = [ int(total_threads / num_jobs) + 1 ] * (total_threads % num_jobs) #these should go to the lower priority markers, as there will be more of them
			fewer_threads_list = [ int(total_threads / num_jobs) ] * (num_jobs - len(more_threads_list))
			return fewer_threads_list + more_threads_list, num_jobs
		return [1] * num_jobs, total_threads

	# ~ print("run_multiple_functions_parallel")		
	from multiprocessing import Pool
	thread_args, no_processes = _distribute_threads_over_jobs(total_threads, len(jobtuple_list))
	arglist = []
	for i in range(len(jobtuple_list)):
		#print("{} + {}".format(tuple(job for job in jobtuple_list[i]), thread_args[i]))
		arglist.append(tuple(job for job in jobtuple_list[i]) + tuple([thread_args[i]]))
	jobpool = Pool(processes = no_processes)
	# ~ print(arglist)
	outfile_list = jobpool.starmap(_run_any_function, arglist)
	# ~ print("---huhu--")
	# ~ print(outfile_list)
	jobpool.close()
	jobpool.join()
	return outfile_list
	
def check_md5file(md5file):
	'''
	assumes filenames/paths in md5file are relative to location of md5file
	Returns a dictionary with filenames as keys and the expected and actual file-hashes as values as well as a bollean indicating a match (True) or mismatch (False)
	'''
	md5location = os.path.dirname(md5file)
	filecheckdict
	allfine = True
	with openfile(md5file) as infile:
		for line in infile:
			tokens=line.strip().split("\t")
			if len(tokens!=2):
				continue
			expected_md5 = tokens[0]
			filename = os.path.join(md5location, tokens[1])
			filehash = calculate_md5hash(filename)
			doesmatch = filehash == expected_md5
			filecheckdict[filename] = {"expected_hash" : expected_md5, "actual_hash" : filehash, "ok" : doesmatch}
	return filecheckdict

def calculate_md5hash(infile):
	"""
	calculate md5-hash of a file
	"""
	import hashlib #for comparing md5-checksums of downloaded databases
	blocksize = 2**20 #chunks to read file in
	with open(infile, "rb") as f:
		filehash = hashlib.md5()
		while True:
			data = f.read(blocksize)
			if not data:
				break
			filehash.update(data)
	return filehash.hexdigest()	

def from_json(jsonfilename):
	import json
	infile = openfile(jsonfilename)
	out_object = json.load(infile)
	infile.close()
	return out_object

def to_json(data_object, jsonfilename):
	import json
	outfile = openfile(jsonfilename, 'wt')
	json.dump(data_object, outfile)
	outfile.close()
	return jsonfilename

def is_emptyfile(filename):
	'''
	Returns True if a file is empty
	Returns False if a file is not empty
	Returns None if a file named <filename> does not exist or isnt a file
	''' 
	if os.path.exists(filename) and os.path.isfile(filename):
		if os.path.getsize(filename) == 0:
			return True
		return False
	
	
def to_pickle(data_object, picklefilename):
	try:
		import cPickle as pickle #this way the faster cPickle gets used IF available, but the slower standard pickle gets used otherwise
	except:
		import pickle
	f = open(picklefilename, 'wb')
	pickle.dump(data_object, f)
	f.close()	

def from_pickle(picklefilename):
	try:
		import cPickle as pickle #this way the faster cPickle gets used IF available, but the slower standard pickle gets used otherwise
	except:
		import pickle
	f = open(picklefilename, 'rb')
	tmp_dict = pickle.load(f)
	f.close()          
	return tmp_dict
	
if __name__ == '__main__':
	pass
