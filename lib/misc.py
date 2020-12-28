#!/usr/bin/env python
""" 
miscellaneous basic functions required by multiple modules 
"""
def openfile(infilename, filemode = "rt"):
	""" a convenience function for creating file handles for compressed as well as uncompressed text files, for reading or writing as needed
	accepts a filename(required) and a filemode (default = "rt") argument.
	filemode may be any of ["w", "wt", "r", "rt"]
	automatically assumes gzip.compression if filename-suffix equals ".gz" 
	"""
	if infilename.endswith(".gz"):
		import gzip
		filehandle = gzip.open(infilename, filemode)
	else:
		filehandle = open(infilename, filemode)
	return filehandle 

def unixzcat(infilelist, outfilename): #my guess is, that this is probably much faster than doing it natively with python...
	pass
	#call zcat to concatenate all files in infilelist to outfilename
	#check if successful
	#delete files in infilelist
	#compress outfile
	#return outfilename

def untar(infilename, targetdir=".", filemode = None): #todo: add delete option to misc.unar(). make it delete tar file when finished
	""" a convenience function for unpacking compressed and ancompressed tar files.
	accepts a filename(required) and an optional filemode (default = None) argument.
	filemode may be any of ["r:", "r:gz", None ]. If filemode == None, it will try to determine filemode based on filename-extension
	".tar" = uncompressed tar --> filemode = "r:", ".tar.gz" = compressed tar --> filemode "r:gz"
	unpacks tar to targetdir and returns a list of unpacked filenames
	"""
	import tarfile
	import os
	import sys
	def track_progress(filelist):
		for f in range(len(filelist)):
			if f % 50 == 0:
				sys.stderr.write("\r extracted {} of {} objects".format(f, len(filelist)))
				sys.stderr.flush()
			yield filelist[f]
		
	assert filemode in ["r:", "r:gz", None], "\nERROR: filemode not allowed : '{}'\n".format(filemode)
	sys.stderr.write("    assessing contents of '{}'\n".format(infilename))
	sys.stderr.flush() 
	if filemode == None:
		if infilename.endswith(".tar.gz"):
			filemode = "r:gz"
		elif infilename.endswith(".tar"):
			filemode = "r:"
		else:
			raise RunTimeError("\nError: can't guess filemode for '{}'. Please provide it!\n".format(infilename))
	infile = tarfile.open(infilename, filemode)
	contentlist = infile.getmembers()
	#print(contentlist)
	print("found {} files in tar".format(len(contentlist)))
	sys.stderr.write("    extracting contents of '{}'\n".format(infilename))
	sys.stderr.flush() 
	sys.stdout.flush() # todo: only for debugging
	infile.extractall(path=targetdir, members = track_progress(contentlist))
	sys.stderr.write("\r finished extracting {}\n".format(infilename))
	infile.close()
	return [ os.path.join(targetdir, f.name) for f in contentlist ] #todo: add check if file or dir

def _run_any_function(modulename, functionname, arg_dict, threads=1):
	module = __import__(modulename)
	function = getattr(module, functionname)
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
		print("_distribute_threads_over_jobs")
		if num_jobs < total_threads:
			more_threads_list = [ int(total_threads / num_jobs) + 1 ] * (total_threads % num_jobs) #these should go to the lower priority markers, as there will be more of them
			fewer_threads_list = [ int(total_threads / num_jobs) ] * (num_jobs - len(more_threads_list))
			return fewer_threads_list + more_threads_list, num_jobs
		return [1] * num_jobs, total_threads

	print("run_multiple_functions_parallel")		
	from multiprocessing import Pool
	thread_args, no_processes = _distribute_threads_over_jobs(total_threads, len(jobtuple_list))
	arglist = []
	for i in range(len(jobtuple_list)):
		print("{} + {}".format(tuple(job for job in jobtuple_list[i]), thread_args[i]))
		arglist.append(tuple(job for job in jobtuple_list[i]) + tuple([thread_args[i]]))
	jobpool = Pool(processes = no_processes)
	outfile_list = jobpool.starmap(_run_any_function, arglist)
	jobpool.close()
	jobpool.join()
	return outfile_list
	
	
if __name__ == '__main__':
	pass
