#!/usr/bin/env python
""" 
miscellaneous basic functions required by multiple modules 
"""
def openfile(infilename, filemode = "rt"):
	""" a convenience function that will be moved to another more general module later"""
	if infilename.endswith(".gz"):
		import gzip
		filehandle = gzip.open(infilename, filemode)
	else:
		filehandle = open(infilename, filemode)
	return filehandle 
