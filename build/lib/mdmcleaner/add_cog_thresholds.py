#!/usr/bin/env python

import sys
import os
import argparse

from mdmcleaner.misc import openfile

myparser=argparse.ArgumentParser("adds threshold cutoffs to hmm files for easier hmmsearch")
myparser.add_argument("-cot", "--cutofftable", action = "store", dest = "cutofftable", default = None, help = "cutofftable")
myparser.add_argument("-def", "--deftable", action = "store", dest = "deftable", default = None, help = "'cog-20.def.tab' from cog ftp server at ncbi")
myparser.add_argument("inmodels", nargs = "+", help = "hmm model files")
args = myparser.parse_args()

def get_cutoff_dict(cutofffilename): #todo lookupfile with cutoffs for ALL used models. 
	"""
	reads cutoff values from cutoff_file into a dictonary
	each model is represented as a seperate line with 4 columns:
		- first column = model name
		- second column = strict cutoff
		- third column = moderate cutoff
		- fourth column = sensitive cutoff
	"""
	cutofffile = openfile(cutofffilename)
	cutoff_dict = {}
	for line in cutofffile:
		if line.startswith("#"):
			continue
		tokens = line.split()
		model = tokens[0]
		strict = float(tokens[1])
		moderate = float(tokens[2])
		sensitive = float(tokens[3])
		cutoff_dict[model] = {"strict" : strict, "moderate" : moderate, "sensitive" : sensitive}
	cutofffile.close()
	return cutoff_dict

def get_naming_dict(deftable):
	infile = openfile(deftable)
	naming_dict = {}
	for line in infile:
		tokens = line.strip("\n").split("\t")
		acc = tokens[0]
		definition = tokens[2]
		#print("--{}__".format(tokens))
		name = tokens[3]
		naming_dict[acc] = { "acc" : acc, "name" : name, "description" : definition }
	return naming_dict

def parse_modelfiles(cutoff_dict=None, namingdict=None):
	hmmfiles = args.inmodels
	for hmmfile in hmmfiles:
		infile = openfile(hmmfile)
		outfile = openfile("edited_" + hmmfile, "wt")
		for line in infile:
			#print(line)
			if line.startswith("NAME "):
				model = line.split()[1]
				if namingdict:
					name = namingdict[model]["name"]
					acc = namingdict[model]["acc"]
					description = namingdict[model]["description"]
					outfile.write("NAME  {}\n".format(name))
					outfile.write("ACC   {}\n".format(acc))
					outfile.write("DESC  {}\n".format(description))
					continue
			outfile.write(line)
			if cutoff_dict and line.startswith("CKSUM "):
				strict = cutoff_dict[model]["strict"]
				moderate = cutoff_dict[model]["moderate"]
				sensitive = cutoff_dict[model]["sensitive"]
				outfile.write("GA    {0:.2f} {0:.2f};\n".format(moderate))
				outfile.write("TC    {0:.2f} {0:.2f};\n".format(strict))
				outfile.write("NC    {0:.2f} {0:.2f};\n".format(sensitive))
	infile.close()
	outfile.close()


def main():
	cutoff_dict, namingdict = None, None
	if args.cutofftable:
		cutoff_dict = get_cutoff_dict(args.cutofftable)
	if args.deftable:
		namingdict = get_naming_dict(args.deftable)
	parse_modelfiles(cutoff_dict, namingdict)
	
	print("finished")

main()
