#!/usr/bin/env python2
from Bio import SeqIO
import argparse
import os.path
import random

def openfasta(fastaname, mode = "rt"):
	import gzip
	if fastaname.endswith(".gz") or fastaname.endswith(".gzip"):
		infile = gzip.open(fastaname, mode)
	else:
		infile = open(fastaname, mode)
	return infile

def readfasta(fastaname):
	infile = openfasta(fastaname)
	incontigs = list(SeqIO.parse(infile, "fasta"))
	infile.close()
	return incontigs

def readwrite_fasta(fastaname, fractions, outputbasename, nonrandom = False):
	incontigs = readfasta(fastaname)
	contigs_per_fraction = len(incontigs)/fractions
	outfilelist = []
	while contigs_per_fraction <= 0:#if there are too few contigs for specified number of fractions --> just make fewer fractions!
		fractions -= 1
		contigs_per_fraction = len(incontigs)/fractions
	for fileindex in range(fractions):
		outfilename = "{}_{}.fasta".format(outputbasename, fileindex + 1)
		outfile = open("{}_{}.fasta".format(outfilename, "w")
		outseqs = []
		while len(incontigs) > 0:
			if nonrandom:
				contindex = 0
			else:
				contindex = random.randint(0,len(incontigs)-1)
			if fileindex < (args.fractions -1):
				if len(outseqs) >= contigs_per_fraction:
					break
			outseqs.append(incontigs.pop(contindex))
		SeqIO.write(outseqs, outfile, "fasta")
		outfile.close()
		outfilelist.append(outfilename)
	return outfilelist

def subdivide(input_fasta, outfileprefix, fractions, nonrandom = False):
	for fasta in input_fasta:
		if args.outfilename == None:
			outfileprefix = fasta
		readwrite_fasta(fasta, fractions, outfileprefix, nonrandom)
	print "done!"

if __name__ == '__main__':
	parser = argparse.ArgumentParser(description = "subdivide contigs in multifastas into specified number of fractions (with random subsampling, but optionally original order can also be kept)")
	parser.add_argument("input_fasta", nargs = "+", help = "Input file(s)")
	parser.add_argument('-nf', '--number_fractions', action = "store", dest = "fractions", type = int, default = 8, help = "Number of fractions to subdivide fasta into")
	parser.add_argument('-o', '--output', action = "store", dest = "outfilename", default = None, help = "basename for output-fastas. index-numbers will be appendedto the names. inputfilename will be used if this argument is left empty")
	parser.add_argument("-nr", '--nonrandom', action = "store_true", dest = "nonrandom", default = False, help = "subdivide contigs in nonrandom order --> will be ordered the same as in the input (default = False)")
	args = parser.parse_args()
	input_fasta = args.input_fasta
	outfileprefix = args.outfilename
	fractions = args.fractions
	subdivide(args.input_fasta, args.outfilename, args.fractions, args.nonrandom)
