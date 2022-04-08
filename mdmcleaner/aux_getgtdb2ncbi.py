#!/usr/bin/env python

import sys
import argparse

from mdmcleaner import misc

myparser = argparse.ArgumentParser("auxilliary script to get gtdb2ncbi taxonomy lookuptable. currently requires tabstopped textfiles for each tab in the gtdb lookup xlsx-tables")
myparser.add_argument("-b", "--bacterial", action = "store", dest = "bacterial", nargs = "+", help = "bacterial lookup tables. should be listed in hierarcical order, with phyum first and species last")
myparser.add_argument("-a", "--archaeal", action = "store", dest = "archaeal", nargs = "+", help = "archaeal lookup tables. should be listed in hierarcical order, with phyum first and species last")
args = myparser.parse_args()

"""
actual plan:
	read in table
	for each gtdb taxon, gather the list of possible ncbi taxa with
	if the bestnacbi taxa matches in
"""

def lca(taxlist):
	#todo: get ncbi accession for each of these. then use getlca function of getdb to get lca
	pass

def read_table(tablename): #todo: maybe keep level info ("g__", "s__" etc)? to make it wasier to look up the exact correct ncbi taxon for taxa with unfortatane double names?
	infile = openfile(tablename)
	lookupdict = {}
	for line in infile:
		if line.startswith("GTDB ") or "Number of genomes" in line:
			continue
		tokens = line.strip().split("\t")
		gtdbtax = tokens[0]
		ncbitokens = tokens[3].replace("%","").split()
		possncbitax = [ (t[i].split("(")[0][3:]), float(t[t+1])) for t in range(0, len(ncbitokens), 2) ] # for each tuple in list, first entry is ncbi-taxname, second is percentage
		if get_ncbitaxid(gtdbtax) == None: #first look if gtdb name exists in ncbi
			if possncbitax[0] == "":
				if len(tokens) == 5 and tokens[4] not in ["-", ""]: 
					ncbitax = tokens[4] #(custom "alternative" column by me, may or may not exist)
				else:
					ncbitax = lca([ tax[0] for tax in possncbitax if (tax[1] > 2 and tax[0] != "" ]) #attempt to do lca classification of all possible included taxa (up to a certain representation, set here arbitrarily to 2%)
			ncbitaxid = get_ncbitaxid(ncbitax, level)
			
		else:
			ncbitaxid = get_ncbitaxid(gtdbtax, level)
		lookupdict[gtdbtax] = {"ncbitaxid" : ncbitaxid, "ncbitax" : ncbitax}
		##plan: first look if gtdb name exists in ncbi
		##   yes: take that (ignore lookup as that may be outdated or based on wrong classifications. this is nly meant or rRNA classification anyway)
	return lookupdict
			

		
