#!/usr/bin/env python
import sys
from collections import namedtuple
import misc
from misc import openfile

taxtuple = namedtuple("taxtuple", "taxid avident avscore")
"""
functions for obtaining least common ancestor (lca) annotations for a list of tax-ids
ALl these functions require an taxdb_object as input.
taxdb-objects already incude a basic strict pairwise lca function. This module provides additional functions for achieving weighted lca annotations for sets of taxids > 2.
"""

def strict_lca(taxdb, blasthitlist):
	"""
	Returns strict lca of a list of blasthits or pre-lca-annotations.
	Each blast hit should be represented as a named tuple with the following fields: (Accession, taxid, identity, score)
	"""
	assert len(blasthitlist) > 0, "\nError, but provide at least one blast hit!\n"
	interim_tax = blasthitlist[0].taxid
	#if len(blasthitlist) == 1: #probably already covered by looping over "range(1, len(blasthitlist))"
	#	return blasthitlist[0].taxid
	for i in range(1,len(blasthitlist):
		interim_tax = taxdb.get_strict_pairwise_lca(interim_tax, blasthitlist[i].taxid)
	interim_score = sum([bh.score for bh in blasthitlist])/len(blasthitlist)
	interim_id = sum([bh.identity for bh in blasthitlist])/len(blasthitlist)
	return taxtuple(taxid = interim_taxid, avident = interim_id, avscore = interim_score)  
	
def weighted_lca(taxdb, blasthitlist, cutoff = 0.95):
	#todo: maybe allow option to use identity as criterium, instead of score?
	"""
	Returns weighted lca of a list of blasthits or pre-lca-annotations.
	Each blast hit should be represented as a named tuple with the following fields: (Accession, taxid, identity, score)
	cutoff should be given as a fraction between 0.5 and 1. A taxon is accepted, as long as the sum of blast scroes supporting this assignment represents more than this fraction of the total sum of scores of all blast hits 
	It is not trecommended to use this function for a cutoff of 1 (strict lca)! For strict lca, please use the function "strict_lca" instead, which should be much faster!
	"""
	assert cutoff >= 0.5 and cutoff <=1, "\nError: cutoff {} not within allowed range (0.5-1.0)!\n"
	if cutoff == 1:
		sys.stderr.write("\nWARNING: it is NOT recommended to use 'weighted_lca' with a cutoff of 1! Try using 'strict_lca' instead!\n")
	#create tempdict
	#todo: dicitonary is NOT the best data structure for this (because deleting keys later on create whole new directories). consider creating a new class for this...
	#todo: is too convoluted! Simplify!
		## current method:
		## get taxonomy path for each hit
		## collect hits of same taxonomy in categories for each taxonomic level in taxpath. Record scores and identitiy-values
		## filter by summed score for each category
		## check if one category represents more that cutoff-fraction (default 95%) of total score, delete all others
		##  if there is only one to sart with -->that is unambigeous --> go a rank furher to determine taxid for next level
		##  if only one remains after filtering --> that is ambigeous but acceptable --> go a rank furher to determine taxid for next level but mark as "ambigeous == True"
		##  if two or more contradicting assignments remain --> DO not iterate further. return shared parent-assignment (current level) as LCA
	tempdict = {}
	for hit in blasthitlist:
		taxpath = [ x[1] for x in taxdb.taxid2taxpath(hit.taxid) ] #todo: add a taxdb-function "simplepath", that only returns the taxids as list. Should speed things up a little?
		parent = None
		for i in range(len(taxpath)):
			if not tempdict[i]:
				tempdict[i] = {}
			taxid = taxpath[i]
			if not taxid in tempdict[i]:
				tempdict[i][taxid] = {"scores" : [hit.score], "identities" : [hit.identity], "parent" : parent}
			else:
				assert parent == tempdict[i][taxid]["parent"], "\nError: conflicting parents for taxid '{}' : '{}' & '{}'. This means your taxonomy-system does not use unique TaxIDs!\n".format(taxid, parent, tempdict[i][taxid]["parent"])
				tempdict[i][taxpath]["scores"].append(hit.score)
				tempdict[i][identities].append[hit.identity)
			parent = taxid
	#evaluate tempdict
	outtaxpath = []
	parentblacklist = []
	taxassignemt = namedtuple("taxassignment", ["taxid", "average_score", "average_ident", "ambigeous"])
	ambigeous = False
	for i in tempdict:#todo: consider using a dedicated object,rather than dictionary?
		found_major_tax = False
		for tax in tempdict[i]
			if tempdict[i][tax]["parent"] in parentblacklist:
				parentblacklist.append(tax)
				del(tempdict[i][tax])
		if len(tempdict[i]) == 1:
			taxid = tempdict[i].keys()[0]
			outtaxpath = taxassignment(taxid, sum(tempdict[i][taxid]["scores"])/len(tempdict[i][taxid]["scores"]), sum(tempdict[i][taxid]["identities"])/len(tempdict[i][taxid]["identities"]), ambigeous)
			continue			
		ambigeous = True
		totalscoresum = sum( [ sum(tempdict[i][x]["scores"]) for x in tempdict[i] ] )
		for tax in tempdict[i]:
			if found_major_tax or sum(tempdict[i][tax]["scores"])/totalscoresum < cutoff::
				parentblacklist.append(tax)
			elif sum(tempdict[i][tax]["scores"])/totalscoresum >= cutoff:
				outtaxpath.append(tax)
				found_major_tax = True
		if not found_major_tax:
			break #if there are multiple contradicting taxonomic assignments, and no weighted major taxon can be determined based on cutoff, then stop here and return current taxon-level as LCA
	return tempdict #todo: only return last lca
		
