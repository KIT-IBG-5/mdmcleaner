#!/usr/bin/env python
import sys
from collections import namedtuple
import misc
from misc import openfile

taxlevels = ["root", "domain", "phylum", "class", "order", "family", "genus", "species"]
taxtuple = namedtuple("taxtuple", "seqid taxid identity score") #use exact same syntax as for input, to enable lca of lca annotations
"""
functions for obtaining least common ancestor (lca) annotations for a list of tax-ids
ALl these functions require an taxdb_object as input.
taxdb-objects already incude a basic strict pairwise lca function. This module provides additional functions for achieving weighted lca annotations for sets of taxids > 2.
"""

def contradicting_taxtuples(taxtuplelistA, taxtuplelistB, return_idents = False):
	"""
	checks if a path matches or not
	inputs are lists of taxtuples, as returned by weighted_lca
	if it does not: it returns the lowest level of mismatch ("e.g. "phylum")
	levels without taxassignment (None) in one of the taxtuplelists are ignored
	returns None if now contradiction is found
	optionally also returns a tuple of the respective average identitied of each assignment on that level (or None, None if no contradiction is found) 
	"""
	#todo: either alloe different kinds of inputs, or normalize lca/taxpath results (e.g. an lca-object)
	#todo: IMPORTANT: for now this assumes only the 7 major taxlevels! make this more flexible for input that doe snot fit this criteria e.g. contains "subfamily")
	if not None in [taxtuplelistA, taxtuplelistB]:
		maxlevel = min(len(taxtuplelistA), len(taxtuplelistB))
		for i in range(maxlevel):
			levelname = taxlevels[i]
			if taxtuplelistA[i].taxid != taxtuplelistB[i].taxid:
				# ~ print("contradiction! {} != {}".format(taxtuplelistA[i].taxid, taxtuplelistB[i].taxid))
				if return_idents:
					return levelname, (taxtuplelistA[i].average_ident, taxtuplelistB[i].average_ident)
				return levelname
	if return_idents:
		return None, None
	return None	

def contradict_taxtuble_taxpath(taxtuplelist, majortaxdict, return_idents = False): #todo: this is very convoluted. standardize taxpath, lca and majrtaxdict data types (create a taxobject or so...)
	if not None in [taxtuplelist, majortaxdict]:
		maxlevel = min(len(taxtuplelist), len(majortaxdict))
		for i in range(maxlevel):
			levelname = taxlevels[i]
			if taxtuplelist[i].taxid != majortaxdict[levelname][0][i]:
				# ~ print(taxtuplelist)
				# ~ print(len(taxtuplelist)
				# ~ print(majortaxdict[0])
				# ~ print(len(majortaxdict[0]))
				# ~ print("contradiction! {} != {}".format(taxtuplelist[i].taxid, majortaxdict[levelname][0][i]))
				if return_idents:
					return levelname, taxtuplelist[i].average_ident
				return levelname
	if return_idents:
		return None, None
	return None


def strict_lca(taxdb, seqid = None, blasthitlist=None, threads=1):
	"""
	Returns strict lca of a list of blasthits or pre-lca-annotations.
	Each blast hit should be represented as a named tuple with the following fields: (Accession, taxid, identity, score)
	"""
	assert blasthitlist != None and len(blasthitlist) > 0, "\nError, but provide at least one blast hit!\n"
	interim_taxid = blasthitlist[0].taxid
	#if len(blasthitlist) == 1: #probably already covered by looping over "range(1, len(blasthitlist))"
	#	return blasthitlist[0].taxid
	for i in range(1,len(blasthitlist)):
		interim_taxid = taxdb.get_strict_pairwise_lca(interim_taxid, blasthitlist[i].taxid)
	interim_score = sum([bh.score for bh in blasthitlist])/len(blasthitlist)
	interim_id = sum([bh.identity for bh in blasthitlist])/len(blasthitlist)
	return taxtuple(seqid = seqid, taxid = interim_taxid, identity = interim_id, score = interim_score) 
	# ~ return "fuck"
	
def weighted_lca(taxdb, seqid = None, blasthitlist=None, cutoff = 0.95, threads=1):
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
	#todo: dicitonary is probably NOT the best data structure for this (because deleting keys later on will create whole new dictionaries). consider creating a new class for this...
	#todo: is too convoluted! Simplify!
		## current method:
		## get taxonomy path for each hit
		## collect hits of same taxonomy in categories for each taxonomic level in taxpath. Record scores and identitiy-values
		## filter by summed score for each category
		## check if one category represents more that cutoff-fraction (default 95%) of total score, delete all others
		##  if there is only one to sart with -->that is unambigeous --> go a rank furher to determine taxid for next level
		##  if only one remains after filtering --> that is ambigeous but acceptable --> go a rank furher to determine taxid for next level but mark as "ambigeous == True"
		##  if two or more contradicting assignments remain --> DO NOT iterate further. return shared parent-assignment (current level) as LCA
	tempdict = {}
	for hit in blasthitlist:
		# ~ print("+"*30)
		# ~ print(hit)
		# ~ print("+"*30)
		taxpath = [ x[1] for x in taxdb.taxid2taxpath(hit.taxid) ] #todo: add a taxdb-function "simplepath", that only returns the taxids as list. Should speed things up a little?
		parent = None
		#print("*"*20)
		#print(tempdict)
		#print(taxpath)
		for i in range(len(taxpath)):
			if not i in tempdict:
				tempdict[i] = {}
			#print("  i={}".format(i))
			#print("  tempdict[i]={}\n*********".format(tempdict[i])) 
			taxid = taxpath[i]
			if not taxid in tempdict[i]:
				tempdict[i][taxid] = {"scores" : [hit.score], "identities" : [hit.identity], "parent" : parent}
			else:
				assert parent == tempdict[i][taxid]["parent"], "\nError: conflicting parents for taxid '{}' : '{}' & '{}'. This means your taxonomy-system does not use unique TaxIDs!\n".format(taxid, parent, tempdict[i][taxid]["parent"])
				tempdict[i][taxid]["scores"].append(hit.score)
				tempdict[i][taxid]["identities"].append(hit.identity)
			parent = taxid
	#evaluate tempdict
	outtaxpath = []
	parentblacklist = []
	taxassignment = namedtuple("taxassignment", ["taxid", "average_score", "average_ident", "ambigeous"]) #todo: replace this with something that is more similar to a "hit tuple". mayble add ambifeous field to those?
	ambigeous = False
	for i in tempdict:#todo: consider using a dedicated object,rather than dictionary?
		# ~ print("*"*60)
		# ~ print("i = {}".format(i))
		found_major_tax = False
		currentlevel_taxa = list(tempdict[i].keys())
		for tax in currentlevel_taxa:
			if tempdict[i][tax]["parent"] in parentblacklist:
				parentblacklist.append(tax)
				del(tempdict[i][tax])
		if len(tempdict[i]) == 1:
			taxid = list(tempdict[i].keys())[0]
			#print("only ONE lineage on level {} --> {}  --> {}".format(i, list(tempdict[i].keys()), taxid))
			taxass=taxassignment(taxid, sum(tempdict[i][taxid]["scores"])/len(tempdict[i][taxid]["scores"]), sum(tempdict[i][taxid]["identities"])/len(tempdict[i][taxid]["identities"]), ambigeous)
			# ~ print(list(tempdict[i].keys()))
			# ~ print("loop1")
			# ~ print("appending this to level {}: {}".format(i, taxid))
			# ~ print("-"*20)
			# ~ print(taxass)
			# ~ print("-"*20)
			outtaxpath.append(taxass)
			# ~ print(outtaxpath)
			continue			
		ambigeous = True
		# ~ print("MULTIPLE LINEAGES FOUND on level {}!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!".format(i))
		# ~ print(list(tempdict[i].keys()))
		totalscoresum = sum( [ sum(tempdict[i][x]["scores"]) for x in tempdict[i] ] )
		for tax in tempdict[i]:
			if found_major_tax or sum(tempdict[i][tax]["scores"])/totalscoresum < cutoff:
				parentblacklist.append(tax)
			elif sum(tempdict[i][tax]["scores"])/totalscoresum >= cutoff:
				taxass=taxassignment(tax, sum(tempdict[i][tax]["scores"])/len(tempdict[i][tax]["scores"]), sum(tempdict[i][tax]["identities"])/len(tempdict[i][tax]["identities"]), ambigeous)
				# ~ print("loop2")
				# ~ print("appending this to level {}: {}".format(i, tax))
				# ~ print("+"*20)
				# ~ print(taxass)
				# ~ print("+"*20)		
				outtaxpath.append(taxass)
				# ~ print(outtaxpath)
				found_major_tax = True
		if not found_major_tax:
			break #if there are multiple contradicting taxonomic assignments, and no weighted major taxon can be determined based on cutoff, then stop here and return current taxon-level as LCA
	return outtaxpath #, tempdict #todo: only return last lca


