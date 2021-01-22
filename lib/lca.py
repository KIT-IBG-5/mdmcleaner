#!/usr/bin/env python
import sys
import misc
from misc import openfile
"""
functions for obtaining least common ancestor (lca) annotations for a list of tax-ids
ALl these functions require an taxdb_object as input.
taxdb-objects already incude a basic strict pairwise lca function. This module provides additional functions for achieving weighted lca annotations for sets of taxids > 2.
"""


#todo: blastparser should return an dictionary of proteins, with protein-ids as keys and lists of hits as values.
#       --> each hit should be represented as a named tuple  with the following fields: (Accession, taxid, identity, score)!
def weighted_lca(taxdb, blasthitlist, cutoff = 0.95):
	#todo: maybe allow option to use identity as criterium, instead of score?
	"""
	Returns weighted lca of a list of blasthit.
	Each blast hit should be represented as a named tuple with the following fields: (Accession, taxid, identity, score)
	cutoff should be given as a fraction between 0.5 and 1. A taxon is accepted, as long as the sum of blast scroes supporting this assignment represents more than this fraction of the total sum of scores of all blast hits 
	It is not trecommended to use this function for a cutoff of 1 (strict lca)! For strict lca, please use the function "strict_lca" instead, which should be much faster!
	"""
	def nested_dict(tempdict, blasthitlist):
		if len(blasthitlist) == 0:
			return tempdict
		if blasthitlist[0] in tempdict:
			 tempdict[
	
	assert cutoff >= 0.5 and cutoff <=1, "\nError: cutoff {} not within allowed range (0.5-1.0)!\n"
	if cutoff == 1:
		sys.stderr.write
	#create tempdict
	tempdict = {}
	for hit in blasthitlist:
		taxpath = [ x[1] for x in taxdb.taxid2path(hit.taxid) ] #todo: add a taxdb-function "simplepath", that only returns the taxids as list. Should speed things up a little?
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
	for i in tempdict:
		## sort by sumemd score
		## keep only top 95% (based on score)
		##  if that is only one -->that is unambigeous --> go a rank furher to determine taxid for next level
		##  if that is two or more --> determine LCA of these. DO not iterate further
		
