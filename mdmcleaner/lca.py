#!/usr/bin/env python
import sys
from collections import namedtuple
from mdmcleaner import misc
from mdmcleaner.misc import openfile

taxlevels = ["root", "domain", "phylum", "class", "order", "family", "genus", "species"]
taxasstuple = namedtuple("taxasstuple", "seqid taxid identity score") #use exact same syntax as for input, to enable lca of lca annotations
species_identity_cutoffs = {	"ssu_rRNA_tax" : 98, \
												"lsu_rRNA_tax" : 96, \
												"totalprots_tax" : 85, \
												"prok_marker_tax" : 90 } #based on AAI distribution observed in "Rodriguez-R, L. M., & Konstantinidis, K. T. (2014). Bypassing cultivation to identify bacterial species. Microbe, 9(3), 111-118." and "https://doi.org/10.1093/nar/gku169"

genus_identity_cutoffs = {	"ssu_rRNA_tax" : 92, \
												"lsu_rRNA_tax" : 85, \
												"totalprots_tax" : 55, \
												 "prok_marker_tax" : 60 } #based on AAI distribution observed in "Rodriguez-R, L. M., & Konstantinidis, K. T. (2014). Bypassing cultivation to identify bacterial species. Microbe, 9(3), 111-118." and "https://doi.org/10.1093/nar/gku169"

order_identity_cutoffs = {	"ssu_rRNA_tax" : 87, \
												"lsu_rRNA_tax" : 83, \
												"totalprots_tax" : 47, \
												 "prok_marker_tax" : 55 } #based on AAI distribution observed in "Rodriguez-R, L. M., & Konstantinidis, K. T. (2014). Bypassing cultivation to identify bacterial species. Microbe, 9(3), 111-118." and "https://doi.org/10.1093/nar/gku169"

phylum_identity_cutoffs = {	"ssu_rRNA_tax" : 77, \
												"lsu_rRNA_tax" : 73, \
												"totalprots_tax" : 40, \
												 "prok_marker_tax" : 45 } #based on AAI distribution observed in "Rodriguez-R, L. M., & Konstantinidis, K. T. (2014). Bypassing cultivation to identify bacterial species. Microbe, 9(3), 111-118." and "https://doi.org/10.1093/nar/gku169"

#todo: apply also orderlevel cutoffs (and domain level cutoffs
"""
functions for obtaining least common ancestor (lca) annotations for a list of tax-ids
ALl these functions require an taxdb_object as input.
taxdb-objects already incude a basic strict pairwise lca function. This module provides additional functions for achieving weighted lca annotations for sets of taxids > 2.
"""

def contradicting_taxasstuples(taxasstuplelistA, taxasstuplelistB, return_idents = False):
	"""
	checks if a path matches or not
	inputs are lists of taxasstuples, as returned by weighted_lca
	if it does not: it returns the lowest level of mismatch ("e.g. "phylum")
	levels without taxassignment (None) in one of the taxasstuplelists are ignored
	returns None if now contradiction is found
	optionally also returns a tuple of the respective average identitied of each assignment on that level (or None, None if no contradiction is found) 
	"""
	#todo: either allow different kinds of inputs, or normalize lca/taxpath results (e.g. an lca-object)
	#todo: IMPORTANT: for now this assumes only the 7 major taxlevels! make this more flexible for input that does not fit this criteria e.g. contains "subfamily")
	if not None in [taxasstuplelistA, taxasstuplelistB]:
		maxlevel = min(len(taxasstuplelistA), len(taxasstuplelistB))
		for i in range(maxlevel):
			levelname = taxlevels[i]
			if taxasstuplelistA[i].taxid != taxasstuplelistB[i].taxid:
				# ~ print("contradiction! {} != {}".format(taxasstuplelistA[i].taxid, taxasstuplelistB[i].taxid))
				if return_idents:
					return levelname, (taxasstuplelistA[i].average_ident, taxasstuplelistB[i].average_ident)
				return levelname
	if return_idents:
		return None, None
	return None	
	
def contradict_taxasstuple_majortaxdict(taxasstuplelist, majortaxdict, return_idents = False): #todo: this is very convoluted. standardize taxpath, lca and majrtaxdict data types (create a taxobject or so...)
	"""
	Checks for contradictions between an input list of taxasstuples (taxasstuples are e.g. created by strict_lca())
	Returns the contradicting taxlevel and optionally also the average identity at that level . Returns None if no contradiction is found
	"""
	if not None in [taxasstuplelist, majortaxdict]:

		maxlevel = min(len(taxasstuplelist), len(majortaxdict))

		for i in range(maxlevel):

			levelname = taxlevels[i]
			if taxasstuplelist[i] is None or majortaxdict[levelname] is None: #todo: apparently should not occur. probably better to fix it in the barrnap- and majortaxdict functions of getmarkers.py, or if that fails in the previous maxlen-check (add a list comprehension that removes all entires with value None). this here is just a workaround...
				break
						
			if taxasstuplelist[i].taxid != majortaxdict[levelname][0][i]:
				if return_idents:
					return levelname, taxasstuplelist[i].average_ident
				return levelname
			sys.stdout.flush()
			
	if return_idents:
		return None, None
	return None


def strict_lca(taxdb, seqid = None, blasthitlist=None, threads=1):
	"""
	Returns strict lca of a list of blasthits or pre-lca-annotations.
	Each blast hit should be represented as a named tuple with the following fields: (Accession/seqid, taxid, identity, score)
	"""
	assert blasthitlist != None and len(blasthitlist) > 0, "\nError, but provide at least one blast hit!\n"
	# ~ print("doing LCA!")
	interim_taxid = blasthitlist[0].taxid
	#if len(blasthitlist) == 1: #probably already covered by looping over "range(1, len(blasthitlist))"
	#	return blasthitlist[0].taxid
	for i in range(1,len(blasthitlist)):
		# ~ print(interim_taxid)
		interim_taxid = taxdb.get_strict_pairwise_lca(interim_taxid, blasthitlist[i].taxid)
	interim_score = sum([bh.score for bh in blasthitlist])/len(blasthitlist)
	interim_id = sum([bh.identity for bh in blasthitlist])/len(blasthitlist)
	# ~ print(interim_taxid)
	# ~ print("============")
	# ~ import pdb; pdb.set_trace()
	return taxasstuple(seqid = seqid, taxid = interim_taxid, identity = interim_id, score = interim_score) 
	# ~ return "fuck"
	
def weighted_lca(taxdb, seqid = None, blasthitlist=None, fractioncutoff = 0.95, taxlevel="totalprots_tax", threads=1, return_contradicting_top2 = False):
	#todo: maybe allow option to use identity as criterium, instead of score?
	"""
	Returns weighted lca of a list of blasthits or pre-lca-annotations.
	Each blast hit should be represented as a named tuple with the following fields: (Accession, taxid, identity, score)
	fractioncutoff should be given as a fraction between 0.5 and 1. A taxon is accepted, as long as the sum of blast scroes supporting this assignment represents more than this fraction of the total sum of scores of all blast hits 
	It is not trecommended to use this function for a cutoff of 1 (strict lca)! For strict lca, please use the function "strict_lca" instead, which should be much faster!
	"""
	top2_contras, top2_contras_avidents, top2_contras_avscores  = None, None, None
	assert fractioncutoff >= 0.5 and fractioncutoff <=1, "\nError: rangecutoff {} not within allowed range (0.5-1.0)!\n"
	if fractioncutoff == 1:
		sys.stderr.write("\nWARNING: it is NOT recommended to use 'weighted_lca' with a fractioncutoff of 1! Try using 'strict_lca' instead!\n")
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
	species_identity_cutoff = species_identity_cutoffs[taxlevel]
	genus_identity_cutoff = genus_identity_cutoffs[taxlevel]
	order_identity_cutoff = order_identity_cutoffs[taxlevel]
	phylum_identity_cutoff = phylum_identity_cutoffs[taxlevel]
	
	tempdict = {}
	for hit in blasthitlist:

		taxpath = [ x[1] for x in taxdb.taxid2taxpath(hit.taxid) ] #todo: add a taxdb-function "simplepath", that only returns the taxids as list. Should speed things up a little?
		parent = None

		for i in range(len(taxpath)):
			if i == 7: #7 = species level
				if hit.identity < species_identity_cutoff: #only assign up to species level, if identity is larger or equal to species_identity_cutoff (default 90% for proteins, 98% for rRNA)
					break
			if i == 6: #6 = genus level
				if hit.identity < genus_identity_cutoff: #only assign up to species level, if identity is larger or equal to genus_identity_cutoff (default 60% for proteins, 92% for rRNA)
					break
			if i == 4: #4 = order level
				if hit.identity < order_identity_cutoff: #only assign up to species level, if identity is larger or equal to order_identity_cutoff (default 55% for proteins, 87% for rRNA)
					break
			if i == 2: #2 = phylum level
				if hit.identity < phylum_identity_cutoff: #only assign up to species level, if identity is larger or equal to species_identity_cutoff (default 45% for proteins, 77% for rRNA)
					break					
			if not i in tempdict:
				tempdict[i] = {}
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
	found_major_tax = False
	for i in tempdict:#todo: consider using a dedicated object,rather than dictionary?
		currentlevel_taxa = list(tempdict[i].keys())
		for tax in currentlevel_taxa:
			if tempdict[i][tax]["parent"] in parentblacklist:
				parentblacklist.append(tax)
				del(tempdict[i][tax])
		if len(tempdict[i]) == 1:
			taxid = list(tempdict[i].keys())[0]
			taxass=taxassignment(taxid, sum(tempdict[i][taxid]["scores"])/len(tempdict[i][taxid]["scores"]), sum(tempdict[i][taxid]["identities"])/len(tempdict[i][taxid]["identities"]), ambigeous)
			outtaxpath.append(taxass)
			continue			
		ambigeous = True
		totalscoresum = sum( [ sum(tempdict[i][x]["scores"]) for x in tempdict[i] ] )
		for tax in tempdict[i]:
			if found_major_tax or sum(tempdict[i][tax]["scores"])/totalscoresum < fractioncutoff:
				parentblacklist.append(tax)
			elif sum(tempdict[i][tax]["scores"])/totalscoresum >= fractioncutoff:
				taxass=taxassignment(tax, sum(tempdict[i][tax]["scores"])/len(tempdict[i][tax]["scores"]), sum(tempdict[i][tax]["identities"])/len(tempdict[i][tax]["identities"]), ambigeous)		
				outtaxpath.append(taxass)
				found_major_tax = True

		if not found_major_tax:
			sorted_taxoptions = sorted(tempdict[i].items(), key = lambda x:-sum(x[1]["scores"])/len(x[1]["scores"]))
			top2_contras = [ t[0] for t in sorted_taxoptions[:2]] #dict.items() returns a list of key/value tuples
			top2_contras_avidents = [ sum(t[1]["identities"])/len(t[1]["identities"]) for t in sorted_taxoptions[:2] ]
			top2_contras_avscores = [ sum(t[1]["scores"])/len(t[1]["scores"]) for t in sorted_taxoptions[:2] ]
			break #if there are multiple contradicting taxonomic assignments, and no weighted major taxon can be determined based on fractioncutoff, then stop here and return current taxon-level as LCA
	#check if annotated to species level, and if that is the case, check if identity abov speciescutoff (default 90% for proteins)
	if return_contradicting_top2:
		return outtaxpath, top2_contras, top2_contras_avidents, top2_contras_avscores
	return outtaxpath #, tempdict #todo: only return last lca


