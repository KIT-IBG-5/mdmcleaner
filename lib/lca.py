#!/usr/bin/env python
import sys
from collections import namedtuple
import misc
from misc import openfile

taxlevels = ["root", "domain", "phylum", "class", "order", "family", "genus", "species"]
taxtuple = namedtuple("taxtuple", "seqid taxid identity score") #use exact same syntax as for input, to enable lca of lca annotations
species_identity_cutoffs = {	"ssu_rRNA_tax" : 98, \
												"lsu_rRNA_tax" : 96, \
												"total_prots_tax" : 85, \
												"prok_marker_tax" : 90 } #based on AAI distribution observed in "Rodriguez-R, L. M., & Konstantinidis, K. T. (2014). Bypassing cultivation to identify bacterial species. Microbe, 9(3), 111-118." and "https://doi.org/10.1093/nar/gku169"

genus_identity_cutoffs = {	"ssu_rRNA_tax" : 92, \
												"lsu_rRNA_tax" : 85, \
												"total_prots_tax" : 55, \
												 "prok_marker_tax" : 60 } #based on AAI distribution observed in "Rodriguez-R, L. M., & Konstantinidis, K. T. (2014). Bypassing cultivation to identify bacterial species. Microbe, 9(3), 111-118." and "https://doi.org/10.1093/nar/gku169"

order_identity_cutoffs = {	"ssu_rRNA_tax" : 87, \
												"lsu_rRNA_tax" : 83, \
												"total_prots_tax" : 47, \
												 "prok_marker_tax" : 55 } #based on AAI distribution observed in "Rodriguez-R, L. M., & Konstantinidis, K. T. (2014). Bypassing cultivation to identify bacterial species. Microbe, 9(3), 111-118." and "https://doi.org/10.1093/nar/gku169"

phylum_identity_cutoffs = {	"ssu_rRNA_tax" : 77, \
												"lsu_rRNA_tax" : 73, \
												"total_prots_tax" : 40, \
												 "prok_marker_tax" : 45 } #based on AAI distribution observed in "Rodriguez-R, L. M., & Konstantinidis, K. T. (2014). Bypassing cultivation to identify bacterial species. Microbe, 9(3), 111-118." and "https://doi.org/10.1093/nar/gku169"

#todo: apply also orderlevel cutoffs (and domain level cutoffs
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
	#todo: either allow different kinds of inputs, or normalize lca/taxpath results (e.g. an lca-object)
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
		# ~ print("not None")
		maxlevel = min(len(taxtuplelist), len(majortaxdict))
		print(maxlevel)
		for i in range(maxlevel):
			# ~ print("{} = {}".format(i, taxtuplelist[i]))
			sys.stdout.flush()
			levelname = taxlevels[i]
			if taxtuplelist[i] is None or majortaxdict[levelname] is None: #todo: apparently should not occur. probably better to fix it in the barrnap- and majortaxdict functions of getmarkers.py, or if that fails in the previous maxlen-check (add a list comprehension that removes all entires with value None). this here is just a workaround...
				break
			# ~ print (taxtuplelist[i].taxid)
			# ~ print(levelname)
			# ~ print(majortaxdict[levelname])
			# ~ print(majortaxdict[levelname][0])
			# ~ print(majortaxdict[levelname][0][i])
						
			if taxtuplelist[i].taxid != majortaxdict[levelname][0][i]:
				# ~ print(taxtuplelist)
				# ~ print(len(taxtuplelist)
				# ~ print(majortaxdict[0])
				# ~ print(len(majortaxdict[0]))
				# ~ print("contradiction! {} != {}".format(taxtuplelist[i].taxid, majortaxdict[levelname][0][i]))
				print("returning {}, {}".format(levelname, taxtuplelist[i].average_ident))
				if return_idents:
					return levelname, taxtuplelist[i].average_ident
				return levelname
			sys.stdout.flush()
	# ~ print("returning None")
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
	
def weighted_lca(taxdb, seqid = None, blasthitlist=None, fractioncutoff = 0.95, taxlevel="total_prots_tax", threads=1):
	#todo: maybe allow option to use identity as criterium, instead of score?
	"""
	Returns weighted lca of a list of blasthits or pre-lca-annotations.
	Each blast hit should be represented as a named tuple with the following fields: (Accession, taxid, identity, score)
	fractioncutoff should be given as a fraction between 0.5 and 1. A taxon is accepted, as long as the sum of blast scroes supporting this assignment represents more than this fraction of the total sum of scores of all blast hits 
	It is not trecommended to use this function for a cutoff of 1 (strict lca)! For strict lca, please use the function "strict_lca" instead, which should be much faster!
	"""
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
	# ~ print("taxlevel {} --> appyling the following cutoffs: fractioncutoff: {}, species: {}, genus: {}".format(taxlevel, fractioncutoff, species_identity_cutoff, genus_identity_cutoff))
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
			if i == 7: #7 = species level
				# ~ print("this was assigned to species level: {}".format(hit))
				if hit.identity < species_identity_cutoff: #only assign up to species level, if identity is larger or equal to species_identity_cutoff (default 90% for proteins, 98% for rRNA)
					# ~ print("identity too low --> ignoring")
					break
				# ~ print("seems ok")
			if i == 6: #6 = genus level
				# ~ print("this was assigned to genus level: {}".format(hit))
				if hit.identity < genus_identity_cutoff: #only assign up to species level, if identity is larger or equal to species_identity_cutoff (default 90% for proteins, 98% for rRNA)
					# ~ print("identity too low --> ignoring")
					break
			if i == 4: #4 = order level
				# ~ print("this was assigned to order level: {}".format(hit))
				if hit.identity < order_identity_cutoff: #only assign up to species level, if identity is larger or equal to species_identity_cutoff (default 90% for proteins, 98% for rRNA)
					# ~ print("identity too low --> ignoring")
					break
			if i == 2: #2 = phylum level
				# ~ print("this was assigned to phylum level: {}".format(hit))
				if hit.identity < phylum_identity_cutoff: #only assign up to species level, if identity is larger or equal to species_identity_cutoff (default 90% for proteins, 98% for rRNA)
					# ~ print("identity too low --> ignoring")
					break		
				# ~ print("seems ok")				
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
			if found_major_tax or sum(tempdict[i][tax]["scores"])/totalscoresum < fractioncutoff:
				parentblacklist.append(tax)
			elif sum(tempdict[i][tax]["scores"])/totalscoresum >= fractioncutoff:
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
			break #if there are multiple contradicting taxonomic assignments, and no weighted major taxon can be determined based on fractioncutoff, then stop here and return current taxon-level as LCA
	#check if annotated to species level, and if that is the case, check if identity abov speciescutoff (default 90% for proteins)
	return outtaxpath #, tempdict #todo: only return last lca


