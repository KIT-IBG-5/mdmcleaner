#!/usr/bin/env python
""" 
creating, reading and parsing of blast files (Diamond and well as blast+).
only creates/handles tabular blast files in "-outfmt 6" format
"""

#assign blasts to contigs-proteins in a nother module.
#  --> in that module create "contigs" classes, that combine all blast hits (and corresponding marker-gene status) for each contig
#		--> in that module combine contig-objects of an assembly dataset into a dataset-object


from collections import namedtuple

from misc import openfile

def read_blast_tsv(infilename, max_evalue = None, min_ident = None):
	""" 
	returns a a list of names tuples, each representing the data in a blast line
	if max_evalue or min_ident are set, lines above or below these cutoffs are ignored
	each blast will later be assigned to a contig and to a "subject-type" (=stype), which can be either of ["16S", "23S", "univ_marker", "bact_marker", "arch_marker", "other"]
	"""
	#note to self: wanted to use namedtuples here, but namedtuples won't work well here, because i need them to be mutable.
	columns = ["query" : 0, "subject" : 1, "ident" : 2, "alignlen" : 3, "qstart": 6, "qend" : 7, "sstart" : 8, "ssend" : 9, "evalue" : 10, "score" : 11, "contig" : None, "stype" : None]
	#blastline = namedtuple("blastline", ["query", "subject", "ident", "alignlen", "qstart", "qend", "sstart", "ssend", "evalue", "score", "contig", "stype"], defaults = [None, None]) #probably won't need start end end coordinates tho...
	infile = openfile(infilename)
	blastlinelist = []
	for line in infile:
		tokens = line.strip().split("\t")
		#bl = blastline(tokens[0], tokens[1], float(tokens[2]), int(tokens[3]), int(tokens[6]), int(tokens[7]),int(tokens[8]), int(tokens[9], float(tokens[10], float(tokens[11]))
		bl = { x : tokens[columns[x]] if type(columns[x]) == int else None for x in columns}
		if max_evalue and bl["evalue"] > max_evalue: #bl.evalue > max_evalue:
			continue
		if min_ident and bl["ident"] < min_ident: #bl.ident < min_ident:
			continue
		blastlinelist.append(bl)
	infile.close()
	return blastlinelist

def contigs2blasthits(blastlinelist, parsetype = "prodigal", lookup_table = None)
	"""
	assigns contigs to blast hits
	parsetype must be one of ["prodigal", "rnammer", "lookup"]
	if parsetype is set to "lookup", a lookup_table must also be provided
	
	"""
	#start of subfunctions
	def _parse_prodigal(blastlinelist)
		import re
		pattern = re.compile("_\d+$")
		for blindex in len(blastlinelist):
			assert re.search(pattern, blastlinelist[blindex].query), "\nERROR: protein {} does not seem to be named according to prodigal standards. You should supply a lookup-table\n".format(columns[0])
			contig = re.sub(pattern, "", blastlinelist[blindex].query)
			blastlinelist
	
	#end of subfunctions
	assert parsetype in ["prodigal", "rnammer", "lookup"], "Do not recognize parsetype {}".format(parsetype))
	if parsetype == "lookup":
		assert lookup_table != None, "if parsetype is set to \"lookup\", a corresponding lookup_table must also be provided")
	if parsetype == prodigal:
		return _parse_prodigal(blastlinelist)

			
	
	
def stype2blasthits(blastlinelist, markerfiles):
	#still a fummy function right now
