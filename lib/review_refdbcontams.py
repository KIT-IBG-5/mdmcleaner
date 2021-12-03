#!/usr/bin/env python
from Bio import SeqIO
'''
much of this will move to the other modules when finished
'''

def get_contig_from_blastdb(seqid, blastdb, blastdbcmd = "blastdbcmd"):
	'''
	this function should go directly to blasthandler.py
	a different functions should also be added to the databaseobject created at getdb, that then calls this blasthandler funciton for the correct blastdb
	seqid may be a single seqid or a list of seqids
	returns a list of Seqrecords
	'''
	import subprocess
	from io import StringIO
	if type(seqid) == list:
		seqid = ",".join(seqid)
	cmd_blastdbcmd = ["blastdbcmd", "-entry", seqid, "-db", blastdb] # todo: add 'blastdbcmd' to dependencies and to configs
	p_blastdbcmd = subprocess.run(cmd_blastdbcmd, stdout = subprocess.PIPE, stderr = subprocess.PIPE, text = True)
	try:
		p_blastdbcmd.check_returncode()
	except Exception:
		sys.stderr.write("\nAn error occured when trying to extract seqids {} from blast database {}\n".format(seqid, blastdb))
		sys.stderr.write("{}\n".format(p_blastdbcmd.stderr))
		raise RuntimeError
	fasta_io = StringIO(p_blastdbcmd.stdout)
	return(list(SeqIO.parse(fasta_io, "fasta")))
	

print(get_contig_from_blastdb(["GCF_000812565.1_NZ_JTAM01000055.1","GCF_001266905.1_NZ_LHYY01000050.1"], "./deltestdb"))
