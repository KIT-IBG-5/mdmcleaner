#!/usr/bin/env python
import sys
import os
import argparse

from _version import __version__

#todo: consider using "fromfile_prefix_chars" from ArgumentParser to optionally read arguments from file
myparser = argparse.ArgumentParser(prog=os.path.basename(sys.argv[0]), description= "identifies and removes potential contamination in draft genomes and metagenomic bins, based on a hierarchically ranked contig classification pipeline")
myparser.add_argument("-v", "--version", action = "version", version='%(prog)s {version}'.format(version=__version__))
