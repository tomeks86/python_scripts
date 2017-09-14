#!/usr/bin/python

from pdb_tools import *

mol=read_pdb(open(sys.argv[1]).readlines())

print mol[1]
