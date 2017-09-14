#!/usr/bin/python

from pdb_tools import *

mol=read_pdb(open(sys.argv[1]).readlines())
N=int(sys.argv[2])

ren_rsn(mol,N)

write_pdb(mol,sys.argv[1])
