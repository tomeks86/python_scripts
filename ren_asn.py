#!/usr/bin/python

from pdb_tools import *

mol=read_pdb(open(sys.argv[1]).readlines())

ren_asn(mol)

write_pdb(mol,sys.argv[1])
