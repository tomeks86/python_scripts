#!/usr/bin/python
# sortuje site 

from pdb_tools import *

mol=read_pdb(open(sys.argv[1]).readlines())

mol=sort_resname_chainid(mol)
#mol=sort_resn(mol,"DPPC")

#if len(sys.argv)>2:
#	out=sys.argv[2]
#else:
#	out="out.pdb"

write_pdb(mol,sys.argv[1])
#write_pdb(mol,"test.pdb")
