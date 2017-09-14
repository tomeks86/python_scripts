#!/usr/bin/python
from pdb_tools import *

mol=read_pdb(open(sys.argv[1]).readlines())

translate_mol(mol,np.array([0,0,float(sys.argv[2])]))

write_pdb(mol,sys.argv[1])
#write_pdb(mol,"test.pdb")
