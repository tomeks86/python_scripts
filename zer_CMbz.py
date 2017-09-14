#!/usr/bin/python
from pdb_tools import *

mol=read_pdb(open(sys.argv[1]).readlines())

cm=CM(mol)
cm[2]=0.	#nie przesuwam nic na z!!!

translate_mol(mol,-cm)

write_pdb(mol,sys.argv[1])
#write_pdb(mol,"test.pdb")
