#!/usr/bin/python

from pdb_tools import *

mol1=read_pdb(open(sys.argv[1]).readlines())
mol2=read_pdb(open(sys.argv[2]).readlines())

suma=[]
add_molecule(mol1,suma)
add_molecule(mol2,suma)

try:
	write_pdb(suma,sys.argv[3])
except IndexError:
	write_pdb(suma,"out.pdb")
