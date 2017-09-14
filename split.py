#!/usr/bin/python

from pdb_tools import *

mol=read_pdb(open(sys.argv[1]).readlines())

ets=[]

for i in range(len(mol)):
	ets.append(mol[i]["site identifier"])
#print set(ets);exit()
for et in set(ets):
	moll=extract(mol,et)
	if et is not None:
		write_pdb(moll,et+".pdb")
	else:
		write_pdb(moll,"blank.pdb")
