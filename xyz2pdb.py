#!/usr/bin/python
from pdb_tools import *
from os.path import exists

if len(sys.argv)<3:
	print "uzycie: xyz2pdb.py plik.xyz resname"
	exit()

finp=sys.argv[1]
resn=sys.argv[2]
if not finp.endswith(".xyz"):
	print "podaj plik xyz!"
	exit()

fout=finp.strip(".xyz")
if exists(fout+".pdb"):
	fout+="2.pdb"
else:
	fout+=".pdb"

mol=read_xyz(open(finp).readlines(),resn)

#print mol[0];exit()
write_pdb(mol,fout)
#write_pdb(mol,"test.pdb")
