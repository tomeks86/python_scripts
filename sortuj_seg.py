#!/usr/bin/python
# sortuje site id (np. TOP BOT)

from pdb_tools import *

mol=read_pdb(open(sys.argv[1]).readlines())

mol=sort_sid(mol,sys.argv[2])

if len(sys.argv)>3:
	out=sys.argv[3]
else:
	out="out.pdb"

write_pdb(mol,out)
