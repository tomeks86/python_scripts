#!/usr/bin/python

import sys

fin=open(sys.argv[1]).readlines()

minn=0
mini="TCL: Minimizing"
for i in range(len(fin)):
	if fin[i].startswith(mini):
		minn=int(fin[i].split()[3])
		break

#print minn

ener="ENERGY:"
for i in range(len(fin)-1,0,-1):
	if fin[i].startswith(ener):
		break

print int(fin[i].split()[1])-minn
