#!/usr/bin/python

from pdb_tools import *

mol=read_pdb(open(sys.argv[1]).readlines())

N=int(sys.argv[2]) # docelowa ilosc czasteczek wody

top=extract(mol,"TOP")
bot=extract(mol,"BOT")
woda=extract(mol,"WT")

Nw=len(woda)/3-N
if Nw<0:
	print "za malo czasteczek...!"
	exit()

idx=range(Nw+N)
rand.shuffle(idx)
idx=idx[:N]

woda_new=[]
for i in idx:
	for j in range(3):
		woda_new.append(woda[3*i+j])

num_woda(woda_new)

write_pdb(top,"topp.pdb")
write_pdb(bot,"bott.pdb")
write_pdb(woda_new,"watt.pdb")

