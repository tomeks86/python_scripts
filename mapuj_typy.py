#!/usr/bin/python

import openbabel as ob
import sys
import numpy as np

table=ob.OBElementTable()

mol=ob.OBMol()
conv=ob.OBConversion()
conv.SetInAndOutFormats('xyz','pdb')
conv.ReadFile(mol,sys.argv[1])

#mapa=open(sys.argv[2])

N=mol.NumAtoms()

#for i in range(N):
#	at=mol.GetAtom(i+1)
#	print i+1, at.GetType()

#mol.write("smi")

#conv.WriteFile(mol,"test.pdb")

matC=np.zeros((N,N))
for i in range(mol.NumBonds()):
	bnd=mol.GetBond(i)
	a=bnd.GetBeginAtomIdx()
	b=bnd.GetEndAtomIdx()
	matC[a-1,b-1]=1
	matC[b-1,a-1]=1

#for i in range(N):
#	for j in range(N):
#		if matC[i,j] != 0:
#			print "%2i" % matC[i,j],
#		else:
#			print "  ",
#	print

for i in range(N):
	props=[]
	at=mol.GetAtom(i+1)
	props.append(table.GetSymbol(at.GetAtomicNum()))
	props.append(at.GetAtomicNum())
	#props.append(at.GetFormalCharge())
	props.append(at.GetValence())
	props.append(at.GetHvyValence())
	props.append(at.GetType())
	props.append(at.GetPartialCharge())
	#props.append(at.isAromatic())
	props.append(at.IsAxial())
	for p in props:
		print p,
	print
