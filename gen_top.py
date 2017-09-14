#!/usr/bin/python

import openbabel as ob
import sys

table=ob.OBElementTable()

mol=ob.OBMol()
conv=ob.OBConversion()
conv.SetInAndOutFormats('xyz','xyz')
conv.ReadFile(mol,sys.argv[1])

conv={"Nam":"NH1","C3":"CT3","HC":"HHH","C2":"C","O2":"O"}

#n=mol.NumAtoms()

#k=0
#l=0
for i in range(1,mol.NumAtoms()+1):
	at=mol.GetAtom(i)
	name=table.GetSymbol(at.GetAtomicNum())
	typ=at.GetType()
	val=at.GetImplicitValence()
	if name=='H' and mol.GetAtom(i-1).GetType()=="Nam":
		typ='H'
	elif name=='H' and k!=0:
		typ="HA"+str(k)
	elif name=='O' and mol.GetAtom(i-1).GetType()=="O2":
		typ="OS"
	elif name=="C" and val==4:
		n=4-at.GetHvyValence()
		if n:
			typ="CT"+str(n)
		else:
			typ="CT"
	if not at.IsHydrogen():
		if name=="C" and val==4:
			k=4-at.GetHvyValence()
		else:
			k=0
	try:
		typ=conv[typ]
	except KeyError:
		pass
	charge=at.GetPartialCharge()
	if charge>=0:
		print "ATOM %-4s %-4s    %.2f" % (name+str(i),typ,charge)
	else:
		print "ATOM %-4s %-4s   %.2f" % (name+str(i),typ,charge)
	
#print "\nBOND ",
j=0
for i in range(mol.NumBonds()):
	if (j%8 == 0):
		print "\nBOND ",
	bond=mol.GetBond(i)
	at1=bond.GetBeginAtom()
	at1name=table.GetSymbol(at1.GetAtomicNum())
	at1idx=at1.GetIdx()
	at2=bond.GetEndAtom()
	at2name=table.GetSymbol(at2.GetAtomicNum())
	at2idx=at2.GetIdx()
	print "%s%i %s%i" % (at1name,at1idx,at2name,at2idx),
	j+=1
	
print
