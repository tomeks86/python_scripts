#!/usr/bin/python

import numpy as np
import sys
import openbabel as ob
import copy
import random as rand
from random import random as rdm
from math import floor,ceil

table=ob.OBElementTable()

defns=["ATOM/HETATM","atom serial number","atom name","alternate location indicator","residue name","chain identifier","residue sequence number","code for insertion of residues","orthogonal coordinates for X (in Angstroms)","orthogonal coordinates for Y (in Angstroms)","orthogonal coordinates for Z (in Angstroms)","occupancy","temperature factor","site identifier","element symbol","charge on the atom"]
formats=["{:6s}","{:5d}","{:>5s}","{:1s}","{:4s}","{:1s}","{:4d}","{:1s}","{:11.3f}","{:8.3f}","{:8.3f}","{:6.2f}","{:6.2f}","{:>10s}","{:>2s}","{:2s}"]
formatx=["{:6s}","{:5x}","{:>5s}","{:1s}","{:4s}","{:1s}","{:4d}","{:1s}","{:11.3f}","{:8.3f}","{:8.3f}","{:6.2f}","{:6.2f}","{:>10s}","{:>2s}","{:2s}"]
# 11.3f mozna zamienic na 10.3f, w razie potrzeby...

def rot_x(phi):
	sph=np.sin(phi)
	cph=np.cos(phi)
	X=np.identity(3)
	X[1,1]=cph
	X[2,2]=cph
	X[2,1]=sph
	X[1,2]=-sph
	return X

def rot_y(phi):
	sph=np.sin(phi)
	cph=np.cos(phi)
	X=np.identity(3)
	X[0,0]=cph
	X[2,2]=cph
	X[0,2]=sph
	X[2,0]=-sph
	return X

def rot_z(phi):
	sph=np.sin(phi)
	cph=np.cos(phi)
	X=np.identity(3)
	X[0,0]=cph
	X[1,1]=cph
	X[1,0]=sph
	X[0,1]=-sph
	return X	

def calcI(mol):
	cm=CM(mol)
	I=np.zeros((3,3))
	N=len(mol)
	for i in range(N):
		at=mol[i]
		m=get_mass(at)
		xyz=get_coord(at)-cm
		I[0,0]+=m*(xyz[1]**2+xyz[2]**2)
		I[1,1]+=m*(xyz[0]**2+xyz[2]**2)
		I[2,2]+=m*(xyz[0]**2+xyz[1]**2)
		I[0,1]-=m*xyz[0]*xyz[1]
		I[0,2]-=m*xyz[0]*xyz[2]
		I[1,2]-=m*xyz[1]*xyz[2]
	I[1,0]=I[0,1]
	I[2,0]=I[0,2]
	I[2,1]=I[1,2]
	return I

def get_mass(at):
	return table.GetMass(table.GetAtomicNum(at["element symbol"]))

def get_coord(at):
	xyz=np.zeros(3)
	xyz[0]=at["orthogonal coordinates for X (in Angstroms)"]
	xyz[1]=at["orthogonal coordinates for Y (in Angstroms)"]
	xyz[2]=at["orthogonal coordinates for Z (in Angstroms)"]
	return xyz

def write_coord(xyz,at):
	at["orthogonal coordinates for X (in Angstroms)"]=xyz[0]
	at["orthogonal coordinates for Y (in Angstroms)"]=xyz[1]
	at["orthogonal coordinates for Z (in Angstroms)"]=xyz[2]

def get_sn(at):
	#if at["residue sequence number"] is not None:
	return at["residue sequence number"]
	#else:
	#	print at;exit()

def write_sn(sn,at):
	at["residue sequence number"]=sn

def get_asn(at):
	return at["atom serial number"]

def write_asn(asn,at):
	at["atom serial number"]=asn

def set_id(mol,ide):
	for i in range(len(mol)):
		mol[i]["site identifier"]=ide

def CM(mol):
	cm=np.zeros(3)
	#M=mol.GetExactMass()
	M=0.
	N=len(mol)
	for i in range(N):
		m=get_mass(mol[i])
		M+=m
		cm+=m*get_coord(mol[i])
	return cm/M

def translate_mol(mol,V):
	N=len(mol)
	for i in range(N):
		at=mol[i]
		at["orthogonal coordinates for X (in Angstroms)"]+=V[0]
		at["orthogonal coordinates for Y (in Angstroms)"]+=V[1]
		at["orthogonal coordinates for Z (in Angstroms)"]+=V[2]
#	return mol

def rotate_mol(mol,X):
	N=len(mol)
	for i in range(N):
		at=mol[i]
		vec=get_coord(at)
		vec=np.dot(X.T,vec)
		write_coord(vec,at)
#	return mol

def copy_mol(mol):
	return copy.deepcopy(mol)

def add_molecule(mol,klast,X=np.identity(3),tr=np.zeros(3)):	#warstwy!
	try:
		rsn=get_sn(klast[-1])
	except IndexError:
		rsn=0
	mol2=copy_mol(mol)
	rotate_mol(mol2,X)
	translate_mol(mol2,tr)
	for i in range(len(mol2)):
		write_sn(get_sn(mol2[i])+rsn,mol2[i])
		klast.append(mol2[i])
#	return klast

def combine_molecule(mol,klast,X=np.identity(3),tr=np.zeros(3)):	#&woda!
	try:
		asn=get_asn(klast[-1])
	except IndexError:
		asn=0
	mol2=copy_mol(mol)
	rotate_mol(mol2,X)
	translate_mol(mol2,tr)
	for i in range(len(mol2)):
		write_asn(get_asn(mol2[i])+asn,mol2[i])
		klast.append(mol2[i])

def read_xyz(lines,resn):
	molecule=[]
	N=int(lines[0].split()[0])
	for i in range(N):
		line=lines[i+2].split()
		at={}
		at["ATOM/HETATM"]="ATOM"
		at["atom serial number"]=i+1
		at["atom name"]=line[0]+str(i+1)
		at["alternate location indicator"]="L"
		at["residue name"]=resn
		at["chain identifier"]="X"
		at["residue sequence number"]=1
		at["code for insertion of residues"]=""
		at["orthogonal coordinates for X (in Angstroms)"]=float(line[1])
		at["orthogonal coordinates for Y (in Angstroms)"]=float(line[2])
		at["orthogonal coordinates for Z (in Angstroms)"]=float(line[3])
		at["occupancy"]=1.
		at["temperature factor"]=0.
		at["site identifier"]="lig"
		at["element symbol"]=line[0]
		at["charge on the atom"]=""
		molecule.append(at)
	return molecule

def read_pdb(lines):
	molecule=[]
	hx=False
	star=False
	for line in lines:
		if "ATOM" not in line:
			continue
		tmp={}
		tmp[defns[0]]=line[:6].strip()
		if not hx:
			try:
				tmp[defns[1]]=int(line[6:11])
			except ValueError:
				hx=True
#				print line[6:12]
				try:
					tmp[defns[1]]=int(line[6:11],16)
				except ValueError:
					star=True
					tmp[defns[1]]=0
		elif not star:
			tmp[defns[1]]=tmp[defns[1]]=int(line[6:11],16)
		else:
			tmp[defns[1]]=0
		tmp[defns[2]]=line[12:16].strip()
		tmp[defns[3]]=line[16:17].strip()
		tmp[defns[4]]=line[17:21].strip()
		tmp[defns[5]]=line[21:22].strip()
		tmp[defns[6]]=int(line[22:26])
		tmp[defns[7]]=line[26:27].strip()
		tmp[defns[8]]=float(line[30:38])
		tmp[defns[9]]=float(line[38:46])
		tmp[defns[10]]=float(line[46:54])
		tmp[defns[11]]=float(line[54:60])
		tmp[defns[12]]=float(line[60:66])
		rest=line[66:].split()
		if len(rest)==1:
			tmp[defns[13]]=rest[0]	#to trzeba czasem zamienic
			tmp[defns[14]]=""		#z tym
			tmp[defns[15]]=""
		elif len(rest)==2:
			tmp[defns[13]]=rest[0]
			tmp[defns[14]]=rest[1]
			tmp[defns[15]]=""
		elif len(rest)==3:
			tmp[defns[13]]=rest[0]
			tmp[defns[14]]=rest[1]
			tmp[defns[15]]=rest[2]
		molecule.append(tmp)
		#print tmp;exit()
	return molecule

def write_pdb(molecule,fout):
	lout=""
	for i in range(len(molecule)):
		for k in range(15):
			#print lout
			if molecule[i]["atom serial number"]>99999:
				lout+=formatx[k].format(molecule[i][defns[k]])
			else:
				#print molecule[i]["atom serial number"]
				lout+=formats[k].format(molecule[i][defns[k]])
		lout+="\n"
	#print lout
	fout=open(fout,'w')
	fout.write(lout)
	fout.close()

#def extr_woda(mol):
#	woda=[]
#	bot=[]
#	top=[]
#	for i in range(len(mol)):
#		if mol[i]["site identifier"].startswith("WT"):
#			woda.append(mol[i])
#		elif mol[i]["site identifier"].startswith("TOP"):
#			top.append(mol[i])
#		elif mol[i]["site identifier"].startswith("BOT"):
#			bot.append(mol[i])
#	return woda,top,bot

def ren_asn(mol):	#przenumerowanie asn w przypadku gwiazdek!
	for i in range(len(mol)):
		mol[i]["atom serial number"]=i+1

def extract(mol,et):
	out=[]
	for i in range(len(mol)):
		if mol[i]["site identifier"].startswith(et):
			out.append(mol[i])
	ren_asn(out)
	return out

def num_woda(mol):	#ponumeruj wode
	N=len(mol)/3
	n_seg=int(ceil(N/9999.))
	licz=1
	for n in range(n_seg-1):
		sg_name="WT"+str(n+1)
		for i in range(1,10000):	#residue sequence number"
			for j in range(3):	#water molecules
				mol[licz-1]["atom serial number"]=licz
				mol[licz-1]["site identifier"]=sg_name
				mol[licz-1]["residue sequence number"]=i
				licz+=1
	#jeszcze n_seg...
	akt=licz
	sg_name="WT"+str(n_seg)
#	print akt, N*3;exit()
	k=1
	for i in range((akt-1)/3,N):
		for j in range(3):
			mol[licz-1]["atom serial number"]=licz
			mol[licz-1]["site identifier"]=sg_name
			mol[licz-1]["residue sequence number"]=k
			licz+=1
		k+=1

def sort_sid(mol,et):	#sortuje site identifier, np. WT1 WT2, TOP
	mol1=[]
	mol2=[]
	for i in range(len(mol)):
		if mol[i]["site identifier"]==et:
			mol1.append(mol[i])
		else:
			mol2.append(mol[i])
	add_molecule(mol2,mol1)
	return mol1

def sort_resn(mol,et):	#sortuje resname
	mol1=[]
	mol2=[]
	for i in range(len(mol)):
		if mol[i]["residue name"]==et:
			mol1.append(mol[i])
		else:
			mol2.append(mol[i])
	add_molecule(mol2,mol1)
	ren_asn(mol1)
	#write_pdb(mol1,"test.pdb");exit()
	return mol1

def ren_rsn(mol,N):	#przenumerowanie, np. mocznik co 8 atomow
	I=1
	akt=0
	k=int(ceil(len(mol)/float(N)))
	for i in range(k):
		for j in range(N):
			mol[akt]["residue sequence number"]=I
			akt+=1
		I+=1

def get_resnames(mol):
	nm=[]
	for i in range(len(mol)):
		resnm=mol[i]["residue name"]
		if not resnm in nm:
			nm.append(resnm)
	nm.sort(reverse=True)
	return nm

def sort_resname_chainid(mol):	#sortuje resname w ramach chain_id
	chids=[]
	for i in range(len(mol)):
		if not mol[i]["chain identifier"] in chids:
			chids.append(mol[i]["chain identifier"])
	#chids.sort()
	mols=[[] for i in range(len(chids))]
	out=[]
	for i in range(len(mol)):
		idakt=mol[i]["chain identifier"]
		mols[chids.index(idakt)].append(mol[i])
	for i in range(len(mols)):
		for rsnm in get_resnames(mols[i]):
			mols[i]=sort_resn(mols[i],rsnm)
			ren_resn(mols[i])
		for j in range(len(mols[i])):
			out.append(mols[i][j])
	ren_asn(out)
	return out

#def ren_resn(mol):	#przenumerowanie reszt w ramach chain -> do mols[i] sort_resname_chainid (po posortowaniu!!!)
#	out=[]
#	rnms=get_resnames(mol)
#	mols=[[] for i in range(len(rnms))]
#	for i in range(len(mol)):
#		mols[rnms.index(mol[i]["residue name"])].append(mol[i])
#	for i in range(len(rnms)):
#		k=0
#		akt=mols[i][0]["residue sequence number"]
#		for j in range(len(mols[i])):
#			if mols[i][j]["residue sequence number"] != akt:
#				akt=mols[i][j]["residue sequence number"]
#				k+=1
#			mols[i][j]["residue sequence number"]=k%9999+1
#			out.append(mols[i][j])
#	return out

def ren_resn(mol):	#przenumerowanie reszt w ramach chain -> do mols[i] sort_resname_chainid (po posortowaniu!!!)
	k=0
	akt=mol[0]["residue sequence number"]
	for i in range(len(mol)):
		if mol[i]["residue sequence number"] != akt:
			akt=mol[i]["residue sequence number"]
			k+=1
		mol[i]["residue sequence number"]=k%9999+1
