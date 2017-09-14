#!/usr/bin/python

from pdb_tools import *
from numpy.linalg import eigh, det

lpar=open(sys.argv[1]).readlines()
fgeom1=open(lpar[0].split()[0]).readlines()
fgeom2=open(lpar[1].split()[0]).readlines()
X,Y=[float(u) for u in lpar[2].split()[0:2]]
Na,Nb,r=[int(u) for u in lpar[3].split()[0:3]]
Z=float(lpar[4].split()[0])
na,nb=[int(u) for u in lpar[5].split()[:2]]

if na+nb != Na*Nb+r:
	print "nie da rady tak rozrzucic, na+nb: %d rozne od Na*Nb+r: %d" % (na+nb,Na*Nb+r)
	exit()

#sel={'a':na, 'b':nb}
#def pick(x):
#	a=x['a']
#	b=x['b']
#	tot=a+b
#	if tot==0:
#		print "wyczerpano ilosc czasteczek!"
#		exit()
#	ind=floor(rand.random()*tot)
#	if ind<a:
#		x['a']-=1
#		return 0
#	else:
#		x['b']-=1
#		return 1

def pick(x):
	if x>=na:
		return 1
	else:
		return 0

mol1=read_pdb(fgeom1)
mol2=read_pdb(fgeom2)

translate_mol(mol1,-CM(mol1))
translate_mol(mol2,-CM(mol2))

#for mol in [mol1,mol2]:
#	I=calcI(mol)
#	D=eigh(I)[1]
#	sgm=np.array([0,0,1,0,-1,0,1,0,0]).reshape((3,3))
#	D=np.dot(D,sgm)
#	if det(D)<0:	#jesli odwrocil chiralnosc
#		D=np.dot(-np.identity(3),D)
#	rotate_mol(mol,D)

mol=[mol1,mol2]

#orig_mol=copy_mol(mol)
#klast=copy_mol(orig_mol)
#add_molecule(orig_mol,klast,np.identity(3),np.array([0,10,15]))
#add_molecule(orig_mol,klast,np.identity(3),np.array([15,10,15]))

#klast=[]
#add_molecule(mol,klast)
#write_pdb(klast,"test.pdb");exit()

dX=X/Na
dY=X/Nb

NN=Na*Nb+r
dNa=0
dNb=0
if r!=0:
	if Nb>Na:
		dNa=1
	else:
		dNb=1
dX=X/(Na+dNa)
dY=Y/(Nb+dNb)

idx=range(na+nb)

rand.shuffle(idx)
k=0
warstwaA=[]
for i in range(Na):
	for j in range(Nb):
		x=i*dX+.5*dX
		y=j*dY+.5*dY
		XYZ=np.dot(rot_z(2.*np.pi*rdm()),rot_x(np.pi))
		#XYZ=rot_z(2.*np.pi*rdm())
		#XYZ=np.dot(rot_x(2.*np.pi*rdm()),np.dot(rot_y(2.*np.pi*rdm()),rot_z(2.*np.pi*rdm())))
		add_molecule(mol[pick(idx[k])],warstwaA,XYZ,np.array([x+(rdm()-.5)*dX,y+(rdm()-.5)*dY,Z+(rdm()-.5)+2]))
		k+=1
if dNa:
	x+=dX
	dY=Y/r
	for i in range(r):
		y=i*dY+.5*dY
		phi=2.*np.pi*rand.random()
		XYZ=np.dot(rot_z(2.*np.pi*rdm()),rot_x(np.pi))
		#XYZ=rot_z(2.*np.pi*rdm())
		#XYZ=np.dot(rot_x(2.*np.pi*rdm()),np.dot(rot_y(2.*np.pi*rdm()),rot_z(2.*np.pi*rdm())))
		add_molecule(mol[pick(idx[k])],warstwaA,XYZ,np.array([x+(rdm()-.5)*dX,y+(rdm()-.5)*dY,Z+(rdm()-.5)+2]))
		k+=1
elif dNb:
	y+=dY
	dX=X/r
	for i in range(r):
		x=i*dX+.5*dX
		phi=2.*np.pi*rand.random()
		XYZ=np.dot(rot_z(2.*np.pi*rdm()),rot_x(np.pi))
		#XYZ=rot_z(2.*np.pi*rdm())
		#XYZ=np.dot(rot_x(2.*np.pi*rdm()),np.dot(rot_y(2.*np.pi*rdm()),rot_z(2.*np.pi*rdm())))
		add_molecule(mol[pick(idx[k])],warstwaA,XYZ,np.array([x+(rdm()-.5)*dX,y+(rdm()-.5)*dY,Z+(rdm()-.5)+2]))
		k+=1
set_id(warstwaA,"TOP")

#druga warstwa

dX=X/(Na+dNa)
dY=Y/(Nb+dNb)

rand.shuffle(idx)
k=0
warstwaB=[]
for i in range(Na):
	for j in range(Nb):
		x=i*dX+.5*dX
		y=j*dY+.5*dY
		XYZ=rot_z(2.*np.pi*rdm())
		#XYZ=np.dot(rot_x(2.*np.pi*rdm()),np.dot(rot_y(2.*np.pi*rdm()),rot_z(2.*np.pi*rdm())))
		add_molecule(mol[pick(idx[k])],warstwaB,XYZ,np.array([x+(rdm()-.5)*dX,y+(rdm()-.5)*dY,-Z+(rdm()-.5)+2]))
		k+=1
if dNa:
	x+=dX
	dY=Y/r
	for i in range(r):
		y=i*dY+.5*dY
		phi=2.*np.pi*rand.random()
		XYZ=rot_z(2.*np.pi*rdm())
		#XYZ=np.dot(rot_x(2.*np.pi*rdm()),np.dot(rot_y(2.*np.pi*rdm()),rot_z(2.*np.pi*rdm())))
		add_molecule(mol[pick(idx[k])],warstwaB,XYZ,np.array([x+(rdm()-.5)*dX,y+(rdm()-.5)*dY,-Z+(rdm()-.5)+2]))
		k+=1
elif dNb:
	y+=dY
	dX=X/r
	for i in range(r):
		x=i*dX+.5*dX
		phi=2.*np.pi*rand.random()
		XYZ=rot_z(2.*np.pi*rdm())
		#XYZ=np.dot(rot_x(2.*np.pi*rdm()),np.dot(rot_y(2.*np.pi*rdm()),rot_z(2.*np.pi*rdm())))
		add_molecule(mol[pick(idx[k])],warstwaB,XYZ,np.array([x+(rdm()-.5)*dX,y+(rdm()-.5)*dY,-Z+(rdm()-.5)+2]))
		k+=1
set_id(warstwaB,"BOT")

#add_molecule(warstwaB,warstwaA)
#print warstwaA[0]["site identifier"];exit()

#translate_mol(warstwaA,np.array([0,0,100]))
#translate_mol(warstwaB,np.array([0,0,100]))

ren_asn(warstwaA)
ren_asn(warstwaB)

write_pdb(warstwaA,"top.pdb")
write_pdb(warstwaB,"bot.pdb")

warstwa=[]
add_molecule(warstwaA,warstwa)
add_molecule(warstwaB,warstwa)
write_pdb(warstwa,"warstwa.pdb")
