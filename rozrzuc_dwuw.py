#!/usr/bin/python

from pdb_tools import *
from numpy.linalg import eigh, det

lpar=open(sys.argv[1]).readlines()
fgeom=open(lpar[0].split()[0]).readlines()
X,Y=[float(u) for u in lpar[1].split()[0:2]]
Na,Nb,r=[int(u) for u in lpar[2].split()[0:3]]
Z=float(lpar[3].split()[0])

mol=read_pdb(fgeom)

translate_mol(mol,-CM(mol))

I=calcI(mol)
D=eigh(I)[1]
sgm=np.array([0,0,1,0,-1,0,1,0,0]).reshape((3,3))
D=np.dot(D,sgm)
if det(D)<0:	#jesli odwrocil chiralnosc
	D=np.dot(-np.identity(3),D)
rotate_mol(mol,D)

klast=[]
NN=Na*Nb+r
dNa=0
dNb=0
if r!=0:
	if X>=Y:
		dNa=1
	else:
		dNb=1
dX=X/(Na+dNa)
dY=Y/(Nb+dNb)

warstwaA=[]
for i in range(Na):
	for j in range(Nb):
		x=i*dX+.5*dX
		y=j*dY+.5*dY
		XYZ=np.dot(rot_z(2.*np.pi*rdm()),rot_x(np.pi))
		#XYZ=rot_z(2.*np.pi*rdm())
		#XYZ=np.dot(rot_x(2.*np.pi*rdm()),np.dot(rot_y(2.*np.pi*rdm()),rot_z(2.*np.pi*rdm())))
		add_molecule(mol,warstwaA,XYZ,np.array([x+(rdm()-.5)*dX,y+(rdm()-.5)*dY,Z+(rdm()-.5)+2]))
if dNa:
	x+=dX
	dY=Y/r
	for i in range(r):
		y=i*dY+.5*dY
		phi=2.*np.pi*rand.random()
		XYZ=np.dot(rot_z(2.*np.pi*rdm()),rot_x(np.pi))
		#XYZ=rot_z(2.*np.pi*rdm())
		#XYZ=np.dot(rot_x(2.*np.pi*rdm()),np.dot(rot_y(2.*np.pi*rdm()),rot_z(2.*np.pi*rdm())))
		add_molecule(mol,warstwaA,XYZ,np.array([x+(rdm()-.5)*dX,y+(rdm()-.5)*dY,Z+(rdm()-.5)+2]))
elif dNb:
	y+=dY
	dX=X/r
	for i in range(r):
		x=i*dX+.5*dX
		phi=2.*np.pi*rand.random()
		XYZ=np.dot(rot_z(2.*np.pi*rdm()),rot_x(np.pi))
		#XYZ=rot_z(2.*np.pi*rdm())
		#XYZ=np.dot(rot_x(2.*np.pi*rdm()),np.dot(rot_y(2.*np.pi*rdm()),rot_z(2.*np.pi*rdm())))
		add_molecule(mol,warstwaA,XYZ,np.array([x+(rdm()-.5)*dX,y+(rdm()-.5)*dY,Z+(rdm()-.5)+2]))
set_id(warstwaA,"TOP")

#druga warstwa
warstwaB=[]
for i in range(Na):
	for j in range(Nb):
		x=i*dX+.5*dX
		y=j*dY+.5*dY
		XYZ=rot_z(2.*np.pi*rdm())
		#XYZ=np.dot(rot_x(2.*np.pi*rdm()),np.dot(rot_y(2.*np.pi*rdm()),rot_z(2.*np.pi*rdm())))
		add_molecule(mol,warstwaB,XYZ,np.array([x+(rdm()-.5)*dX,y+(rdm()-.5)*dY,-Z+(rdm()-.5)+2]))
if dNa:
	x+=dX
	dY=Y/r
	for i in range(r):
		y=i*dY+.5*dY
		phi=2.*np.pi*rand.random()
		XYZ=rot_z(2.*np.pi*rdm())
		#XYZ=np.dot(rot_x(2.*np.pi*rdm()),np.dot(rot_y(2.*np.pi*rdm()),rot_z(2.*np.pi*rdm())))
		add_molecule(mol,warstwaB,XYZ,np.array([x+(rdm()-.5)*dX,y+(rdm()-.5)*dY,-Z+(rdm()-.5)+2]))
elif dNb:
	y+=dY
	dX=X/r
	for i in range(r):
		x=i*dX+.5*dX
		phi=2.*np.pi*rand.random()
		XYZ=rot_z(2.*np.pi*rdm())
		#XYZ=np.dot(rot_x(2.*np.pi*rdm()),np.dot(rot_y(2.*np.pi*rdm()),rot_z(2.*np.pi*rdm())))
		add_molecule(mol,warstwaB,XYZ,np.array([x+(rdm()-.5)*dX,y+(rdm()-.5)*dY,-Z+(rdm()-.5)+2]))
set_id(warstwaB,"BOT")

#add_molecule(warstwaB,warstwaA)
#print warstwaA[0]["site identifier"];exit()
#translate_mol(warstwaA,np.array([0,0,100]))
#translate_mol(warstwaB,np.array([0,0,100]))

write_pdb(warstwaA,"top.pdb")
write_pdb(warstwaB,"bot.pdb")
