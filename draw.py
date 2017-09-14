#!/usr/bin/python

import numpy as np
from numpy.linalg import norm, inv, det, eigh, eig
import vtk
import openbabel as ob
from optparse import OptionParser
import Tkinter as tk

usage = "usage: %prog [options] "
parser = OptionParser(usage=usage)

parser.add_option("-i", "--inputfile", dest="inputfile",
                  help="xyz file",
                  action="store", type="string",default="none")
parser.add_option("-a", "--arrows", dest="arrows",
                  help="arrow file",
                  action="store", type="string", default="none")
parser.add_option("-b", "--arrows2", dest="arrows2",
                  help="arrow file",
                  action="store", type="string", default="none")
parser.add_option("-n", "--nobonds", dest="nobonds",
                  help="without bonds",
                  action="store_true", default=False)
parser.add_option("-c", "--cell", dest="cell",
                  help="draws unit cell",
                  action="store_true", default=False)
parser.add_option("-p", "--prime", dest="prime",
                  help="if cell parameters in line 2 of xyz file, but start from position 2",
                  action="store_true", default=False)
parser.add_option("-t", "--schema", dest="schema",
                  help="to change color schema (to white background)",
                  action="store_true", default=False)
parser.add_option("-s", "--scales", dest="scales",
                   help="sets scales for beta vectorial, dipole moment, electric field",
                   action="store", type="string", default=".5,.5,.5")
parser.add_option("-x", "--axes", dest="axes",
                   help="draws local coordinate axes sytem at selected position",
                   action="store", type="string", default="none")
parser.add_option("-r", "--polar", dest="polar",
                   help="draws elipsoids of atomic polarizabilities",
                   action="store", type="string", default="none")
parser.add_option("-v", "--vector", dest="vector",
                   help="saves vector image svg",
                   action="store", type="string", default="none")

(options, args) = parser.parse_args()

if options.inputfile=="none":
	parser.error("No input file given")

schema=[(0,0,0),(1,1,1)]	#background, lines
if options.schema:
	schema=[(1,1,1),(0,0,0)]

scales=[.5,.5,.5]
no={"b":0,"m":1,"m1":1,"m2":1,"m3":1,"d":1,"d1":1,"d2":1,"d3":1,"e":2,"e1":2,"e2":2,"f":2,"f1":2,"f2":2}
if options.scales <> ".5,.5,.5":
	sk=options.scales.split(',')
	if len(sk)==1:
		scales[0]=float(sk[0])
	elif len(sk)==2:
		scales[0]=float(sk[0])
		scales[1]=float(sk[1])
	elif len(sk)==3:
		scales[0]=float(sk[0])
		scales[1]=float(sk[1])
		scales[2]=float(sk[2])
	else:
		parser.error("incorrect scales, have in mind the sequence: beta vectorial, dipole moment, electric field")

table=ob.OBElementTable()
mol=ob.OBMol()
conv=ob.OBConversion()
conv.SetInAndOutFormats('xyz','xyz')
if not conv.ReadFile(mol,options.inputfile):
	raise Exception('Unable to read file')

def drawAtom(atom,darken=False):
	scale=.1
	ball=vtk.vtkSphereSource()
	ball.SetPhiResolution(50)
	ball.SetThetaResolution(50)
	x=atom.GetX()
	y=atom.GetY()
	z=atom.GetZ()
	ball.SetCenter(x,y,z)
	ball.SetRadius(table.GetVdwRad(atom.GetAtomicNum())*scale)
	ballMapper=vtk.vtkPolyDataMapper()
	ballMapper.SetInput(ball.GetOutput())
	ballActor=vtk.vtkActor()
	ballActor.SetMapper(ballMapper)
	color=list(table.GetRGB(atom.GetAtomicNum()))
	if darken:
		for i in range(3):
			color[i]*=.5
	ballActor.GetProperty().SetColor(color)
	ballActor.GetProperty().SetOpacity(20.)
	return ballActor

def simeq(a,b,eps=1e-7):
	if abs(a-b)<eps:
		return True
	else:
		return False

def transformer(obj,translation=np.zeros(3),scaling=np.ones(3),angle=0.,rotationaxis=np.zeros(3)):
	trans=vtk.vtkTransform()
	trans.PostMultiply()
	trans.Scale(scaling)
	trans.RotateWXYZ(angle,rotationaxis[0],rotationaxis[1],rotationaxis[2])
	trans.Translate(translation)
	tf = vtk.vtkTransformPolyDataFilter()
	tf.SetInput(obj.GetOutput())
	tf.SetTransform(trans)
	tf.Update()
	return tf.GetOutput()

def genEllipse(xyz,eps,scale=.07):
	ball=vtk.vtkSphereSource()
	ball.SetPhiResolution(20)
	ball.SetThetaResolution(20)
	E,X=eigh(eps)
	X=X.T
	if det(X)<0: #abs((np.trace(X)-1.)/2.)>1:
		X*=-1
	phi=np.arccos((np.trace(X)-1.)/2.)
	nx=np.sqrt((X[0,0]-np.cos(phi))/(1.-np.cos(phi)))
	ny=np.sqrt((X[1,1]-np.cos(phi))/(1.-np.cos(phi)))
	nz=np.sqrt((X[2,2]-np.cos(phi))/(1.-np.cos(phi)))
	x12 = lambda n1,n2,n3: n3*np.sin(phi)+n1*n2*(1.-np.cos(phi))
	x13 = lambda n1,n2,n3: -n2*np.sin(phi)+n1*n3*(1.-np.cos(phi))
	x23 = lambda n1,n2,n3: n1*np.sin(phi)+n2*n3*(1.-np.cos(phi))
	#print x13(nx,ny,nz), X[0,2]
	#for i in range(3):
	#	print norm(X[i,:])
	if simeq(X[0,1],x12(nx,ny,nz)) and simeq(X[0,2],x13(nx,ny,nz)) and simeq(X[1,2],x23(nx,ny,nz)):
		#print 1
		pass
	elif simeq(X[0,1],x12(nx,ny,-nz)) and simeq(X[0,2],x13(nx,ny,-nz)) and simeq(X[1,2],x23(nx,ny,-nz)):
		nz*=-1
		#print 2
	elif simeq(X[0,1],x12(nx,-ny,nz)) and simeq(X[0,2],x13(nx,-ny,nz)) and simeq(X[1,2],x23(nx,-ny,nz)):
		ny*=-1
		#print 3
	elif simeq(X[0,1],x12(-nx,ny,nz)) and simeq(X[0,2],x13(-nx,ny,nz)) and simeq(X[1,2],x23(-nx,ny,nz)):
		nx*=-1
		#print 4
	elif simeq(X[0,1],x12(nx,-ny,-nz)) and simeq(X[0,2],x13(nx,-ny,-nz)) and simeq(X[1,2],x23(nx,-ny,-nz)):
		ny*=-1
		nz*=-1
		#print 5
	elif simeq(X[0,1],x12(-nx,ny,-nz)) and simeq(X[0,2],x13(-nx,ny,-nz)) and simeq(X[1,2],x23(-nx,ny,-nz)):
		nx*=-1
		nz*=-1
		#print 6
	elif simeq(X[0,1],x12(-nx,-ny,nz)) and simeq(X[0,2],x13(-nx,-ny,nz)) and simeq(X[1,2],x23(-nx,-ny,nz)):
		nx*=-1
		ny*=-1
		#print 7
	elif simeq(X[0,1],x12(-nx,-ny,-nz)) and simeq(X[0,2],x13(-nx,-ny,-nz)) and simeq(X[1,2],x23(-nx,-ny,-nz)):
		nx*=-1
		ny*=-1
		nz*=-1
		#print 8
	else:
		print 'error in ellipses'; exit()
	polydata=transformer(ball,xyz,scale*E,180./np.pi*phi,(nx,ny,nz))
	ballMapper=vtk.vtkPolyDataMapper()
	ballMapper.SetInput(polydata)
	ballActor=vtk.vtkActor()
	ballActor.SetMapper(ballMapper)
	ballActor.GetProperty().SetColor(0.,0.,.5)
	return ballActor

def drawBond(bond,darken=False):
	scale=.07
	bnd=vtk.vtkCylinderSource()
	bnd.SetRadius(.1)
	bnd.SetResolution(50)
	at1=np.zeros((3))
	at2=np.zeros((3))
	At1=bond.GetBeginAtom()
	At2=bond.GetEndAtom()
	at1[:]=[At1.GetX(),At1.GetY(),At1.GetZ()]
	at2[:]=[At2.GetX(),At2.GetY(),At2.GetZ()]
	rvdW1=scale*table.GetVdwRad(At1.GetAtomicNum())
	rvdW2=scale*table.GetVdwRad(At2.GetAtomicNum())
	len=norm((at1-at2))-rvdW1-rvdW2
	bnd.SetHeight(len)
	n=(at1-at2)/norm((at1-at2))
	#print n
	pos=at2+n*len/2.+n*rvdW2#(at1+at2)/2.
	bndMapper=vtk.vtkPolyDataMapper()
	bndMapper.SetInput(bnd.GetOutput())
	bndActor=vtk.vtkActor()
	bndActor.SetMapper(bndMapper)
	if darken:
		bndActor.GetProperty().SetColor(.95,.95,.95)
	else:
		bndActor.GetProperty().SetColor(1.,1.,1.)
	bndActor=orientAct(bndActor,*n)
	bndActor.SetPosition(pos)
	return bndActor

def drawArrow(x123,l123,et,scale=.5):
	l=norm(l123)
	n=l123/l
	scale=l*scale
	arr=vtk.vtkArrowSource()
	arr=vtk.vtkArrowSource()
	arr.SetShaftResolution(50)
	arr.SetTipResolution(50)
	arr.SetTipRadius(.1/scale)
	arr.SetShaftRadius(.04/scale)
	arr.SetTipLength(.2/scale)
	arrMapper=vtk.vtkPolyDataMapper()
	arrMapper.SetInput(arr.GetOutput())
	arrActor=vtk.vtkActor()
	arrActor.SetMapper(arrMapper)
	arrActor.SetScale(scale)
	if et=="e":
		arrActor.GetProperty().SetColor((1,.37,.16))
	elif et=="e1":
		arrActor.GetProperty().SetColor((1,.25,.10))
	elif et=="e2":
		arrActor.GetProperty().SetColor((1,.4,.22))
	elif et in ("d","m"):
		arrActor.GetProperty().SetColor((.6,0.,.6))
	elif et in ("d1","m1"):
		arrActor.GetProperty().SetColor((.6,0.1,.52))
	elif et in ("d2","m2"):
		arrActor.GetProperty().SetColor((.6,0.25,.6))
	elif et in ("d3","m3"):
		arrActor.GetProperty().SetColor((.6,0.4,.68))
	elif et is "b":
		arrActor.GetProperty().SetColor((0.,.64,0.))
	else:
		arrActor.GetProperty().SetColor(.5,.5,.5)
	arrActor=orientArrow(arrActor,*n)
	arrActor.SetPosition(*x123)
	return arrActor

def drawSample(pos,color,scale=.5,len=1.):
	arr=vtk.vtkArrowSource()
	arr=vtk.vtkArrowSource()
	arr.SetShaftResolution(50)
	arr.SetTipResolution(50)
	arr.SetTipRadius(.1/scale/len)
	arr.SetShaftRadius(.04/scale/len)
	arr.SetTipLength(.2/scale/len)
	arrMapper=vtk.vtkPolyDataMapper()
	arrMapper.SetInput(arr.GetOutput())
	arrFollower=vtk.vtkFollower()
	arrFollower.SetMapper(arrMapper)
	arrFollower=orientArrow(arrFollower,1,0,0)
	arrFollower.SetPosition(pos)
	arrFollower.SetScale(scale*len)
	arrFollower.GetProperty().SetColor(color)
	return arrFollower

def genLine(start,end,color="none"):
	line=vtk.vtkLineSource()
	line.SetPoint1(start)
	line.SetPoint2(end)
	lineMapper=vtk.vtkPolyDataMapper()
	lineMapper.SetInput(line.GetOutput())
	lineActor=vtk.vtkActor()
	lineActor.SetMapper(lineMapper)
	if color is not "none":
		lineActor.GetProperty().SetColor(color)
	return lineActor

def orientAct(actor,n1,n2,n3):
	phi1,phi1=(0.,0.)
	if abs(n2)>.0001 or abs(n3)>.0001:
		phi1=np.arccos(n2/np.sqrt(n2**2.+n3**2.))
		if n3<0:
			phi1*=-1
	phi2=np.arccos(n2*np.cos(phi1)+n3*np.sin(phi1))
	if n1>0:
		phi2*=-1
	phi1*=180./np.pi
	phi2*=180./np.pi
	#print phi1,phi2
	actor.RotateX(phi1)
	actor.RotateZ(phi2)
	return actor

def orientArrow(actor,n1,n2,n3):
	phi1,phi2=(0.,0.)
	if abs(n1)>.0001 or abs(n3)>.0001:
		phi1=np.arccos(n1/np.sqrt(n1**2.+n3**2.))
		if n3>0:
			phi1*=-1
	phi2=np.arccos(n1*np.cos(phi1)-n3*np.sin(phi1))
	if n2<0:
		phi2*=-1
	phi1*=180./np.pi
	phi2*=180./np.pi
	actor.RotateY(phi1)
	actor.RotateZ(phi2)
	return actor
	
def genText(what,pos,color="none",scale=.5):
	text=vtk.vtkVectorText()
	text.SetText(what)
	textMapper=vtk.vtkPolyDataMapper()
	textMapper.SetInput(text.GetOutput())
	textFollower=vtk.vtkFollower()
	textFollower.SetMapper(textMapper)
	textFollower.AddPosition(pos)
	if color=="none":
		textFollower.GetProperty().SetColor((1,1,1))
	else:
		textFollower.GetProperty().SetColor(color)
	textFollower.SetScale(scale)
	return textFollower

def czytArr(fin):
	fin=fin.readlines()
	pos,arr,et=[[],[],[]]
	for i in range(len(fin)):
		temp1=np.zeros(3)
		temp2=np.zeros(3)
		line=fin[i].split()
		try:
			temp1[:]=[float(u) for u in line[:3]]
			temp2[:]=[float(u) for u in line[3:6]]
		except IndexError:
			print "end of file?"
			break
		except ValueError:
			print "end of file?"
			break
		try:
			temp3=line[6]
		except IndexError:
			temp3="none"
		pos.append(temp1)
		arr.append(temp2)
		et.append(temp3)
	return pos,arr,et

renWin=vtk.vtkRenderWindow()
renWin.SetSize(500,500)
ren1=vtk.vtkRenderer()
ren1.SetBackground(schema[0])

iren=vtk.vtkRenderWindowInteractor()
iren.SetRenderWindow(renWin)

#mol_objects=vtk.vtkAppendPolyData()
#print mol_objects.GetOutput()

if options.schema=="2":
	for i in range(mol.NumAtoms()):
		at=drawAtom(mol.GetAtom(i+1),True)
		ren1.AddActor(at)
		#mol_objects.AddInput(at.GetMapper().GetInput())
	if not options.nobonds:
		for i in range(mol.NumBonds()):
			bond=drawBond(mol.GetBond(i),True)
			ren1.AddActor(bond)
			#mol_objects.AddInput(bond.GetMapper().GetInput())
else:
	for i in range(mol.NumAtoms()):
		ren1.AddActor(drawAtom(mol.GetAtom(i+1)))
	if not options.nobonds:
		for i in range(mol.NumBonds()):
			ren1.AddActor(drawBond(mol.GetBond(i)))

#print mol_objects.GetOutput();exit()
if options.arrows is not "none":
	(xyz,magn,et)=czytArr(open(options.arrows))
	for i in range(len(xyz)):
		ren1.AddActor(drawArrow(xyz[i],magn[i],et[i],scales[no[et[i]]]))
	diptext=genText("dipole / 4 D",(-8.5,-.4,0.),scale=.25,color=schema[1])
	dip=drawSample((-8.5,-.5,0.),(.6,0.,.6),scale=scales[1],len=4)
	eltext=genText("field / 4 GV/m",(-8.5,-1.4,0.),scale=.25,color=schema[1])
	elf=drawSample((-8.5,-1.5,0.),(1,.37,.16),scale=scales[2],len=4)
	bettext=genText("beta / 2000 au",(-8.5,-.9,0.),scale=.25,color=schema[1])
	bet=drawSample((-8.5,-1.,0.),(0.,.64,0.),scale=scales[0],len=2000.)
	ren1.AddActor(dip)
	dip.SetCamera(ren1.GetActiveCamera())
	ren1.AddActor(elf)
	elf.SetCamera(ren1.GetActiveCamera())
	ren1.AddActor(bet)
	bet.SetCamera(ren1.GetActiveCamera())
	ren1.AddActor(diptext)
	diptext.SetCamera(ren1.GetActiveCamera())
	ren1.AddActor(eltext)
	eltext.SetCamera(ren1.GetActiveCamera())
	ren1.AddActor(bettext)
	bettext.SetCamera(ren1.GetActiveCamera())
	camera=ren1.GetActiveCamera()
	camera.SetParallelProjection(True)
	
if options.arrows2 is not "none":
	(xyz,magn,et)=czytArr(open(options.arrows2))
	for i in range(len(xyz)):
		ren1.AddActor(drawArrow(xyz[i],magn[i],et[i],scales[no[et[i]]]))
	diptext=genText("dipole / 4 D",(-8.5,-.4,0.),scale=.25,color=schema[1])
	dip=drawSample((-8.5,-.5,0.),(.6,0.,.6),scale=scales[1],len=2)
	eltext=genText("field / 0.5 GV/m",(-8.5,-1.4,0.),scale=.25,color=schema[1])
	elf=drawSample((-8.5,-1.5,0.),(1,.37,.16),scale=scales[2],len=.5)
	#bettext=genText("beta / 2000 au",(-8.5,-.9,0.),scale=.25,color=schema[1])
	#bet=drawSample((-8.5,-1.,0.),(0.,.64,0.),scale=scales[0],len=2000.)
	ren1.AddActor(dip)
	dip.SetCamera(ren1.GetActiveCamera())
	ren1.AddActor(elf)
	elf.SetCamera(ren1.GetActiveCamera())
	#ren1.AddActor(bet)
	#bet.SetCamera(ren1.GetActiveCamera())
	ren1.AddActor(diptext)
	diptext.SetCamera(ren1.GetActiveCamera())
	ren1.AddActor(eltext)
	eltext.SetCamera(ren1.GetActiveCamera())
	#ren1.AddActor(bettext)
	#bettext.SetCamera(ren1.GetActiveCamera())
	camera=ren1.GetActiveCamera()
	camera.SetParallelProjection(True)

if options.cell:
	par=open(options.inputfile).readlines()[1].split()
	try:
		if options.prime:
			a,b,c=[float(u) for u in par[2:5]]
			al,be,ga=[float(u)*np.pi/180. for u in par[5:8]]
		else:
			a,b,c=[float(u) for u in par[:3]]
			al,be,ga=[float(u)*np.pi/180. for u in par[3:6]]
	except IndexError:
		parser.error("lattice parameters not found or incorrect")
	except ValueError:
		parser.error("lattice parameters not found or incorrect")
		print par
	csa=np.cos(al)
	csb=np.cos(be)
	csc=np.cos(ga)
	snc=np.sin(ga)
	G=np.zeros((3,3))
	G[0,0]=a**2.
	G[1,1]=b**2.
	G[2,2]=c**2.
	G[0,1]=a*b*csc
	G[1,0]=G[0,1]
	G[0,2]=a*c*csb
	G[2,0]=G[0,2]
	G[1,2]=b*c*csa
	G[2,1]=G[1,2]
	G_r=inv(G)
	V=np.sqrt(float(det(G)))
	X=np.zeros((3,3))
	X[0,0]=1./a/snc
	X[1,0]=-csc/b/snc
	X[1,1]=1./b
	X[:,2]=V/a/b/snc*G_r[:,2]
	X_inv=inv(X).T
	def transf(l,X_inv):
		return np.dot(l,X_inv)
	ren1.AddActor(genLine((0,0,0),transf((1,0,0),X_inv),(1,0,0)))
	ren1.AddActor(genLine((0,0,0),transf((0,1,0),X_inv),(0,1,0)))
	ren1.AddActor(genLine((0,0,0),transf((0,0,1),X_inv),(0,0,1)))
	ren1.AddActor(genLine(transf((1,1,1),X_inv),transf((1,1,0),X_inv),schema[1]))
	ren1.AddActor(genLine(transf((1,1,1),X_inv),transf((1,0,1),X_inv),schema[1]))
	ren1.AddActor(genLine(transf((1,1,1),X_inv),transf((0,1,1),X_inv),schema[1]))
	ren1.AddActor(genLine(transf((1,0,0),X_inv),transf((1,0,1),X_inv),schema[1]))
	ren1.AddActor(genLine(transf((1,0,0),X_inv),transf((1,1,0),X_inv),schema[1]))
	ren1.AddActor(genLine(transf((0,1,0),X_inv),transf((0,1,1),X_inv),schema[1]))
	ren1.AddActor(genLine(transf((0,1,0),X_inv),transf((1,1,0),X_inv),schema[1]))
	ren1.AddActor(genLine(transf((0,0,1),X_inv),transf((1,0,1),X_inv),schema[1]))
	ren1.AddActor(genLine(transf((0,0,1),X_inv),transf((0,1,1),X_inv),schema[1]))
	sc1=1.2
	sc2=.2
	lbl_a=genText("a",transf((sc1*.9,sc2*-.1,sc2*-.1),X_inv),schema[1])
	lbl_b=genText("b",transf((sc2*-.1,sc1*.9,sc2*-.1),X_inv),schema[1])
	lbl_c=genText("c",transf((sc2*-.1,sc2*-.1,sc1*.9),X_inv),schema[1])
	for act in (lbl_a, lbl_b, lbl_c):
		ren1.AddActor(act)
		act.SetCamera(ren1.GetActiveCamera())
	camera=ren1.GetActiveCamera()
	camera.SetParallelProjection(True)
	camera.SetViewUp(transf((0,1,0),X_inv))

if options.axes is not "none":
	pos=[float(u) for u in options.axes.split()[:3]]
	pos=np.array(pos)
	xax=drawArrow(pos,[2,0,0],'x',1)
	xtx=genText('X',pos+[2.2,-.2,-.2],scale=.4,color=schema[1])
	#xtx=genText('X',pos+[2.2,-.4,-.2],scale=.4,color=schema[1])
	yax=drawArrow(pos,[0,2,0],'x',1)
	ytx=genText('Y',pos+[-.2,2.5,-.2],scale=.4,color=schema[1])
	#ytx=genText('Y',pos+[-.2,2.2,-.2],scale=.4,color=schema[1])
	zax=drawArrow(pos,[0,0,2],'x',1)
	ztx=genText('Z',pos+[-.2,-.2,3.2],scale=.4,color=schema[1])
	#ztx=genText('Z',pos+[-.2,-.2,2.2],scale=.4,color=schema[1])
	for ac in (xax,xtx,yax,ytx,zax,ztx):
		ren1.AddActor(ac)
		if ac in (xtx,ytx,ztx):
			ac.SetCamera(ren1.GetActiveCamera())

if options.polar is not "none":
	f_pol=open(options.polar).readlines()
	if len(f_pol[0].split())==9:
		for i in range(1,mol.NumAtoms()+1):
			atom=mol.GetAtom(i)
			x,y,z=(atom.GetX(),atom.GetY(),atom.GetZ())
			a11,a22,a33,a12,a13,a23=[float(u) for u in f_pol[i-1].split()[3:9]]
			alf=np.zeros((3,3))
			alf[0,0]=a11
			alf[1,1]=a22
			alf[2,2]=a33
			alf[0,1]=a12
			alf[1,0]=a12
			alf[0,2]=a13
			alf[2,0]=a13
			alf[1,2]=a23
			alf[2,1]=a23
			ren1.AddActor(genEllipse((x,y,z),alf))
	else:
		for i in range(len(f_pol)):
			line=f_pol[i].split()
			x,y,z=[float(u) for u in line[:3]]
			a11,a22,a33,a12,a13,a23=[float(u) for u in f_pol[i].split()[6:12]]
			alf=np.zeros((3,3))
			alf[0,0]=a11
			alf[0,0]=a11
			alf[1,1]=a22
			alf[2,2]=a33
			alf[0,1]=a12
			alf[1,0]=a12
			alf[0,2]=a13
			alf[2,0]=a13
			alf[1,2]=a23
			alf[2,1]=a23
			ren1.AddActor(genEllipse((x,y,z),alf))#*5))
#exit()

#ren1.ResetCameraClippingRange()
ren1.ResetCamera()
camera=ren1.GetActiveCamera()
camera.SetParallelProjection(True)
#ren1.RemoveAllLights()
camera.Azimuth(90)
#camera.SetViewUp((0,0,1))
#print camera.GetFocalPoint()
#print camera.GetOrientationWXYZ()
#print camera.GetPosition()
#camera.SetViewAngle(90)
#print camera.GetViewAngle()

def capt(exp,fname,i):
	fn=fname+"_"+str(i)
	exp.SetFilePrefix(fn)
	exp.Write()
	print "screen no", i, "captured"
	i+=1
	return i

if options.vector is not "none":
	i=0
	fname=options.vector
	exp=vtk.vtkGL2PSExporter()#	FileName=options.vector+"svg")
	#exp.SetFileFormatToSVG()
	exp.SetFileFormatToEPS()
	#exp.SetFileFormatToPDF()
	exp.DrawBackgroundOff()
	exp.SimpleLineOffsetOff()
	#exp.Write3DPropsAsRasterImageOn()
	#exp.OcclusionCullOff()
	exp.SetSortToSimple()
	#exp.BestRootOff()
	#exp.PS3ShadingOff()
	#exp.SimpleLineOffsetOn()
	exp.SetRenderWindow(renWin)
	exp.SetFilePrefix(options.vector)
	#exp.Write()
	##root=tk.Tk()
	picker=vtk.vtkPicker()
	def annotatePick(object,event):
		global picker,capture,i
		i=capt(exp,fname,i)
	picker.AddObserver("EndPickEvent",annotatePick)
	iren.SetPicker(picker)
	renWin.AddRenderer(ren1)
	renWin.Render()
	iren.Start()
	exp.SetFilePrefix(fname)
	exp.Write()
	##root.bind_all('<Return>',capt(exp,fname,i))
	#root.withdraw()
	##root.mainloop()
	#fout=open(options.vector+".eps","w")
	#fout.close()
	#w=vtk.vtkXMLPolyDataWriter()
	#w.SetInput(mol_objects.GetOutput())
	#filename=options.vector+".vtp"#w.GetDefaultFileExtension()
	#w.SetFileName(filename)
	#w.Write()
else:
	renWin.AddRenderer(ren1)
	renWin.Render()
	iren.Start()
