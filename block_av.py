#!/usr/bin/python

import numpy as np
import matplotlib.pyplot as plt
import sys

def calc_sgm(x):
	sx=0.
	sxkw=0.
	n=len(x)
	for i in range(n):
		sx+=x[i]
		sxkw+=x[i]**2.
	sx/=n
	sxkw/=n
	sgm=np.sqrt((sxkw-sx**2.)/(n-1))
	sgmerr=sgm/np.sqrt(2*(n-1))
	return (sgm,sgmerr)

def block(x):
	n=len(x)
	n2=n//2
	xnw=[]
	for i in range(n2):
		xnw.append((x[2*i]+x[2*i+1])/2.)
	xnw=np.array(xnw)
	return xnw

fin=open(sys.argv[1]).readlines()
#psign=int(sys.argv[2])

n=len(fin)
try:
	nav=int(sys.argv[2])
except IndexError:
	print "average over all points"
	nav=n
except TypeError:
	print "average over all points"
	nav=n

x=[]
if nav>n:
	nav=n
	print "uwaga n=", nav
for i in range(n-nav,n):
	x.append(float(fin[i].split()[1]))
#I=0
#for i in range(n):
#	if int(fin[i].split()[0])<psign:
#		continue
#	else:
#		I=i
#		break
#for i in range(I,n):
#	x.append(float(fin[i].split()[1]))

x=np.array(x)
xavg=np.average(x)

kmax=int(np.log(nav)/np.log(2.))-2

ks=[]
sgm=[]
sgmerr=[]
for k in range(kmax):
	ks.append(k)
	sgm_sgmerr=calc_sgm(x)
	sgm.append(sgm_sgmerr[0])
	sgmerr.append(sgm_sgmerr[1])
	if k==0:
		#print "%.2f %.2f %.2f" % (xavg, sgm_sgmerr[0], sgm_sgmerr[1])
		print "%.2f %.2f" % (xavg, sgm_sgmerr[0])
	x=block(x)

#plt.errorbar(ks,sgm,yerr=sgmerr,fmt='b.')
#plt.show()

