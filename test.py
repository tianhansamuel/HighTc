from numpy import *
from pylab import *
from math import *
import matplotlib
import numpy as np
from matplotlib import cm
from matplotlib import pyplot as plt
from numpy import linalg as LA
from math import exp


NN=400
NB=15
N1=100

J1=-1*0.67
Delta=4*0.67
GG=40
dt=20
seuil=0.1

ff1=open('sc_temp.dat','w')

ff2=open('EE.dat','w')

MM=zeros((2,2),complex)

#############################################  bare propagator 


def tb(px,py):
    tb=-2*(cos(px)+cos(py))+cos(px)*cos(py)
    return tb

def TM(alpha,kx):
   MM[0,0]=1
   MM[0,1]=alpha*(cos(kx)+1j*sin(kx))
   MM[1,0]=alpha*(cos(kx)-1j*sin(kx))
   MM[1,1]=alpha**2
   return MM

def SS(qx,qy):
   sq=120*exp(-2*((qx-pi)**2+(qy-pi)**2))
   return sq

S1=0
for ix in range(NB):
   qx=2*pi*ix/NB 
   for iy in range(NB):
    qy=2*pi*iy/NB 
    S1+=SS(qx,qy)/(NB**2)
print S1

def G0(omega,kx,ky,T,alpha):
    GG=0
    for it in range(NN):
     om=(2*it+1)*pi*T
     Sig=0
     for ia in range(N1):
      GA=alpha*(1.5-1.0*ia/N1)
      Sig+=tb(kx+pi,ky+pi)**2*(1-abs(ia/N1-0.5))/(omega-1j*om-tb(kx+pi,ky+pi)-GA)/N1
     if tb(kx+pi,ky+pi)<=tb(kx,ky):
        GG+=1.0/(omega-tb(kx,ky)-1j*om-Sig)*1.0/NN
     else:        GG+=1.0/(omega-tb(kx,ky)-1j*om)*1.0/NN
    return GG



for io in range(NN):
    omega=-5+10.0*io/NN
    print io
    ff1.write(str(omega)+'\t')
    ff1.write(str(G0(omega,pi,0,0.01,0.2).real)+'\t')
    ff1.write(str(G0(omega,pi,0,0.01,0.2).imag)+'\t')
    ff1.write(str(G0(omega,pi,0,0.01,0).real)+'\t')
    ff1.write(str(G0(omega,pi,0,0.01,0).imag)+'\t')
    ff1.write(str(G0(omega,pi/2,pi/2,0.01,0.1).real)+'\t')
    ff1.write(str(G0(omega,pi/2,pi/2,0.01,0.1).imag)+'\n')
