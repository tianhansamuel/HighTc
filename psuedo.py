from numpy import *
from pylab import *
from math import *
import matplotlib
import numpy as np
from matplotlib import cm
from matplotlib import pyplot as plt
from numpy import linalg as LA
from math import exp


NN=5000
NB=200
DL=5.0
ff1=open('psdg.dat','w')

ff2=open('oomm.dat','w')

#############################################  bare propagator 

def GG(omega,DD,TT,la):
    GS=0
    for it in range(NB):
       om=(2*it+1)*TT
       GS+=1.0/(omega+1j*om-la**2/(omega-DD+1j*om))/NB
    return GS

def GAP(DD,TT,la):
   OO=0
   for im in range(NB):
     omega=(2.0*im/NB-1)*DL
     if (GG(omega,DD,TT,la).imag-GG(omega+DL/NB,DD,TT,la).imag)*(GG(omega-DL/NB,DD,TT,la).imag-GG(omega,DD,TT,la).imag)<0:
       OO=omega
   return OO

print GAP(10,1,1)
print GAP(10,5,1)


def TC(DD,la):
   tc=0
   for it in range(1,NN):
    TT=100.0*it/NN
    if GAP(DD,TT,la)<0.05:
       tc=TT
       break
   return tc

for ii in range(NN):
    omega=6.0*ii/NN-1
    G1=-GG(omega,1,0.1,1).imag
    G2=-GG(omega,1,0.3,1).imag
    G3=-GG(omega,1,0.5,1).imag
    G4=-GG(omega,1,0.8,1).imag
    G5=-GG(omega,1,1,1).imag
    G6=-GG(omega,1,1.2,1).imag
    ff1.write(str(omega)+'\t')
    ff1.write(str(G1)+'\t')
    ff1.write(str(G2)+'\t')
    ff1.write(str(G3)+'\t')
    ff1.write(str(G4)+'\t')
    ff1.write(str(G5)+'\t')
    ff1.write(str(G6)+'\n')

