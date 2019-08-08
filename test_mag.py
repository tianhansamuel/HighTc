from numpy import *
from pylab import *
from math import *
import matplotlib
import numpy as np
from matplotlib import cm
from matplotlib import pyplot as plt
from numpy import linalg as LA
from math import exp


NN=500
NB=10000
N1=100

J1=-1*0.67
Delta=4*0.67
GG=40
dt=20
seuil=0.1

ff1=open('dd.dat','w')

ff2=open('aa.dat','w')

MM=zeros((2,2),complex)

#############################################  bare propagator 

def om(alpha,dd):
    EG=-(1-4*alpha**2)*(1-dd)**2/((1+4*alpha**2))+dd*(1-dd)**2*(7-7*alpha)/(2*(1+4*alpha**2))-dd*(1-dd)**2*7*alpha/(1+4*alpha**2)
    return EG

def om1(alpha,dd):
    EG=-(1-4*alpha**2)*(1-dd)**2/((1+4*alpha**2))+dd*(1-dd)**2*(3-3*alpha)/(2*(1+4*alpha**2))-dd*(1-dd)**2*3*alpha/(1+4*alpha**2)
    return EG


def AMAX(dd):
    EG=100
    aa=0
    for ia in range(NB):
       alpha=0.3*ia/NB
       if om(alpha,dd)<EG:
          aa=alpha
          EG=om(alpha,dd)
    return (aa,EG)

def AMAX1(dd):
    EG=100
    aa=0
    for ia in range(NB):
       alpha=0.2*ia/NB
       if om1(alpha,dd)<EG:
          aa=alpha
          EG=om1(alpha,dd)
    return (aa,EG)


for ia in range(NN):
    aa=0.25*ia/NN
    E0=om(aa,0.1)
    E1=om(aa,0.2)
    E2=om(aa,0.3)
    E3=om(aa,0.4)
    E4=om(aa,0.5)
    EE0=om1(aa,0.1)
    EE1=om1(aa,0.2)
    EE2=om1(aa,0.3)
    EE3=om1(aa,0.4)
    EE4=om1(aa,0.5)
    ff1.write(str(aa)+'\t')
    ff1.write(str(E0)+'\t')
    ff1.write(str(E1)+'\t')
    ff1.write(str(E2)+'\t')
    ff1.write(str(E3)+'\t')
    ff1.write(str(E4)+'\t')
    ff1.write(str(EE0)+'\t')
    ff1.write(str(EE1)+'\t')
    ff1.write(str(EE2)+'\t')
    ff1.write(str(EE3)+'\t')
    ff1.write(str(EE4)+'\n')



for ii in range(NN):
   de=1.0*ii/NN
   am=AMAX(de)[0]
   am1=AMAX1(de)[0]
   ff2.write(str(de)+'\t')
   ff2.write(str(am)+'\t')
   ff2.write(str(am1)+'\n')

