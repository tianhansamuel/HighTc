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

ff1=open('dd_new.dat','w')

ff2=open('aa_new.dat','w')

MM=zeros((2,2),complex)

#############################################  bare propagator 


def tb(px,py):
    tb=-2*(cos(px)+cos(py))+cos(px)*cos(py)
    return tb

def om(alpha,dd):
    EG1=-(1-4*alpha**2)*(1-dd)**2/((1+4*alpha**2))
    EG2=alpha**2*dd*((1-dd)**3*3*alpha/(1+4*alpha**2)+3*dd*(1-dd)**2*2*alpha/(1+4*alpha**2)+3*dd**2*(1-dd)*alpha/(1+4*alpha**2))
    EG3=alpha**2*dd*((1-dd)**4*4*alpha/(1+4*alpha**2)+4*dd*(1-dd)**3*3*alpha/(1+4*alpha**2)+6*dd**2*(1-dd)**2*2*alpha/(1+4*alpha**2)+4*dd**3*(1-dd)*alpha/(1+4*alpha**2))*2
    EG4=dd*((1-dd)**3*3*(1-alpha)/(1+4*alpha**2)+3*dd*(1-dd)**2*2*(1-alpha)/(1+4*alpha**2)+3*dd**2*(1-dd)*(1-alpha)/(1+4*alpha**2))
    EG5=dd*((1-dd)**4*4*(1-alpha)/(1+4*alpha**2)+4*dd*(1-dd)**3*3*(1-alpha)/(1+4*alpha**2)+6*dd**2*(1-dd)**2*2*(1-alpha)/(1+4*alpha**2)+4*dd**3*(1-dd)*(1-alpha)/(1+4*alpha**2))*2
    EG=EG1+EG2+EG4
    return EG



def AMAX(dd):
    EG=100
    aa=0
    for ia in range(NB):
       alpha=1.0*ia/NB
       if om(alpha,dd)<EG:
          aa=alpha
          EG=om(alpha,dd)
    return (aa,EG)

def AMAX1(dd):
    EG=100
    aa=0
    for ia in range(NB):
       alpha=0.5*ia/NB
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
    ff1.write(str(aa)+'\t')
    ff1.write(str(E0)+'\t')
    ff1.write(str(E1)+'\t')
    ff1.write(str(E2)+'\t')
    ff1.write(str(E3)+'\t')
    ff1.write(str(E4)+'\n')



for ii in range(NN):
   de=1.0*ii/NN
   am=AMAX(de)[0]
   ff2.write(str(de)+'\t')
   ff2.write(str(am)+'\n')

