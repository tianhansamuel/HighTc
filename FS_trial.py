from numpy import *
from pylab import *
from math import *
import matplotlib
import numpy as np
from matplotlib import cm
from matplotlib import pyplot as plt
from numpy import linalg as LA
from math import exp


NN=20
NB=400
DL=5.0
ff1=open('trial.dat','w')

def tb(px,py,t,t1,t2,t3,t4):
    EE=0.5*t*(cos(px)+cos(py))+t1*cos(px)*cos(py)+0.5*t2*(cos(2*px)+cos(2*py))+0.5*t3*(cos(2*px)*cos(py)+cos(px)*cos(2*py))+t4*cos(2*px)*cos(2*py)
    return EE

def vf(px,t,t1,t2,t3,t4):
    v1=(tb(px+0.001,pi-px-0.001,t,t1,t2,t3,t4)-tb(px,pi-px,t,t1,t2,t3,t4))/0.001
    return v1

def QM(kx):
    qm=0
    for ix in range(NB):
      qq=2*pi*ix/NB
      if xi(kx-qq)>xi(kx):
         qm=qq
         break
    return qm

for ix in range(NB+1):
   kx=pi*ix/NB-0.5*pi
   e1=tb(kx,pi-kx,-0.5908,0.0962,-0.1306,-0.0507,0.0939)-tb(0,pi,-0.5908,0.0962,-0.1306,-0.0507,0.0939)
   VV=vf(kx,-0.5908,0.0962,-0.1306,-0.0507,0.0939)
   ff1.write(str(kx/(pi))+'\t')
   ff1.write(str(VV)+'\t')
   ff1.write(str(e1)+'\n')



