from numpy import *
from pylab import *
from math import *
import matplotlib
import numpy as np
from matplotlib import cm
from matplotlib import pyplot as plt
from numpy import linalg as LA
from math import exp


NN=1000
NB=100
DL=5.0
ff1=open('delta.dat','w')

ff2=open('criticality.dat','w')

ff3=open('criticality1.dat','w')

dd=1


def tb(px,py,t,t1,t2,t3,t4):
    EE=0.5*t*(cos(px)+cos(py))+t1*cos(px)*cos(py)+0.5*t2*(cos(2*px)+cos(2*py))+0.5*t3*(cos(2*px)*cos(py)+cos(px)*cos(2*py))+t4*cos(2*px)*cos(2*py)
    return EE

def RE(qx,qy,J):
   OM=2*J*(2-cos(qx)-cos(qy))
   return OM

def dd(EF):
    SS=0
    for ix in range(NN):
      for iy in range(NN):
         px=2*pi*ix/NN-pi
         py=2*pi*iy/NN-pi
         if tb(px,py,-0.5908,0.0962,-0.1306,-0.0507,0.0939)<EF:
            SS+=1.0/(NN**2)
    return SS

def EE(px,py,qx,qy,J):
    E1=tb(px+qx,py+qy,-0.5908,0.0962,-0.1306,-0.0507,0.0939)+RE(qx,qy,J)
    return E1
      



for ix in range(NB+1):
   for iy in range(NB+1):
    qx=-pi+2*pi*ix/NB
    qy=-pi+2*pi*iy/NB
    ff2.write(str(EE(0,4*pi/5,qx,qy,0.2)-tb(0,4*pi/5,-0.5908,0.0962,-0.1306,-0.0507,0.0939))+'\t')
    ff3.write(str(EE(pi/2,pi/2,qx,qy,0.2)-tb(pi/2,pi/2,-0.5908,0.0962,-0.1306,-0.0507,0.0939))+'\t')
   ff2.write('\n')
   ff3.write('\n')


