from numpy import *
from math import *


NN=30
NB=200
DL=5.0
ff1=open('FS12.dat','w')


def tb(px,py,t,t1,t2,t3,t4):
    EE=0.5*t*(cos(px)+cos(py))+t1*cos(px)*cos(py)+0.5*t2*(cos(2*px)+cos(2*py))+0.5*t3*(cos(2*px)*cos(py)+cos(px)*cos(2*py))+t4*cos(2*px)*cos(2*py)
    return EE

def OM(qx,qy):
    EO=0.04*(0.5-0.25*(cos(qx)+cos(qy)))
    return EO

def OM1(qx,qy):
    EO1=0.04*(1.5+0.25*(cos(qx)+cos(qy)))
    return EO1


def CC(px,py,t,t1,t2,t3,t4):
     c1=0
     if -tb(px,py,t,t1,t2,t3,t4)>=-tb(px+pi,py+pi,t,t1,t2,t3,t4)+0.04:
              c1=1
     return c1

for ix in range(NB+1):
  print ix
  for iy in range(NB+1):
      kx=-pi+2*pi*ix/NB
      ky=-pi+2*pi*iy/NB
      CT=CC(kx,ky,-0.5908,0.0962,-0.1306,-0.0507,0.0939)
      ff1.write(str(CT)+'\t')
  ff1.write('\n')



