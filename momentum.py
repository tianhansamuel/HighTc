from numpy import *
from math import *


NN=300
N1=100
NB=200
ff1=open('mm100_04.dat','w')

dd=1



def tb(px,py,t,t1,t2,t3,t4):
    EE=0.5*t*(cos(px)+cos(py))+t1*cos(px)*cos(py)+0.5*t2*(cos(2*px)+cos(2*py))+0.5*t3*(cos(2*px)*cos(py)+cos(px)*cos(2*py))+t4*cos(2*px)*cos(2*py)
    return EE

def intense(qx,qy,t,t1,t2,t3,t4,EF):
  ww=0.04
  SS=0
  for ix in range(N1):
    for iy in range(N1):
     px=2*ix*pi/N1
     py=2*iy*pi/N1
     if tb(px,py,t,t1,t2,t3,t4)-tb(px+pi,py+pi,t,t1,t2,t3,t4)>ww and abs(tb(px,py,t,t1,t2,t3,t4)-EF)<0.05 and tb(px,py,t,t1,t2,t3,t4)-tb(px+pi+qx,py+pi+qy,t,t1,t2,t3,t4)>ww and tb(px+qx,py+qy,t,t1,t2,t3,t4)<tb(px,py,t,t1,t2,t3,t4):
       SS+=1
  return SS*1.0

def DE(px,py,t,t1,t2,t3,t4,nd):
    omega=tb(px+pi,py+pi,t,t1,t2,t3,t4)-tb(px,py,t,t1,t2,t3,t4)+AMAX(nd)[0]
    return omega


def om(alpha,dd):
    EG=-(1-12*alpha**2)/(1+12*alpha**2)*(1-dd)**2+dd*(1-dd)**4*3*(1-alpha)/(1+12*alpha**2)+8*dd*(1-dd)**5*(1-alpha)/(1+12*alpha**2)
    return EG

def AMAX(dd):
    EG=100
    aa=0
    for ia in range(NB):
       alpha=0.3*ia/NB
       if om(alpha,dd)<EG:
          aa=alpha
          EG=om(alpha,dd)
    RES=4*aa/(1-12*aa**2)*0.129
    return (RES,aa)



for iix in range(N1+1):
  print iix
  for iiy in range(N1+1):
   px=4*iix*pi/N1-2*pi
   py=4*iiy*pi/N1-2*pi
   gg=intense(px,py,-0.5908,0.0962,-0.1306,-0.0507,0.0939,-0.04)
   ff1.write(str(gg)+'\t')
  ff1.write('\n')


