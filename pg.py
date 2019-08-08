from numpy import *
from math import *


NN=300
N1=100
NB=1000
ff1=open('pseudo1.dat','w')

dd=1



def tb(px,py,t,t1,t2,t3,t4):
    EE=0.5*t*(cos(px)+cos(py))+t1*cos(px)*cos(py)+0.5*t2*(cos(2*px)+cos(2*py))+0.5*t3*(cos(2*px)*cos(py)+cos(px)*cos(2*py))+t4*cos(2*px)*cos(2*py)
    return EE

def intense(t,t1,t2,t3,t4,EF):
  ww=0.04
  dd=10
  kx=0
  ky=0
  for ix in range(NB):
    for iy in range(NB):
     px=0.5*ix*pi/NB
     py=iy*pi/NB
     if  abs(tb(px,py,t,t1,t2,t3,t4)-EF)+abs(tb(px,py,t,t1,t2,t3,t4)-tb(px+pi,py+pi,t,t1,t2,t3,t4)-ww)<dd:
       dd=abs(tb(px,py,t,t1,t2,t3,t4)-EF)+abs(tb(px,py,t,t1,t2,t3,t4)-tb(px+pi,py+pi,t,t1,t2,t3,t4)-ww)
       kx=px
       ky=py
  return (kx,ky)


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


def dd(EF,t,t1,t2,t3,t4):
      NF=0
      for ix in range(NB+1):
        for iy in range(NB+1):
         kx=-pi+2*pi*ix/NB
         ky=-pi+2*pi*iy/NB
         if tb(kx,ky,t,t1,t2,t3,t4)<=EF:
           NF+=1
      dl=NF*1.0/(NB+1)**2
      return dl

for iix in range(NN+1):
   print iix
   EE=-0.2*iix/NN
   DEL=dd(EE,-0.5908,0.0962,-0.1306,-0.0507,0.0939)
   kx=intense(-0.5908,0.0962,-0.1306,-0.0507,0.0939,EE)[0]
   ky=intense(-0.5908,0.0962,-0.1306,-0.0507,0.0939,EE)[1]
   dop=2*(0.5-DEL+0.1)
   aa=AMAX(dop)[1]
   PG=0.129*aa*(cos(pi+2*kx)+1)/((1-12*aa**2))
   ff1.write(str(EE)+'\t')
   ff1.write(str(DEL)+'\t')
   ff1.write(str(kx)+'\t')
   ff1.write(str(ky)+'\t')
   ff1.write(str(PG)+'\n')

