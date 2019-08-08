from numpy import *
from math import *


NN=300
NB=300
NN1=600
NR=100

J1=-1*0.67
Delta=4*0.67
GG=40
dt=20
seuil=0.1

ff2=open('tc3.dat','w')

MM=zeros((2,2),complex)

#############################################  bare propagator 


def tb(px,py,t,t1,t2,t3,t4):
    EE=0.5*t*(cos(px)+cos(py))+t1*cos(px)*cos(py)+0.5*t2*(cos(2*px)+cos(2*py))+0.5*t3*(cos(2*px)*cos(py)+cos(px)*cos(2*py))+t4*cos(2*px)*cos(2*py)
    return EE

def DE(px,py,t,t1,t2,t3,t4,ww):
    omega=tb(px+pi,py+pi,t,t1,t2,t3,t4)-tb(px,py,t,t1,t2,t3,t4)+ww
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

def TANH(e,Temp):
    if e!=0:
       cc=tanh(e*1.0/Temp)/e
    else: cc=1.0/Temp
    return cc

def GAMMA(kx,ky,T,t,t1,t2,t3,t4,aa1,ww,dop):
    DF=DE(kx,ky,t,t1,t2,t3,t4,ww)
    xi=tb(kx,ky,t,t1,t2,t3,t4)
    xi1=tb(kx+pi,ky+pi,t,t1,t2,t3,t4)
    if DF**2-(xi1-xi)**2!=0:
     SS=0.25*(cos(kx)-cos(ky))**2*(1-12*aa1**2)*pi/(1+12*aa1**2)*(1-dop)**2*xi1**2*DF/(DF**2-(xi1-xi)**2)*(TANH(DF,T)-TANH(xi1-xi,T))
    else: SS=-0.25*(cos(kx)-cos(ky))**2*(1-12*aa1**2)*pi/(1+12*aa1**2)*(1-dop)**2*xi1**2*DF/(DF+(xi1-xi))*1.0/(DF*(xi1-xi))
    return SS


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

def CCT(EF,T,t,t1,t2,t3,t4,ww,aa,dop):
    NF=0
    CT=0
    for ix in range(NN):
      for iy in range(NN):
        kx=-pi+2*pi*ix/NN
        ky=-pi+2*pi*iy/NN    
        if abs(tb(kx,ky,t,t1,t2,t3,t4)-EF)<0.03 and tb(kx,ky,t,t1,t2,t3,t4)-tb(kx+pi,ky+pi,t,t1,t2,t3,t4)<ww:
           NF+=1
           CT+=GAMMA(kx, ky, T, t,t1,t2,t3,t4,aa,ww,dop)
    CC=CT*1.0/NF
    return CC


def TCD(EF,t,t1,t2,t3,t4,ww,aa,dop):
    tc=0
    for tt1 in range(NN1):
       TTC=0.3-0.3*tt1/NN1
       if abs(CCT(EF,TTC,t,t1,t2,t3,t4,ww,aa,dop))>1:
         tc=TTC
         break
    return tc


for ii in range(NN):
   print ii
   EF=-0.2*ii/NN
   dop=2*(0.5-dd(EF,-0.5908,0.0962,-0.1306,-0.0507,0.0939)+0.1)
   ww1=AMAX(dop)[1]
   aa1=AMAX(dop)[0]
   tc=TCD(EF,-0.5908,0.0962,-0.1306,-0.0507,0.0939,ww1,aa1,dop)
   ff2.write(str(EF)+'\t')
   ff2.write(str(dop)+'\t')
   ff2.write(str(tc)+'\n')


