from numpy import *
from math import *


NN=200
NB=100
NR=100

J1=-1*0.67
Delta=4*0.67
GG=40
dt=20
seuil=0.1

ff2=open('tc_topo.dat','w')


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
       alpha=0.1*ia/NB
       if om(alpha,dd)<EG:
          aa=alpha
          EG=om(alpha,dd)
    RES=4*aa/(1-12*aa**2)*0.129
    return (RES,aa)

def TT(e,Temp):
    if e!=0:
       cc=tanh(e*1.0/(2*Temp))/e
    else: cc=1.0/(2*Temp)
    return cc

def GAMMA(kx,ky,dd,T,t,t1,t2,t3,t4,aa,ww):
    DF=DE(kx,ky,t,t1,t2,t3,t4,ww)
    xi=tb(kx,ky,t,t1,t2,t3,t4)
    xip=tb(kx+pi,ky+pi,t,t1,t2,t3,t4)
    SS=0
    if ww**2-(xip-xi)**2!=0:
     SS=0.25*(cos(kx)-cos(ky))**2*(1-12*aa**2)/(1+12*aa**2)*(1-dd)**2*(xip)**2*ww/(ww**2-(xip-xi)**2)*(TT(ww,T)-TT(xip-xi,T))
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

def CCT(EF,T,t,t1,t2,t3,t4,dop,aa,ww):
    NF=0
    CT=0
    for ix in range(NN):
      for iy in range(NN):
        kx=-pi+2*pi*ix/NN
        ky=-pi+2*pi*iy/NN    
        if abs(tb(kx,ky,t,t1,t2,t3,t4)-EF)<0.01 and tb(kx,ky,t,t1,t2,t3,t4)-tb(kx+pi,ky+pi,t,t1,t2,t3,t4)<ww:
           NF+=1
           CT+=GAMMA(kx, ky, dop, T, t,t1,t2,t3,t4,aa,ww)
    CC=CT*1.0/NF
    CC1=CT*1.0/(NN**2)
    return (CC,CC1)


def SUS(kx,ky,T,t,t1,t2,t3,t4):
    EF=tb(kx,ky,t,t1,t2,t3,t4)
    dop=2*(0.5-dd(EF,t,t1,t2,t3,t4)+0.1)
    ww=AMAX(dop)[0]
    aa=AMAX(dop)[1]
    gg=GAMMA(kx, ky, ww, T, t,t1,t2,t3,t4,aa,ww)
    return gg

def TC(EF,t,t1,t2,t3,t4):
    tc=0
    for tt in range(NN):
       TT=0.5-0.5*tt/NN
       if abs(CCT(EF,TT,t,t1,t2,t3,t4))<1:
         tc=TT
         break
    return tc
for ix in range(NR+1):
  print ix
  for iy in range(NR+1):
    kx=2*pi*ix/NR
    ky=2*pi*iy/NR
    EF=tb(kx,ky,-0.5908,0.0962,-0.1306,-0.0507,0.0939)
    dop=2*(0.5-dd(EF,-0.5908,0.0962,-0.1306,-0.0507,0.0939)+0.1)
    ww=AMAX(dop)[1]
    aa=AMAX(dop)[0]
    if tb(kx,ky,-0.5908,0.0962,-0.1306,-0.0507,0.0939)-tb(kx+pi,ky+pi,-0.5908,0.0962,-0.1306,-0.0507,0.0939)<ww:
     gg=SUS(kx,ky,0.005,-0.5908,0.0962,-0.1306,-0.0507,0.0939)
     if gg<0:
       ss=gg
     else: ss=0
    else: ss=0
    ff2.write(str(ss)+'\t')
  ff2.write('\n')
