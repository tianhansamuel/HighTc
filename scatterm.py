from numpy import *
from math import *


NN=100
N1=100
NB=100
ff1=open('scatter0_0.03_re.dat','w')
ff2=open('scatter0_0.03_imag.dat','w')

TEMP=0.02

eta1=0.01

eta=0.05
dd=1

ed=0.03


def tb(px,py,t,t1,t2,t3,t4):
    EE=0.5*t*(cos(px)+cos(py))+t1*cos(px)*cos(py)+0.5*t2*(cos(2*px)+cos(2*py))+0.5*t3*(cos(2*px)*cos(py)+cos(px)*cos(2*py))+t4*cos(2*px)*cos(2*py)
    return EE

def GG(omega,px,py,t,t1,t2,t3,t4,aa,ww,dop):
   GG=1.0/(omega-tb(px,py,t,t1,t2,t3,t4)-Self(omega,px,py,t,t1,t2,t3,t4,aa,ww,dop)+eta1*1j) 
   return GG


e1=tb(pi,0.25*pi,-0.5908,0.0962,-0.1306,-0.0507,0.0939)

def GG0(omega,px,py,t,t1,t2,t3,t4):
   GG=1.0/(omega-tb(px,py,t,t1,t2,t3,t4)+eta*1j) 
   return GG

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

def Self(omega,px,py,t,t1,t2,t3,t4,aa,ww,dop):
  SS=0
  NS=1
  if tb(px,py,t,t1,t2,t3,t4)-tb(px+pi,py+pi,t,t1,t2,t3,t4)>ww:
   for ix in range(N1):
    for iy in range(N1):
     qx=2*ix*pi/N1
     qy=2*iy*pi/N1
     wwq=ww/4*(6+cos(qx)+cos(qy))
     if tb(px,py,t,t1,t2,t3,t4)-tb(px+pi+qx,py+pi+qy,t,t1,t2,t3,t4)>ww and tb(px+qx,py+qy,t,t1,t2,t3,t4)<tb(px,py,t,t1,t2,t3,t4):
       SS+=pi/2*(1-12*aa**2)/(1+12*aa**2)*(1-dop)**2*(tb(px+qx,py+qy,t,t1,t2,t3,t4))**2/(omega-tb(px+pi+qx,py+pi+qy,t,t1,t2,t3,t4)-wwq+eta*1j)
       NS+=1
  SS=SS/NS
  return SS

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


def AA(omega,kx,ky,t,t1,t2,t3,t4,VS,aa,ww,dop):
   SS=0
   SS1=0
   for ix in range(NN):
     for iy in range(NN):
        qx=-pi+2*pi*ix/NN
        qy=-pi+2*pi*iy/NN
        SS+=1.0/(NN**2)*GG(omega,qx,qy,t,t1,t2,t3,t4,aa,ww,dop)*GG(omega,kx+qx,ky+qy,t,t1,t2,t3,t4,aa,ww,dop)
        SS1+=1.0/(NN**2)*GG(omega,qx,qy,t,t1,t2,t3,t4,aa,ww,dop)
   return SS*(VS-SS1)


for ix in range(NN):
 print ix
 for iy in range(NN):
   kx=2*pi*ix/NN
   ky=2*pi*iy/NN
   if tb(kx,ky,-0.5908,0.0962,-0.1306,-0.0507,0.0939)-tb(kx+pi,ky+pi,-0.5908,0.0962,-0.1306,-0.0507,0.0939)>0.04:
    TT=GG(e1+ed,kx,ky,-0.5908,0.0962,-0.1306,-0.0507,0.0939,0.07,0.04,0.14).real
    TT1=GG(e1+ed,kx,ky,-0.5908,0.0962,-0.1306,-0.0507,0.0939,0.07,0.04,0.14).imag
   else:
    TT=GG(e1+ed,kx,ky,-0.5908,0.0962,-0.1306,-0.0507,0.0939,0.07,0.04,0.14).real
    TT1=GG(e1+ed,kx,ky,-0.5908,0.0962,-0.1306,-0.0507,0.0939,0.07,0.04,0.14).imag
   ff1.write(str(TT)+'\t')
   ff2.write(str(TT1)+'\t')
 ff1.write('\n')
 ff2.write('\n')


### 0. 0.25pi  2. 0.10pi
