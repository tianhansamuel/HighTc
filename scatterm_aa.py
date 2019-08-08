from numpy import *
from math import *


NN=100
N1=100
NB=100
ff1=open('cc_re.dat','w')
ff4=open('cc_imag.dat','w')
ff2=open('scattercc_imag.dat','w')
ff3=open('scattercc_real.dat','w')

TEMP=0.02

eta1=0.01

eta=0.1
dd=1

ed=0


def tb(px,py,t,t1,t2,t3,t4):
    EE=0.5*t*(cos(px)+cos(py))+t1*cos(px)*cos(py)+0.5*t2*(cos(2*px)+cos(2*py))+0.5*t3*(cos(2*px)*cos(py)+cos(px)*cos(2*py))+t4*cos(2*px)*cos(2*py)
    return EE

def GG(omega,px,py,t,t1,t2,t3,t4,aa,ww,dop):
   GG=1.0/(omega-tb(px,py,t,t1,t2,t3,t4)+eta1*1j) 
   return GG


e1=tb(pi,0.1*pi,-0.5908,0.0962,-0.1306,-0.0507,0.0939)

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
   for jx1 in range(NN):
     for jy1 in range(NN):
        qx=-pi+2*pi*jx1/NN
        qy=-pi+2*pi*jy1/NN
        SS+=1.0/(NN**2)*GG(omega,qx,qy,t,t1,t2,t3,t4,aa,ww,dop)*GG(omega,kx+qx,ky+qy,t,t1,t2,t3,t4,aa,ww,dop)
        SS1+=1.0/(2*pi*NN)**2*GG(omega,qx,qy,t,t1,t2,t3,t4,aa,ww,dop)
   return SS*1.0/(1.0/VS-SS1)

ll=1
for ix in range(NN+1):
 print ix
 for iy in range(NN+1):
   kx=-pi+2*pi*ix*1.0/NN
   ky=-pi+2*pi*iy*1.0/NN
   CC=AA(e1+ed,kx,ky,-0.5908,0.0962,-0.1306,-0.0507,0.0939,0.1,0.07,0.04,0.14)
   T=GG(e1+ed,kx,ky,-0.5908,0.0962,-0.1306,-0.0507,0.0939,0.07,0.04,0.14)
   TT=T.real
   TT0=T.imag
   TT1=CC.imag
   TT2=CC.real
   ff1.write(str(TT)+'\t')
   ff4.write(str(TT0)+'\t')

   ff2.write(str(TT1)+'\t')
   ff3.write(str(TT2)+'\t')
 ff1.write('\n')
 ff2.write('\n')
 ff3.write('\n')
 ff4.write('\n')
### 0. 0.25pi  2. 0.17pi 3. 0.15pi
