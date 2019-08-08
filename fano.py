from numpy import *
from math import *


NN=100
N1=100
NB=100
ff1=open('spech.dat','w')

dd=1

def tb(px,py,t,t1,t2,t3,t4):
    EE=0.5*t*(cos(px)+cos(py))+t1*cos(px)*cos(py)+0.5*t2*(cos(2*px)+cos(2*py))+0.5*t3*(cos(2*px)*cos(py)+cos(px)*cos(2*py))+t4*cos(2*px)*cos(2*py)
    return EE

def GG(omega,px,py,t,t1,t2,t3,t4):
   GG=1.0/(omega-tb(px,py,t,t1,t2,t3,t4)-Self(omega,px,py,t,t1,t2,t3,t4)+0.01*1j) 
   return GG


def GG0(omega,px,py,t,t1,t2,t3,t4):
   GG=1.0/(omega-tb(px,py,t,t1,t2,t3,t4)+0.05*1j) 
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

def Self(omega,px,py,t,t1,t2,t3,t4):
  dop=2*(0.5-dd(tb(px,py,t,t1,t2,t3,t4),t,t1,t2,t3,t4)+0.1)
  ww=AMAX(dop)[1]
  SS=0
  NS=1
  if tb(px,py,t,t1,t2,t3,t4)-tb(px+pi,py+pi,t,t1,t2,t3,t4)>ww:
   for ix in range(N1):
    for iy in range(N1):
     qx=2*ix*pi/N1
     qy=2*iy*pi/N1
     wwq=ww/4*(6+cos(qx)+cos(qy))
     if tb(px,py,t,t1,t2,t3,t4)-tb(px+pi+qx,py+pi+qy,t,t1,t2,t3,t4)>ww and tb(px+qx,py+qy,t,t1,t2,t3,t4)<tb(px,py,t,t1,t2,t3,t4):
       SS+=pi/2*(tb(px+qx,py+qy,t,t1,t2,t3,t4))**2/(omega-tb(px+pi+qx,py+pi+qy,t,t1,t2,t3,t4)-wwq+0.05*1j)
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


def gn(omega2,EF,t,t1,t2,t3,t4):
  GR=0
  NFG=0
  for ix in range(NN):
   for iy in range(NN):
     kx=-pi+2*pi*ix/NN
     ky=-pi+2*pi*iy/NN    
     if abs(tb(kx,ky,t,t1,t2,t3,t4)-EF)<0.02:
       NFG+=1
       GR+=GG(omega2,kx,ky,t,t1,t2,t3,t4)
  GR=GR/NFG
  return GR


for ii in range(NN):
    print ii
    omega1=-0.35+0.7*ii/NN
    G1=gn(omega1-0.04,-0.04,-0.5908,0.0962,-0.1306,-0.0507,0.0939).imag
    G2=gn(omega1-0.06,-0.06,-0.5908,0.0962,-0.1306,-0.0507,0.0939).imag
    G3=gn(omega1-0.08,-0.08,-0.5908,0.0962,-0.1306,-0.0507,0.0939).imag
    G4=gn(omega1-0.1,-0.1,-0.5908,0.0962,-0.1306,-0.0507,0.0939).imag
    G8=gn(omega1-0.11,-0.11,-0.5908,0.0962,-0.1306,-0.0507,0.0939).imag
    G5=gn(omega1-0.12,-0.12,-0.5908,0.0962,-0.1306,-0.0507,0.0939).imag
    G6=gn(omega1-0.14,-0.14,-0.5908,0.0962,-0.1306,-0.0507,0.0939).imag
    G7=gn(omega1-0.16,-0.16,-0.5908,0.0962,-0.1306,-0.0507,0.0939).imag
    ff1.write(str(omega1)+'\t')
    ff1.write(str(G1)+'\t')
    ff1.write(str(G2)+'\t')
    ff1.write(str(G3)+'\t')
    ff1.write(str(G4)+'\t')
    ff1.write(str(G8)+'\t')
    ff1.write(str(G5)+'\t')
    ff1.write(str(G6)+'\t')
    ff1.write(str(G7)+'\n')


