from numpy import *
from math import *


NN=50
N1=50
NB=300
ff1=open('test_single_2.dat','w')

dd=1

def tb(px,py,t,t1,t2,t3,t4):
    EE=0.5*t*(cos(px)+cos(py))+t1*cos(px)*cos(py)+0.5*t2*(cos(2*px)+cos(2*py))+0.5*t3*(cos(2*px)*cos(py)+cos(px)*cos(2*py))+t4*cos(2*px)*cos(2*py)
    return EE

def EB(px,py,t,t1,t2,t3,t4,DD):
    EE=sqrt(tb(px,py,t,t1,t2,t3,t4)**2+0.5*DD**2*(cos(px)-cos(py))**2)
    return EE

def GG(omega,px,py,t,t1,t2,t3,t4,aa,dop,mu):
   GG=1.0/(omega-tb(px,py,t,t1,t2,t3,t4)+mu-Self(omega,px,py,t,t1,t2,t3,t4,aa,dop,mu)+0.05*1j) 
   return GG

def Self(omega,px,py,t,t1,t2,t3,t4,aa,dop,mu):
  SS=0
  NS=1
  ww=4*aa/(1-12*aa**2)*0.129
  if tb(px,py,t,t1,t2,t3,t4)-tb(px+pi,py+pi,t,t1,t2,t3,t4)>ww:
   for ix in range(N1):
    for iy in range(N1):
     qx=2*ix*pi/N1
     qy=2*iy*pi/N1
     wwq=ww/4*(6+cos(qx)+cos(qy))
     if tb(px,py,t,t1,t2,t3,t4)-tb(px+pi+qx,py+pi+qy,t,t1,t2,t3,t4)>ww and tb(px+qx,py+qy,t,t1,t2,t3,t4)<=tb(px,py,t,t1,t2,t3,t4):
       SS+=pi*(1-12*aa**2)/(1+12*aa**2)*(1-dop)**2*(tb(px+qx,py+qy,t,t1,t2,t3,t4))**2/(omega-tb(px+pi+qx,py+pi+qy,t,t1,t2,t3,t4)+mu-wwq+0.05*1j)
       NS+=1
  SS=SS/NS
  return SS



def gn(omega2,EF,t,t1,t2,t3,t4,aa,dop,mu):
  GR=0
  NFG=0
  for ix in range(NN):
   for iy in range(NN):
     kx=-pi+2*pi*ix/NN
     ky=-pi+2*pi*iy/NN    
     if abs(tb(kx,ky,t,t1,t2,t3,t4)-EF)<0.01:
       NFG+=1
       GR+=-2*GG(omega2,kx,ky,t,t1,t2,t3,t4,aa,dop,mu)
  GR=GR/NFG
  return (GR,NFG)

for ii in range(NB):
    print ii
    omega1=-0.35+0.7*ii/NB
    G1=gn(omega1,-0.04,-0.5908,0.0962,-0.1306,-0.0507,0.0939,0.060,0.08,-0.04)
    G2=gn(omega1,-0.06,-0.5908,0.0962,-0.1306,-0.0507,0.0939,0.065,0.10,-0.06)
    G3=gn(omega1,-0.08,-0.5908,0.0962,-0.1306,-0.0507,0.0939,0.068,0.12,-0.08)
    G11=gn(omega1,-0.09,-0.5908,0.0962,-0.1306,-0.0507,0.0939,0.070,0.14,-0.09)
    G4=gn(omega1,-0.10,-0.5908,0.0962,-0.1306,-0.0507,0.0939,0.072,0.14,-0.10)
    G8=gn(omega1,-0.11,-0.5908,0.0962,-0.1306,-0.0507,0.0939,0.075,0.16,-0.11)
    G10=gn(omega1,-0.115,-0.5908,0.0962,-0.1306,-0.0507,0.0939,0.074,0.17,-0.115)
    G5=gn(omega1,-0.12,-0.5908,0.0962,-0.1306,-0.0507,0.0939,0.072,0.18,-0.12)
    G9=gn(omega1,-0.13,-0.5908,0.0962,-0.1306,-0.0507,0.0939,0.070,0.19,-0.13)
    G6=gn(omega1,-0.14,-0.5908,0.0962,-0.1306,-0.0507,0.0939,0.068,0.20,-0.14)
    G7=gn(omega1,-0.16,-0.5908,0.0962,-0.1306,-0.0507,0.0939,0.065,0.22,-0.16)

    GG1=G1[0].imag
    GG2=G2[0].imag
    GG3=G3[0].imag
    GG11=G11[0].imag
    GG4=G4[0].imag
    GG8=G8[0].imag
    GG10=G10[0].imag
    GG5=G5[0].imag
    GG9=G9[0].imag
    GG6=G6[0].imag
    GG7=G7[0].imag

    N1=G1[1]
    N2=G2[1]
    N3=G3[1]
    N11=G11[1]
    N4=G4[1]
    N8=G8[1]
    N10=G10[1]
    N5=G5[1]
    N9=G9[1]
    N6=G6[1]
    N7=G7[1]
    ff1.write(str(omega1)+'\t')
    ff1.write(str(GG1)+'\t')
    ff1.write(str(GG2)+'\t')
    ff1.write(str(GG3)+'\t')
    ff1.write(str(GG11)+'\t')
    ff1.write(str(GG4)+'\t')
    ff1.write(str(GG8)+'\t')
    ff1.write(str(GG10)+'\t')
    ff1.write(str(GG5)+'\t')
    ff1.write(str(GG9)+'\t')
    ff1.write(str(GG6)+'\t')
    ff1.write(str(GG7)+'\t')

    ff1.write(str(N1)+'\t')
    ff1.write(str(N2)+'\t')
    ff1.write(str(N3)+'\t')
    ff1.write(str(N11)+'\t')
    ff1.write(str(N4)+'\t')
    ff1.write(str(N8)+'\t')
    ff1.write(str(N10)+'\t')
    ff1.write(str(N5)+'\t')
    ff1.write(str(N9)+'\t')
    ff1.write(str(N6)+'\t')
    ff1.write(str(N7)+'\n')
