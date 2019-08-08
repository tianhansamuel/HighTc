from numpy import *
from math import *


NN=300
N1=300
NB=300
ff1=open('test1.dat','w')

dd=1

def tb(px,py,t,t1,t2,t3,t4):
    EE=0.5*t*(cos(px)+cos(py))+t1*cos(px)*cos(py)+0.5*t2*(cos(2*px)+cos(2*py))+0.5*t3*(cos(2*px)*cos(py)+cos(px)*cos(2*py))+t4*cos(2*px)*cos(2*py)
    return EE

def GG(omega,px,py,t,t1,t2,t3,t4):
   GG=1.0/(omega-tb(px,py,t,t1,t2,t3,t4)-Self(omega,px,py,t,t1,t2,t3,t4)+0.01*1j) 
   return GG


def GB(omega,px,py,t,t1,t2,t3,t4,DD):
    GBB=(omega+tb(px,py,t,t1,t2,t3,t4)+0.01*1j)/((omega+0.01*1j)**2-tb(px,py,t,t1,t2,t3,t4)**2-(DD*(cos(px)-cos(py)))**2)
    return GBB


def gn(omega2,EF,t,t1,t2,t3,t4,DD):
  GR=0
  NFG=0
  for ix in range(NN):
   for iy in range(NN):
     kx=-pi+2*pi*ix/NN
     ky=-pi+2*pi*iy/NN    
     if abs(tb(kx,ky,t,t1,t2,t3,t4)-EF)<0.01:
       NFG+=1
       GR+=GB(omega2,kx,ky,t,t1,t2,t3,t4,DD)
  GR=GR/NFG
  return GR


for ii in range(NN):
    print ii
    omega1=-0.4+0.8*ii/NN
    G1=gn(omega1-0.04,-0.04,-0.5908,0.0962,-0.1306,-0.0507,0.0939,0.1).imag
    G2=gn(omega1-0.06,-0.06,-0.5908,0.0962,-0.1306,-0.0507,0.0939,0.1).imag
    G3=gn(omega1-0.08,-0.08,-0.5908,0.0962,-0.1306,-0.0507,0.0939,0.1).imag
    G4=gn(omega1-0.10,-0.10,-0.5908,0.0962,-0.1306,-0.0507,0.0939,0.1).imag
    G8=gn(omega1-0.11,-0.11,-0.5908,0.0962,-0.1306,-0.0507,0.0939,0.1).imag
    G5=gn(omega1-0.12,-0.12,-0.5908,0.0962,-0.1306,-0.0507,0.0939,0.1).imag
    G6=gn(omega1-0.14,-0.14,-0.5908,0.0962,-0.1306,-0.0507,0.0939,0.1).imag
    G7=gn(omega1-0.16,-0.16,-0.5908,0.0962,-0.1306,-0.0507,0.0939,0.1).imag
    ff1.write(str(omega1)+'\t')
    ff1.write(str(G1)+'\t')
    ff1.write(str(G2)+'\t')
    ff1.write(str(G3)+'\t')
    ff1.write(str(G4)+'\t')
    ff1.write(str(G8)+'\t')
    ff1.write(str(G5)+'\t')
    ff1.write(str(G6)+'\t')
    ff1.write(str(G7)+'\n')


