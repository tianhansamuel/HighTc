from numpy import *
from math import *


NN=300
N1=200
NB=200
ff1=open('spec_ll1.dat','w')
ff2=open('spec_ll2.dat','w')
ff3=open('spec_ll3.dat','w')
ff4=open('spec_ll4.dat','w')

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
       SS+=pi/2*(1-12*ww**2)/(1+12*ww**2)*(1-dop)**2*(tb(px+qx,py+qy,t,t1,t2,t3,t4))**2/(omega-tb(px+pi+qx,py+pi+qy,t,t1,t2,t3,t4)-wwq+0.05*1j)
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



def GTT(omega,px,py,t,t1,t2,t3,t4,T):
   Gt=0
   for ii in range(N1):
    ot=(2*ii+1)*T
    Gt+=GG(omega+ot,px,py,t,t1,t2,t3,t4)/N1
   return GG

TEMP=0.0005

for ii in range(NN):
   print "1,"
   print ii
   omega=-0.5*ii/NN+0.2
   ee1=tb(0.753982236862,2.54469004941,-0.5908,0.0962,-0.1306,-0.0507,0.0939)
   TT=GG(omega+ee1,0.753982236862,2.54469004941,-0.5908,0.0962,-0.1306,-0.0507,0.0939)*1.0/(exp(1.0*(omega+ee1)/TEMP)+1)
   ff1.write(str(omega)+'\t')
   ff1.write(str(TT.imag)+'\t')
   ff1.write(str(TT.real)+'\n')

for ii in range(NN):
   print "2,"
   print ii
   omega=-0.5*ii/NN+0.2
   ee1=tb(0.806342114421,3.06828882501,-0.5908,0.0962,-0.1306,-0.0507,0.0939)
   TT=GG(omega+ee1,0.806342114421,3.06828882501,-0.5908,0.0962,-0.1306,-0.0507,0.0939)*1.0/(exp(1.0*(omega+ee1)/TEMP)+1)
   ff2.write(str(omega)+'\t')
   ff2.write(str(TT.imag)+'\t')
   ff2.write(str(TT.real)+'\n')


for ii in range(NN):
   print "3,"
   print ii
   omega=-0.5*ii/NN+0.2
   ee1=tb(0.806342114421,3.13112067808,-0.5908,0.0962,-0.1306,-0.0507,0.0939)
   TT=GG(omega+ee1,0.806342114421,3.13112067808,-0.5908,0.0962,-0.1306,-0.0507,0.0939)*1.0/(exp(1.0*(omega+ee1)/TEMP)+1)
   ff3.write(str(omega)+'\t')
   ff3.write(str(TT.imag)+'\t')
   ff3.write(str(TT.real)+'\n')

for ii in range(NN):
   print "4,"
   print ii
   omega=-0.5*ii/NN+0.2
   ee1=tb(0.774926187885,2.74365758414,-0.5908,0.0962,-0.1306,-0.0507,0.0939)
   TT=GG(omega+ee1,0.774926187885,2.74365758414,-0.5908,0.0962,-0.1306,-0.0507,0.0939)*1.0/(exp(1.0*(omega+ee1)/TEMP)+1)
   ff4.write(str(omega)+'\t')
   ff4.write(str(TT.imag)+'\t')
   ff4.write(str(TT.real)+'\n')
