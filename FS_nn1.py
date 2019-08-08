from numpy import *
from math import *


NN=300
N1=200
NB=200
ff1=open('specT_h1.dat','w')
ff2=open('specT_h2.dat','w')
ff3=open('specT_h3.dat','w')
ff4=open('specT_h4.dat','w')
ff5=open('specT_h5.dat','w')
ff6=open('specT_h6.dat','w')
ff7=open('specT_h7.dat','w')
ff8=open('specT_h8.dat','w')
ff9=open('specT_h9.dat','w')

dd=1

eta=0.05
eta1=0.01
TEMP=0.02

def tb(px,py,t,t1,t2,t3,t4):
    EE=0.5*t*(cos(px)+cos(py))+t1*cos(px)*cos(py)+0.5*t2*(cos(2*px)+cos(2*py))+0.5*t3*(cos(2*px)*cos(py)+cos(px)*cos(2*py))+t4*cos(2*px)*cos(2*py)
    return EE

def GG(omega,px,py,t,t1,t2,t3,t4):
   GG=1.0/(omega-tb(px,py,t,t1,t2,t3,t4)-Self(omega,px,py,t,t1,t2,t3,t4)+eta1*1j) 
   return GG


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

def Self(omega,px,py,t,t1,t2,t3,t4):
  dop=2*(0.5-dd(tb(px,py,t,t1,t2,t3,t4),t,t1,t2,t3,t4)+0.1)
  aa=AMAX(dop)[1]
  ww=AMAX(dop)[0]
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


def GTT(omega,px,py,t,t1,t2,t3,t4,T):
   Gt=0
   for ii in range(N1):
    ot=(2*ii+1)*T
    Gt+=GG(omega+ot,px,py,t,t1,t2,t3,t4)/N1
   return GG


for ii in range(NN):
   print "1,"
   print ii
   omega=-0.5*ii/NN+0.2
   ee1=tb(1.57393791945,0.911061869541,-0.5908,0.0962,-0.1306,-0.0507,0.0939)
   TT=GG(omega+ee1,1.57393791945,0.911061869541,-0.5908,0.0962,-0.1306,-0.0507,0.0939)*1.0/(exp(1.0*(omega+ee1)/TEMP)+1)
   ff1.write(str(omega)+'\t')
   ff1.write(str(TT.imag)+'\t')
   ff1.write(str(TT.real)+'\n')


for ii in range(NN):
   print "2,"
   print ii
   omega=-0.5*ii/NN+0.2
   ee1=tb(1.70274321825,0.860796387084,-0.5908,0.0962,-0.1306,-0.0507,0.0939)
   TT=GG(omega+ee1,1.70274321825,0.860796387084,-0.5908,0.0962,-0.1306,-0.0507,0.0939)*1.0/(exp(1.0*(omega+ee1)/TEMP)+1)
   ff2.write(str(omega)+'\t')
   ff2.write(str(TT.imag)+'\t')
   ff2.write(str(TT.real)+'\n')

for ii in range(NN):
   print "3,"
   print ii
   omega=-0.5*ii/NN+0.2
   ee1=tb(2.01533168728,0.782256570744,-0.5908,0.0962,-0.1306,-0.0507,0.0939)
   TT=GG(omega+ee1,2.01533168728,0.782256570744,-0.5908,0.0962,-0.1306,-0.0507,0.0939)*1.0/(exp(1.0*(omega+ee1)/TEMP)+1)
   ff3.write(str(omega)+'\t')
   ff3.write(str(TT.imag)+'\t')
   ff3.write(str(TT.real)+'\n')

for ii in range(NN):
   print "4,"
   print ii
   omega=-0.5*ii/NN+0.2
   ee1=tb(2.35619449019,0.750840644208,-0.5908,0.0962,-0.1306,-0.0507,0.0939)
   TT=GG(omega+ee1,2.35619449019,0.750840644208,-0.5908,0.0962,-0.1306,-0.0507,0.0939)*1.0/(exp(1.0*(omega+ee1)/TEMP)+1)
   ff4.write(str(omega)+'\t')
   ff4.write(str(TT.imag)+'\t')
   ff4.write(str(TT.real)+'\n')



for ii in range(NN):
   print "5,"
   print ii
   omega=-0.5*ii/NN+0.2
   ee1=tb(2.561968809,0.757123829515,-0.5908,0.0962,-0.1306,-0.0507,0.0939)
   TT=GG(omega+ee1,2.561968809,0.757123829515,-0.5908,0.0962,-0.1306,-0.0507,0.0939)*1.0/(exp(1.0*(omega+ee1)/TEMP)+1)
   ff5.write(str(omega)+'\t')
   ff5.write(str(TT.imag)+'\t')
   ff5.write(str(TT.real)+'\n')



for ii in range(NN):
   print "6,"
   print ii
   omega=-0.5*ii/NN+0.2
   ee1=tb(2.70334047841,0.769690200129,-0.5908,0.0962,-0.1306,-0.0507,0.0939)
   TT=GG(omega+ee1,2.70334047841,0.769690200129,-0.5908,0.0962,-0.1306,-0.0507,0.0939)*1.0/(exp(1.0*(omega+ee1)/TEMP)+1)
   ff6.write(str(omega)+'\t')
   ff6.write(str(TT.imag)+'\t')
   ff6.write(str(TT.real)+'\n')

for ii in range(NN):
   print "7,"
   print ii
   omega=-0.5*ii/NN+0.2
   ee1=tb(2.90911479722,0.794822941358,-0.5908,0.0962,-0.1306,-0.0507,0.0939)
   TT=GG(omega+ee1,2.90911479722,0.794822941358,-0.5908,0.0962,-0.1306,-0.0507,0.0939)*1.0/(exp(1.0*(omega+ee1)/TEMP)+1)
   ff7.write(str(omega)+'\t')
   ff7.write(str(TT.imag)+'\t')
   ff7.write(str(TT.real)+'\n')

for ii in range(NN):
   print "8,"
   print ii
   omega=-0.5*ii/NN+0.2
   ee1=tb(3.00022098418,0.804247719319,-0.5908,0.0962,-0.1306,-0.0507,0.0939)
   TT=GG(omega+ee1,3.00022098418,0.804247719319,-0.5908,0.0962,-0.1306,-0.0507,0.0939)*1.0/(exp(1.0*(omega+ee1)/TEMP)+1)
   ff8.write(str(omega)+'\t')
   ff8.write(str(TT.imag)+'\t')
   ff8.write(str(TT.real)+'\n')

for ii in range(NN):
   print "9,"
   print ii
   omega=-0.5*ii/NN+0.2
   ee1=tb(3.14002185726,0.810530904626,-0.5908,0.0962,-0.1306,-0.0507,0.0939)
   TT=GG(omega+ee1,3.14002185726,0.810530904626,-0.5908,0.0962,-0.1306,-0.0507,0.0939)*1.0/(exp(1.0*(omega+ee1)/TEMP)+1)
   ff9.write(str(omega)+'\t')
   ff9.write(str(TT.imag)+'\t')
   ff9.write(str(TT.real)+'\n')

