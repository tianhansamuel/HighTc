from numpy import *
from math import *


NN=100
N1=100
NB=100
ff1=open('ref2_0.05re.dat','w')
ff2=open('ref2_0.05imag.dat','w')

ff3=open('refn2_0.05re.dat','w')
ff4=open('refn2_0.05imag.dat','w')

TEMP=0.02

eta1=0.01

eta=0.1
dd=1

ed=0.05


def tb(px,py,t,t1,t2,t3,t4):
    EE=0.5*t*(cos(px)+cos(py))+t1*cos(px)*cos(py)+0.5*t2*(cos(2*px)+cos(2*py))+0.5*t3*(cos(2*px)*cos(py)+cos(px)*cos(2*py))+t4*cos(2*px)*cos(2*py)
    return EE

def GB(omega,px,py,t,t1,t2,t3,t4,DD,aa,ww,dop):
   GM=zeros((2,2),'complex')
   if tb(px,py,t,t1,t2,t3,t4)-tb(px+pi,py+pi,t,t1,t2,t3,t4)>ww:
    GM[0][0]=(omega+e1-tb(px,py,t,t1,t2,t3,t4)-Self(omega+e1,px,py,t,t1,t2,t3,t4,aa,ww,dop)+eta1*1j)
    GM[1][1]=(omega-e1+tb(px,py,t,t1,t2,t3,t4)+Self(omega+e1,px,py,t,t1,t2,t3,t4,aa,ww,dop)+eta1*1j) 
   else:
    GM[0][0]=(omega+e1-tb(px,py,t,t1,t2,t3,t4)+eta1*1j)
    GM[1][1]=(omega-e1+tb(px,py,t,t1,t2,t3,t4)+eta1*1j) 
    GM[0][1]=DD
    GM[1][0]=DD 
   return linalg.inv(GM)

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

def density(EF,t,t1,t2,t3,t4):
      NF=0
      for ix in range(NB+1):
        for iy in range(NB+1):
         kx=-pi+2*pi*ix/NB
         ky=-pi+2*pi*iy/NB
         if tb(kx,ky,t,t1,t2,t3,t4)<=EF:
           NF+=1
      dl=NF*1.0/(NB+1)**2
      return dl   

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


e1=tb(pi,0.2*pi,-0.5908,0.0962,-0.1306,-0.0507,0.0939)
dope=2*(0.5-density(e1,-0.5908,0.0962,-0.1306,-0.0507,0.0939)+0.1)
ww1=AMAX(dope)[1]
aa1=AMAX(dope)[0]





for ix in range(NN+1):
 for iy in range(NN+1):
   print iy
   kx=-pi+2*pi*ix*1.0/NN
   ky=-pi+2*pi*iy*1.0/NN
   CC=GB(ed,kx,ky,-0.5908,0.0962,-0.1306,-0.0507,0.0939,0,aa1,ww1,dope)
   CC1=GB(-ed,kx,ky,-0.5908,0.0962,-0.1306,-0.0507,0.0939,0,aa1,ww1,dope)

   ff1.write(str(CC[0][0].real)+'\t')
   ff1.write(str(CC[0][1].real)+'\t')
   ff1.write(str(CC[1][0].real)+'\t')
   ff1.write(str(CC[1][1].real)+'\t')

   ff2.write(str(CC[0][0].imag)+'\t')
   ff2.write(str(CC[0][1].imag)+'\t')
   ff2.write(str(CC[1][0].imag)+'\t')
   ff2.write(str(CC[1][1].imag)+'\t')

   ff3.write(str(CC1[0][0].real)+'\t')
   ff3.write(str(CC1[0][1].real)+'\t')
   ff3.write(str(CC1[1][0].real)+'\t')
   ff3.write(str(CC1[1][1].real)+'\t')

   ff4.write(str(CC1[0][0].imag)+'\t')
   ff4.write(str(CC1[0][1].imag)+'\t')
   ff4.write(str(CC1[1][0].imag)+'\t')
   ff4.write(str(CC1[1][1].imag)+'\t')

 ff1.write('\n')
 ff2.write('\n')
 ff3.write('\n')
 ff4.write('\n')

### 1. 0.1pi 3mv 2. 0.2pi 3mv 3. 0.25pi 5mv 4. 0.28pi 5mv 5.0.3pi 5mv 6.0.32pi 9mv 7. 0.35pi 9mv
