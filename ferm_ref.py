from numpy import *
from math import *


NN=1000
N1=1000
NB=1000

def tb(px,py,t,t1,t2,t3,t4):
    EE=0.5*t*(cos(px)+cos(py))+t1*cos(px)*cos(py)+0.5*t2*(cos(2*px)+cos(2*py))+0.5*t3*(cos(2*px)*cos(py)+cos(px)*cos(2*py))+t4*cos(2*px)*cos(2*py)
    return EE

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


e1=tb(pi,0.1*pi,-0.5908,0.0962,-0.1306,-0.0507,0.0939)
dope1=(0.5-density(e1,-0.5908,0.0962,-0.1306,-0.0507,0.0939)+0.1)

e2=tb(pi,0.2*pi,-0.5908,0.0962,-0.1306,-0.0507,0.0939)
dope2=(0.5-density(e2,-0.5908,0.0962,-0.1306,-0.0507,0.0939)+0.1)

e3=tb(pi,0.25*pi,-0.5908,0.0962,-0.1306,-0.0507,0.0939)
dope3=(0.5-density(e3,-0.5908,0.0962,-0.1306,-0.0507,0.0939)+0.1)

e4=tb(pi,0.28*pi,-0.5908,0.0962,-0.1306,-0.0507,0.0939)
dope4=(0.5-density(e4,-0.5908,0.0962,-0.1306,-0.0507,0.0939)+0.1)

e5=tb(pi,0.3*pi,-0.5908,0.0962,-0.1306,-0.0507,0.0939)
dope5=(0.5-density(e5,-0.5908,0.0962,-0.1306,-0.0507,0.0939)+0.1)

e6=tb(pi,0.32*pi,-0.5908,0.0962,-0.1306,-0.0507,0.0939)
dope6=(0.5-density(e6,-0.5908,0.0962,-0.1306,-0.0507,0.0939)+0.1)

e7=tb(pi,0.35*pi,-0.5908,0.0962,-0.1306,-0.0507,0.0939)
dope7=(0.5-density(e7,-0.5908,0.0962,-0.1306,-0.0507,0.0939)+0.1)

print dope1
print dope2
print dope3
print dope4
print dope5
print dope6
print dope7

### 1. 0.1pi 3mv 2. 0.2pi 3mv 3. 0.25pi 5mv 4. 0.28pi 5mv 5.0.3pi 5mv 6.0.32pi 9mv 7. 0.35pi 9mv
