from numpy import *
from math import *


NN=100
N1=100
NB=100
ff1=open('sc_cc_re.dat','w')
ff4=open('sc_cc_imag.dat','w')
ff2=open('scatter_cc_imag.dat','w')
ff3=open('scatter_cc_real.dat','w')

TEMP=0.02

eta1=0.01

eta=0.1
dd=1

ed=0


def tb(px,py,t,t1,t2,t3,t4):
    EE=0.5*t*(cos(px)+cos(py))+t1*cos(px)*cos(py)+0.5*t2*(cos(2*px)+cos(2*py))+0.5*t3*(cos(2*px)*cos(py)+cos(px)*cos(2*py))+t4*cos(2*px)*cos(2*py)
    return EE

e1=tb(pi,0.1*pi,-0.5908,0.0962,-0.1306,-0.0507,0.0939)

def GB(omega,px,py,t,t1,t2,t3,t4,DD):
   GM=zeros((2,2),'complex')
   GM[0][0]=(omega+e1-tb(px,py,t,t1,t2,t3,t4)+eta1*1j)
   GM[1][1]=(omega-e1+tb(px,py,t,t1,t2,t3,t4)+eta1*1j) 
   GM[0][1]=DD
   GM[1][0]=DD 
   return linalg.inv(GM)


def G0(omega,t,t1,t2,t3,t4,DD):
   SS1=0
   for jx1 in range(NN):
     for jy1 in range(NN):
        px=-pi+2*pi*jx1/NN
        py=-pi+2*pi*jy1/NN
        SS1+=GB(omega,px,py,t,t1,t2,t3,t4,DD)*1.0/(2*pi*NN)**2
   return SS1


def TT(omega,t,t1,t2,t3,t4,DD,VS):
   TV=zeros((2,2),'complex')
   TV[0][0]=1.0/VS
   TV[1][1]=-1.0/VS
   TV-=G0(omega,t,t1,t2,t3,t4,DD)
   return linalg.inv(TV)



def AA(omega,kx,ky,t,t1,t2,t3,t4,DD,VS):
   SS1=0
   SS2=0
   for jx1 in range(NN):
     for jy1 in range(NN):
        px=-pi+2*pi*jx1/NN
        py=-pi+2*pi*jy1/NN
        A1=GB(omega,px,py,t,t1,t2,t3,t4,DD)
        A2=GB(omega,px+kx,py+ky,t,t1,t2,t3,t4,DD)
        MM=matmul(A1, TT(omega,t,t1,t2,t3,t4,DD,VS))
        M1=matmul(MM,A2)

        B1=GB(omega,px,py,t,t1,t2,t3,t4,DD)
        B2=GB(omega,px+kx,py+ky,t,t1,t2,t3,t4,DD)
        MMT=matmul(B1, TT(-omega,t,t1,t2,t3,t4,DD,VS))
        MT1=matmul(MMT,B2)

        SS1+=1.0/(2*pi*NN)**2*M1
        SS2+=1.0/(2*pi*NN)**2*MT1
   return SS1[0][0]+SS2[1][1]
print GB(ed,pi,0.1*pi,-0.5908,0.0962,-0.1306,-0.0507,0.0939,0.03)
ll=1
for ix in range(NN+1):
 for iy in range(NN+1):
   print iy
   kx=-pi+2*pi*ix*1.0/NN
   ky=-pi+2*pi*iy*1.0/NN
   CC=(GB(ed,kx,ky,-0.5908,0.0962,-0.1306,-0.0507,0.0939,0.03)[0][0].imag+GB(-ed,kx,ky,-0.5908,0.0962,-0.1306,-0.0507,0.0939,0.03)[1][1].imag)
   CC1=AA(ed,kx,ky,-0.5908,0.0962,-0.1306,-0.0507,0.0939,0.03,0.1)[0][0].imag
   ff1.write(str(CC)+'\t')
   ff2.write(str(CC1)+'\t')
 ff1.write('\n')
 ff2.write('\n')

### 0. 0.25pi  2. 0.17pi 3. 0.15pi
