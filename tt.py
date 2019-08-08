from numpy import *
from math import *


NN=100
N1=100
NB=100
eta1=0.01

sscc=loadtxt('sc3_05re.dat', dtype =float)  
ssii=loadtxt('sc3_05imag.dat', dtype =float)   

sscc1=loadtxt('scn3_05re.dat', dtype =float)  
ssii1=loadtxt('scn3_05imag.dat', dtype =float)   


SS1=shape(sscc)
VS=0.1

ff1=open('ssc_305.dat','w')

ff2=open('impsc_305.dat','w')

def tb(px,py,t,t1,t2,t3,t4):
    EE=0.5*t*(cos(px)+cos(py))+t1*cos(px)*cos(py)+0.5*t2*(cos(2*px)+cos(2*py))+0.5*t3*(cos(2*px)*cos(py)+cos(px)*cos(2*py))+t4*cos(2*px)*cos(2*py)
    return EE

e1=tb(pi,0.25*pi,-0.5908,0.0962,-0.1306,-0.0507,0.0939)


SS1=0
for jx1 in range(NN):
     for jy1 in range(NN):
        SS1+=(sscc[jx1][jy1]+1j*ssii[jx1][jy1])*1.0/(2*pi*NN)**2



TV=zeros((2,2),'complex')
TV[0][0]=1.0/VS
TV[1][1]=-1.0/VS
TV-=SS1
TT=linalg.inv(TV)

def AA(ix1,iy1):
   SS=0
   SS1=0
   for ix in range(NN+1):
     for iy in range(NN+1):
        iix=int((ix+ix1-0.5*NN)%NN)
        iiy=int((iy+iy1-0.5*NN)%NN)
        G1=zeros((2,2),'complex')
        G2=zeros((2,2),'complex')
        G1[0][0]=sscc[ix][4*iy]+1j*ssii[ix][4*iy]
        G1[1][0]=sscc[ix][4*iy+1]+1j*ssii[ix][4*iy+1]
        G1[0][1]=sscc[ix][4*iy+2]+1j*ssii[ix][4*iy+2]
        G1[1][1]=sscc[ix][4*iy+3]+1j*ssii[ix][4*iy+3]

        G2[0][0]=sscc1[iix][4*iiy]+1j*ssii1[iix][4*iiy]
        G2[1][0]=sscc1[iix][4*iiy+1]+1j*ssii1[iix][4*iiy+1]
        G2[0][1]=sscc1[iix][4*iiy+2]+1j*ssii1[iix][4*iiy+2]
        G2[1][1]=sscc1[iix][4*iiy+3]+1j*ssii1[iix][4*iiy+3]
        M1=matmul(G1, TT)
        SS+=1.0/(NN**2)*matmul(M1, G2)
   AS=SS[0][0].imag+SS[1][1].imag
   return AS

ll=10

for ix in range(NN+1):
  print ix
  for iy in range(NN+1):
      Am=ssii[ix%NN][4*(iy%NN)]+ssii1[ix%NN][4*(iy%NN)+3]
      Bm=AA(ix%NN,iy%NN)
      ff1.write(str(Am)+'\t')
      ff2.write(str(Bm)+'\t')
  ff1.write('\n')
  ff2.write('\n')

