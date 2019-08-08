from numpy import *
from math import *


NN=100
N1=100
NB=100


sscc=loadtxt('scatterm_4_re.dat', dtype =float)  
ssii=loadtxt('scatterm_4_imag.dat', dtype =float) 

SS1=shape(sscc)

eta1=0.01

eta=0.05
dd=1

ff1=open('scatterm4_0.dat','w')


MM=zeros((2*NN,2*NN),'complex')

for ix in range(NN):
  for iy in range(NN):
   for ix1 in range(NN):
     for iy1 in range(NN):
        iix=int((ix+ix1-NN/2)%NN)
        iiy=int((iy+iy1-NN/2)%NN)
        SS+=1.0/(NN**2)*(sscc[ix][iy]+1j*ssii[ix][iy])*(sscc[iix][iiy]+1j*ssii[iix][iiy])
        SS1+=1.0/(2*pi*NN)**2*(sscc[ix][iy]+1j*ssii[ix][iy])
   return SS*1.0/(1.0/VS-SS1).imag

for ix in range(NN+1):
  print ix
  for iy in range(NN+1):
      Am=AA(ix,iy,0.1)
      ff1.write(str(Am.imag)+'\t')
  ff1.write('\n')

