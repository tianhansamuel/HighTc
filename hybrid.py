from numpy import *
from math import *


NN=100
N1=100
NB=100

ll=1
sscc=loadtxt('scatterm_1_re.dat', dtype =float)  
ssii=loadtxt('scatterm_1_imag.dat', dtype =float)   

sif=loadtxt('scatterm1_0.dat', dtype =float)   


SS1=shape(sscc)
eta1=0.01

eta=0.05
dd=1

ff1=open('scatter1_0.dat','w')


for ix in range(NN+1):
  iix=ix%NN
  for iy in range(NN+1):
      iiy=iy%NN
      Am=ssii[iix][iiy]+ll*sif[iix][iiy]
      ff1.write(str(Am)+'\t')
  ff1.write('\n')

