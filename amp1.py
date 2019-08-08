from numpy import *
from math import *


NN=100
N1=100
NB=100


sscc=loadtxt('scatter4_m_0_re.dat', dtype =float)
sscc1=loadtxt('scatter3_m_0_re.dat', dtype =float)   
sscc2=loadtxt('scatter2_m_0_re.dat', dtype =float)   
sscc3=loadtxt('scatter1_m_0_re.dat', dtype =float)   
sscc4=loadtxt('scatter0_m_0_re.dat', dtype =float)   
sscc5=loadtxt('scatter5_m_0_re.dat', dtype =float)

ssii=loadtxt('scatter4_m_0_imag.dat', dtype =float)
ssii1=loadtxt('scatter3_m_0_imag.dat', dtype =float)   
ssii2=loadtxt('scatter2_m_0_imag.dat', dtype =float)   
ssii3=loadtxt('scatter1_m_0_imag.dat', dtype =float)
ssii4=loadtxt('scatter0_m_0_imag.dat', dtype =float)   
ssii5=loadtxt('scatter5_m_0_imag.dat', dtype =float)


SS1=shape(sscc)

eta1=0.01

eta=0.05
dd=1

ff1=open('scam4_0.dat','w')
ff2=open('scam3_0.dat','w')
ff3=open('scam2_0.dat','w')
ff4=open('scam1_0.dat','w')
ff5=open('scam0_0.dat','w')
ff6=open('scam5_0.dat','w')


def AA(ix1,iy1,VS):
   SS=0
   SS1=0
   for ix in range(NN):
     for iy in range(NN):
        iix=int((ix+ix1)%NN)
        iiy=int((iy+iy1)%NN)
        SS+=1.0/(NN**2)*(sscc[ix][iy]+1j*ssii[ix][iy])*(sscc[iix][iiy]+1j*ssii[iix][iiy])
        SS1+=1.0/(2*pi*NN)**2*(sscc[ix][iy]+1j*ssii[ix][iy])
   return SS*1.0/(1.0/VS-SS1).imag


def AA1(ix1,iy1,VS):
   SS=0
   SS1=0
   for ix in range(NN):
     for iy in range(NN):
        iix=int((ix+ix1)%NN)
        iiy=int((iy+iy1)%NN)
        SS+=1.0/(NN**2)*(sscc1[ix][iy]+1j*ssii1[ix][iy])*(sscc1[iix][iiy]+1j*ssii1[iix][iiy])
        SS1+=1.0/(2*pi*NN)**2*(sscc1[ix][iy]+1j*ssii1[ix][iy])
   return SS*1.0/(1.0/VS-SS1).imag

def AA2(ix1,iy1,VS):
   SS=0
   SS1=0
   for ix in range(NN):
     for iy in range(NN):
        iix=int((ix+ix1)%NN)
        iiy=int((iy+iy1)%NN)
        SS+=1.0/(NN**2)*(sscc2[ix][iy]+1j*ssii2[ix][iy])*(sscc2[iix][iiy]+1j*ssii2[iix][iiy])
        SS1+=1.0/(2*pi*NN)**2*(sscc2[ix][iy]+1j*ssii2[ix][iy])
   return SS*1.0/(1.0/VS-SS1).imag

def AA3(ix1,iy1,VS):
   SS=0
   SS1=0
   for ix in range(NN):
     for iy in range(NN):
        iix=int((ix+ix1)%NN)
        iiy=int((iy+iy1)%NN)
        SS+=1.0/(NN**2)*(sscc3[ix][iy]+1j*ssii3[ix][iy])*(sscc3[iix][iiy]+1j*ssii3[iix][iiy])
        SS1+=1.0/(2*pi*NN)**2*(sscc3[ix][iy]+1j*ssii3[ix][iy])
   return SS*1.0/(1.0/VS-SS1).imag

def AA4(ix1,iy1,VS):
   SS=0
   SS1=0
   for ix in range(NN):
     for iy in range(NN):
        iix=int((ix+ix1)%NN)
        iiy=int((iy+iy1)%NN)
        SS+=1.0/(NN**2)*(sscc4[ix][iy]+1j*ssii4[ix][iy])*(sscc4[iix][iiy]+1j*ssii4[iix][iiy])
        SS1+=1.0/(2*pi*NN)**2*(sscc4[ix][iy]+1j*ssii4[ix][iy])
   return SS*1.0/(1.0/VS-SS1).imag


def AA5(ix1,iy1,VS):
   SS=0
   SS1=0
   for ix in range(NN):
     for iy in range(NN):
        iix=int((ix+ix1)%NN)
        iiy=int((iy+iy1)%NN)
        SS+=1.0/(NN**2)*(sscc5[ix][iy]+1j*ssii5[ix][iy])*(sscc5[iix][iiy]+1j*ssii5[iix][iiy])
        SS1+=1.0/(2*pi*NN)**2*(sscc5[ix][iy]+1j*ssii5[ix][iy])
   return SS*1.0/(1.0/VS-SS1).imag

for ix in range(NN+1):
  print ix
  for iy in range(NN+1):
      Am=AA(ix,iy,0.1)
      Am1=AA1(ix,iy,0.1)
      Am2=AA2(ix,iy,0.1)
      Am3=AA3(ix,iy,0.1)
      Am4=AA4(ix,iy,0.1)
      Am5=AA5(ix,iy,0.1)
      ff1.write(str(Am.imag)+'\t')
      ff2.write(str(Am1.imag)+'\t')
      ff3.write(str(Am2.imag)+'\t')
      ff4.write(str(Am3.imag)+'\t')
      ff5.write(str(Am4.imag)+'\t')
      ff6.write(str(Am5.imag)+'\t')
  ff1.write('\n')
  ff2.write('\n')
  ff3.write('\n')
  ff4.write('\n')
  ff5.write('\n')
  ff6.write('\n')

