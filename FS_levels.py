from numpy import *
from math import *

NM=1000
NN=30
NB=500
DL=5.0
ff1=open('FS_10.dat','w')
ff2=open('FS_14.dat','w')
ff3=open('FS_17.dat','w')
ff4=open('FS_19.dat','w')
ff5=open('FS_20.dat','w')
ff6=open('FS_22.dat','w')
ff7=open('FS_24.dat','w')

def tb(px,py,t,t1,t2,t3,t4):
    EE=0.5*t*(cos(px)+cos(py))+t1*cos(px)*cos(py)+0.5*t2*(cos(2*px)+cos(2*py))+0.5*t3*(cos(2*px)*cos(py)+cos(px)*cos(2*py))+t4*cos(2*px)*cos(2*py)
    return EE

for ix in range(NB+1):
  for iy in range(NB+1):
    px=-pi+2*pi*ix/NB
    py=-pi+2*pi*iy/NB    
    ee=tb(px,py,-0.5908,0.0962,-0.1306,-0.0507,0.0939)
    if abs(ee+0.0362)<0.003 and tb(px,py,-0.5908,0.0962,-0.1306,-0.0507,0.0939)-tb(px+pi,py+pi,-0.5908,0.0962,-0.1306,-0.0507,0.0939)<0.035:
     ff1.write(str(px)+'\t')
     ff1.write(str(py)+'\n')
    if abs(ee+0.06974)<0.003 and tb(px,py,-0.5908,0.0962,-0.1306,-0.0507,0.0939)-tb(px+pi,py+pi,-0.5908,0.0962,-0.1306,-0.0507,0.0939)<0.035:
     ff2.write(str(px)+'\t')
     ff2.write(str(py)+'\n')
    if abs(ee+0.09236)<0.003 and tb(px,py,-0.5908,0.0962,-0.1306,-0.0507,0.0939)-tb(px+pi,py+pi,-0.5908,0.0962,-0.1306,-0.0507,0.0939)<0.035:
     ff3.write(str(px)+'\t')
     ff3.write(str(py)+'\n')
    if abs(ee+0.10614)<0.003 and tb(px,py,-0.5908,0.0962,-0.1306,-0.0507,0.0939)-tb(px+pi,py+pi,-0.5908,0.0962,-0.1306,-0.0507,0.0939)<0.035:
     ff4.write(str(px)+'\t')
     ff4.write(str(py)+'\n')
    if abs(ee+0.11238)<0.003 and tb(px,py,-0.5908,0.0962,-0.1306,-0.0507,0.0939)-tb(px+pi,py+pi,-0.5908,0.0962,-0.1306,-0.0507,0.0939)<0.035:
     ff5.write(str(px)+'\t')
     ff5.write(str(py)+'\n')
    if abs(ee+0.1233)<0.003 and tb(px,py,-0.5908,0.0962,-0.1306,-0.0507,0.0939)-tb(px+pi,py+pi,-0.5908,0.0962,-0.1306,-0.0507,0.0939)<0.035:
     ff6.write(str(px)+'\t')
     ff6.write(str(py)+'\n')
    if abs(ee)<0.003 and tb(px,py,-0.5908,0.0962,-0.1306,-0.0507,0.0939)-tb(px+pi,py+pi,-0.5908,0.0962,-0.1306,-0.0507,0.0939)<0.035:
     ff7.write(str(px)+'\t')
     ff7.write(str(py)+'\n')

