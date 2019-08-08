from numpy import *
from pylab import *
from math import *
import matplotlib
import numpy as np
from matplotlib import cm
from matplotlib import pyplot as plt
from numpy import linalg as LA
from math import exp


NN=1000
NB=10000
DL=5.0
ff1=open('FS10.dat','w')
ff2=open('FS12.dat','w')
ff3=open('FS14.dat','w')
ff4=open('FS15.dat','w')
ff5=open('FS17.dat','w')
ff6=open('FS18.dat','w')
ff7=open('FS20.dat','w')

dd=1


def tb(px,py,t,t1,t2,t3,t4):
    EE=0.5*t*(cos(px)+cos(py))+t1*cos(px)*cos(py)+0.5*t2*(cos(2*px)+cos(2*py))+0.5*t3*(cos(2*px)*cos(py)+cos(px)*cos(2*py))+t4*cos(2*px)*cos(2*py)
    return EE

def dd(EF):
    SS=0
    for ix in range(NN):
      for iy in range(NN):
         px=2*pi*ix/NN-pi
         py=2*pi*iy/NN-pi
         if tb(px,py,-0.5908,0.0962,-0.1306,-0.0507,0.0939)<EF:
            SS+=1.0/(NN**2)
    return SS


for ix in range(NN):
  for iy in range(NN):
    px=2*pi*ix/NN-pi
    py=2*pi*iy/NN-pi
    if abs(tb(px,py,-0.5908,0.0962,-0.1306,-0.0507,0.0939)-0.064056)<0.01:   #### 0.1
       ff1.write(str(px)+'\t')
       ff1.write(str(py)+'\n')
    if abs(tb(px,py,-0.5908,0.0962,-0.1306,-0.0507,0.0939)-0.1209504)<0.01:  #### 0.15
       ff2.write(str(px)+'\t')
       ff2.write(str(py)+'\n')
    if abs(tb(px,py,-0.5908,0.0962,-0.1306,-0.0507,0.0939)-0.1562016)<0.01:  #### 0.17
       ff3.write(str(px)+'\t')
       ff3.write(str(py)+'\n')
    if abs(tb(px,py,-0.5908,0.0962,-0.1306,-0.0507,0.0939)-0.1682544)<0.01:  ##### 0.18
       ff4.write(str(px)+'\t')
       ff4.write(str(py)+'\n')
    if abs(tb(px,py,-0.5908,0.0962,-0.1306,-0.0507,0.0939)-0.1801776)<0.01:  ##### 0.2
       ff5.write(str(px)+'\t')
       ff5.write(str(py)+'\n')
    if abs(tb(px,py,-0.5908,0.0962,-0.1306,-0.0507,0.0939)-0.1801776)<0.01:  ##### 0.2
       ff6.write(str(px)+'\t')
       ff6.write(str(py)+'\n')
    if abs(tb(px,py,-0.5908,0.0962,-0.1306,-0.0507,0.0939)-0.1801776)<0.01:  ##### 0.2
       ff7.write(str(px)+'\t')
       ff7.write(str(py)+'\n')



