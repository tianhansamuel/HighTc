from numpy import *
from pylab import *
from math import *
import matplotlib
import numpy as np
from matplotlib import cm
from matplotlib import pyplot as plt
from numpy import linalg as LA
from math import exp


NN=400
NB=15
N1=100

J1=-1*0.67
Delta=4*0.67
GG=40
dt=20
seuil=0.1

ff1=open('dd.dat','w')

MM=zeros((2,2),complex)

#############################################  bare propagator 


def tb(px,py):
    tb=-2*(cos(px)+cos(py))+cos(px)*cos(py)
    return tb

def om(px,py,alpha):
    if tb(px+pi,py+pi)<=tb(px,py):
     omega=tb(px,py)+tb(px+pi,py+pi)+sqrt((tb(px,py)+tb(px+pi,py+pi))**2-4*tb(px,py)*tb(px+pi,py+pi)+4*tb(px+pi,py+pi)**2-4*tb(px,py)*alpha)
     omega1=tb(px,py)+tb(px+pi,py+pi)-sqrt((tb(px,py)+tb(px+pi,py+pi))**2-4*tb(px,py)*tb(px+pi,py+pi)+4*tb(px+pi,py+pi)**2-4*tb(px,py)*alpha)
    else: omega=omega1=2*tb(px,py)
    return (0.5*omega,0.5*omega1)

def KY(EF,kx):
    AA=-(EF+2*cos(kx))/(2-cos(kx))
    if abs(AA)<=1:
       ky=arccos(AA)
    else: ky=-1000
    return ky

EE=0.5

for io in range(NN):
    kx=2*pi*io/NN-pi
    ky=KY(EE,kx)
    if ky!=-1000:
     ff1.write(str(kx)+'\t')
     ff1.write(str(ky)+'\t')
     ff1.write(str(cos(kx)-cos(ky))+'\t')
     ff1.write(str(om(kx,ky,0.1)[0]-EE)+'\t')
     ff1.write(str(om(kx,ky,0.1)[1]-EE)+'\n')








