from numpy import *
from pylab import *
from math import *
import matplotlib
import numpy as np
from matplotlib import cm
from matplotlib import pyplot as plt
from numpy import linalg as LA
from math import exp

NM=5000
NN=30
NB=200
DL=5.0
ff1=open('sc.dat','w')


def tb(px,py,t,t1,t2,t3,t4):
    EE=0.5*t*(cos(px)+cos(py))+t1*cos(px)*cos(py)+0.5*t2*(cos(2*px)+cos(2*py))+0.5*t3*(cos(2*px)*cos(py)+cos(px)*cos(2*py))+t4*cos(2*px)*cos(2*py)
    return EE


def FF(kx,ky,T):
    e=tb(kx+pi,ky+pi,-0.5908,0.0962,-0.1306,-0.0507,0.0939)-tb(kx,ky,-0.5908,0.0962,-0.1306,-0.0507,0.0939)+0.04
    if e!=0
        ta=tanh(e/T)/e
    else: ta=1.0/T
    return ta.real


def FF1(kx,ky,T):
    e=tb(kx+pi,ky+pi,-0.5908,0.0962,-0.1306,-0.0507,0.0939)-tb(kx,ky,-0.5908,0.0962,-0.1306,-0.0507,0.0939)
    if e!=0
        ta=tanh(e/T)/e
    else: ta=1.0/T
    return ta.real

def SS(kx,ky,qx,qy,T):
    e=tb(kx+qx,ky+qy,-0.5908,0.0962,-0.1306,-0.0507,0.0939)-tb(kx,ky,-0.5908,0.0962,-0.1306,-0.0507,0.0939)+0.04*(0.5-0.25*cos(qx)-0.25*cos(qy))
    if e!=0
        ta=tanh(e/T)/e
    else: ta=1.0/T
    return ta.real

def SSP(kx,ky,qx,qy,T):
    e=tb(kx+qx+pi,ky+qy+pi,-0.5908,0.0962,-0.1306,-0.0507,0.0939)-tb(kx,ky,-0.5908,0.0962,-0.1306,-0.0507,0.0939)+0.04*(1.5+0.25*cos(qx)+0.25*cos(qy))
    if e!=0
        ta=tanh(e/T)/e
    else: ta=1.0/T
    return ta.real


def SS1(kx,ky,qx,qy,T):
    e=tb(kx+qx,ky+qy,-0.5908,0.0962,-0.1306,-0.0507,0.0939)-tb(kx,ky,-0.5908,0.0962,-0.1306,-0.0507,0.0939)
    if e!=0
        ta=tanh(e/T)/e
    else: ta=1.0/T
    return ta.real

def SSP1(kx,ky,qx,qy,T):
    e=tb(kx+qx+pi,ky+qy+pi,-0.5908,0.0962,-0.1306,-0.0507,0.0939)-tb(kx,ky,-0.5908,0.0962,-0.1306,-0.0507,0.0939)
    if e!=0
        ta=tanh(e/T)/e
    else: ta=1.0/T
    return ta.real

def DEL(kx,ky):
   DD=cos(kx)-cos(ky)
   return DD


def CC(kx,ky,T):
    xi=2*tb(kx,ky,-0.5908,0.0962,-0.1306,-0.0507,0.0939)**2*(tb(kx+pi,ky+pi,-0.5908,0.0962,-0.1306,-0.0507,0.0939)-tb(kx,ky,-0.5908,0.0962,-0.1306,-0.0507,0.0939)+0.04)/((tb(kx+pi,ky+pi,-0.5908,0.0962,-0.1306,-0.0507,0.0939)-tb(kx,ky,-0.5908,0.0962,-0.1306,-0.0507,0.0939)+0.04)**2-(tb(kx+pi,ky+pi,-0.5908,0.0962,-0.1306,-0.0507,0.0939)-tb(kx,ky,-0.5908,0.0962,-0.1306,-0.0507,0.0939))**2)*(FF(kx,ky,T)-FF1(kx,ky,T))*DEL(kx,ky)**2
    return xi

def CC1(kx,ky,qx,qy,T):
    xi=2*tb(kx,ky,-0.5908,0.0962,-0.1306,-0.0507,0.0939)**2*(tb(kx+qx,ky+qy,-0.5908,0.0962,-0.1306,-0.0507,0.0939)-tb(kx,ky,-0.5908,0.0962,-0.1306,-0.0507,0.0939)+0.04*(0.5-0.25*cos(px)-0.25*cos(py)))/((tb(kx+qx,ky+qy,-0.5908,0.0962,-0.1306,-0.0507,0.0939)-tb(kx,ky,-0.5908,0.0962,-0.1306,-0.0507,0.0939)+0.04*(0.5-0.25*cos(px)-0.25*cos(py)))**2-(tb(kx+qx,ky+qy,-0.5908,0.0962,-0.1306,-0.0507,0.0939)-tb(kx,ky,-0.5908,0.0962,-0.1306,-0.0507,0.0939))**2)*(SS(kx,ky,qx,qy,T)-SS1(kx,ky,qx,qy,T))*DEL(kx,ky)*DEL(kx+qx,ky+qy)
    return xi

def CC2(kx,ky,qx,qy,T):
    xi=-2*tb(kx,ky,-0.5908,0.0962,-0.1306,-0.0507,0.0939)**2*(tb(kx+qx+pi,ky+qy+pi,-0.5908,0.0962,-0.1306,-0.0507,0.0939)-tb(kx,ky,-0.5908,0.0962,-0.1306,-0.0507,0.0939)+0.04*(1.5+0.25*cos(px)+0.25*cos(py)))/((tb(kx+qx+pi,ky+qy+pi,-0.5908,0.0962,-0.1306,-0.0507,0.0939)-tb(kx,ky,-0.5908,0.0962,-0.1306,-0.0507,0.0939)+0.04*(1.5+0.25*cos(px)+0.25*cos(py)))**2-(tb(kx+qx+pi,ky+qy+pi,-0.5908,0.0962,-0.1306,-0.0507,0.0939)-tb(kx,ky,-0.5908,0.0962,-0.1306,-0.0507,0.0939))**2)*(SSP(kx,ky,qx,qy,T)-SSP1(kx,ky,qx,qy,T))*DEL(kx,ky)*DEL(kx+qx,ky+qy)
    return xi

