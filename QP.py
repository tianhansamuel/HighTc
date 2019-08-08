from numpy import *
from pylab import *
from math import *
import matplotlib
import numpy as np
from matplotlib import cm
from matplotlib import pyplot as plt
from numpy import linalg as LA
from math import exp


NN=5000
NB=100
DL=5.0
ff1=open('coherence.dat','w')


#############################################  bare propagator 

def GG(omega,DD,TT,la,cc):
    GS=0
    for it in range(NB):
       om=(2*it+1)*TT
       SIG=0
       for ii1 in range(NB):
          EE=DD-cc*ii1*1.0/NB
          SIG+=1.0/(omega+1j*om+EE)*1.0/(NB*cc)
       GS+=1.0/(omega-1+1j*om-la**2*SIG)*1.0/NB
    return GS


for ii in range(NN):
    omega=6.0*ii/NN
    print ii
    G1=-GG(omega,1,0.01,2,0.5).imag
    G2=-GG(omega,1,0.05,2,0.5).imag
    G3=-GG(omega,1,0.1,2,0.5).imag
    G4=-GG(omega,1,0.2,2,0.5).imag
    G5=-GG(omega,1,0.8,2,0.5).imag
    G6=-GG(omega,1,1,2,0.5).imag
    ff1.write(str(omega)+'\t')
    ff1.write(str(G1)+'\t')
    ff1.write(str(G2)+'\t')
    ff1.write(str(G3)+'\t')
    ff1.write(str(G4)+'\t')
    ff1.write(str(G5)+'\t')
    ff1.write(str(G6)+'\n')

