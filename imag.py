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
ff1=open('imag.dat','w')

dd=1

for il in range(NN):
    tt=10.0*il/NN
    if tt<=3:
     gamma=tt/(3+(3-tt))
     gamma1=tt/(3+2*(3-tt))
     gamma2=tt/(3+3*(3-tt))
    else: 
      gamma=tt*1.0/3
      gamma1=tt*1.0/3
      gamma2=tt*1.0/3
    ff1.write(str(tt)+'\t')
    ff1.write(str(gamma)+'\t')
    ff1.write(str(gamma1)+'\t')
    ff1.write(str(gamma2)+'\n')
