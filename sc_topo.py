from numpy import *
from math import *


NN=200
NB=50

J1=-1*0.67
Delta=4*0.67
GG=40
dt=20
seuil=0.1

ff1=open('sc_temp1.dat','w')


MM=zeros((2,2),complex)

#############################################  bare propagator 


def tb(px,py,t,t1,t2,t3,t4):
    EE=0.5*t*(cos(px)+cos(py))+t1*cos(px)*cos(py)+0.5*t2*(cos(2*px)+cos(2*py))+0.5*t3*(cos(2*px)*cos(py)+cos(px)*cos(2*py))+t4*cos(2*px)*cos(2*py)
    return EE

def DE(px,py,t,t1,t2,t3,t4,ww):
    omega=tb(px+pi,py+pi,t,t1,t2,t3,t4)-tb(px,py,t,t1,t2,t3,t4)+ww
    return omega


def om(alpha,dd):
    EG=-(1-12*alpha**2)/(1+12*alpha**2)*(1-dd)**2+dd*(1-dd)**4*3*(1-alpha)/(1+12*alpha**2)+8*dd*(1-dd)**5*(1-alpha)/(1+12*alpha**2)
    return EG

def AMAX(dd):
    EG=100
    aa=0
    for ia in range(NB):
       alpha=0.3*ia/NB
       if om(alpha,dd)<EG:
          aa=alpha
          EG=om(alpha,dd)
    RES=4*aa/(1-12*aa**2)*0.129
    return (RES,aa)

def TT(e,Temp):
    if e!=0:
       cc=tanh(e*1.0/(2*Temp))/e
    else: cc=1.0/Temp
    return cc

def GAMMA(kx,ky,dd,T,t,t1,t2,t3,t4,aa,ww):
    DF=DE(kx,ky,t,t1,t2,t3,t4,ww)
    xi=tb(kx,ky,t,t1,t2,t3,t4)
    xip=tb(kx+pi,ky+pi,t,t1,t2,t3,t4)
    SS=0.5*(cos(kx)-cos(ky))**2*(1-12*aa**2)/(1+12*aa**2)*(1-dd)**2*(xip)**2*ww/(ww**2-(xip-xi)**2)*(TT(ww,T)-TT(xip-xi,T))
    return SS


def dd(EF,t,t1,t2,t3,t4):
      NF=0
      for ix in range(NB+1):
        for iy in range(NB+1):
         kx=-pi+2*pi*ix/NB
         ky=-pi+2*pi*iy/NB
         if tb(kx,ky,t,t1,t2,t3,t4)<=EF:
           NF+=1
      dl=NF*1.0/(NB+1)**2
      return dl       

def CCT(EF,T,t,t1,t2,t3,t4,dop,aa,ww):
    NF=0
    CT=0
    for ix in range(NN):
      for iy in range(NN):
        kx=-pi+2*pi*ix/NN
        ky=-pi+2*pi*iy/NN    
        if abs(tb(kx,ky,t,t1,t2,t3,t4)-EF)<0.01 and tb(kx,ky,t,t1,t2,t3,t4)-tb(kx+pi,ky+pi,t,t1,t2,t3,t4)<ww:
           NF+=1
           CT+=GAMMA(kx, ky, dop, T, t,t1,t2,t3,t4,aa,ww)
    CC=CT*1.0/NF
    CC1=CT*1.0/(NN**2)
    return (CC,CC1)


for it in range(1,NN):
    print it
    T=0.04*it/NN
    EF=0.02
    dop=2*(0.5-dd(EF,-0.5908,0.0962,-0.1306,-0.0507,0.0939)+0.1)
    aa=AMAX(dop)[1]
    ww=AMAX(dop)[0]
    print ww
    gg0=CCT(EF,T,-0.5908,0.0962,-0.1306,-0.0507,0.0939,dop,aa,ww)[0]
    ggt0=CCT(EF,T,-0.5908,0.0962,-0.1306,-0.0507,0.0939,dop,aa,ww)[1]


    EF=0
    dop=2*(0.5-dd(EF,-0.5908,0.0962,-0.1306,-0.0507,0.0939)+0.1)
    aa=AMAX(dop)[1]
    ww=AMAX(dop)[0]
    print ww
    gg=CCT(EF,T,-0.5908,0.0962,-0.1306,-0.0507,0.0939,dop,aa,ww)[0]
    ggt=CCT(EF,T,-0.5908,0.0962,-0.1306,-0.0507,0.0939,dop,aa,ww)[1]

    EF=-0.02
    dop=2*(0.5-dd(EF,-0.5908,0.0962,-0.1306,-0.0507,0.0939)+0.1)
    aa=AMAX(dop)[1]
    ww=AMAX(dop)[0]
    print ww
    gg1=CCT(EF,T,-0.5908,0.0962,-0.1306,-0.0507,0.0939,dop,aa,ww)[0]
    gg1t=CCT(EF,T,-0.5908,0.0962,-0.1306,-0.0507,0.0939,dop,aa,ww)[1]

    EF=-0.04
    dop=2*(0.5-dd(EF,-0.5908,0.0962,-0.1306,-0.0507,0.0939)+0.1)
    aa=AMAX(dop)[1]
    ww=AMAX(dop)[0]
    print ww
    gg2=CCT(EF,T,-0.5908,0.0962,-0.1306,-0.0507,0.0939,dop,aa,ww)[0]
    gg2t=CCT(EF,T,-0.5908,0.0962,-0.1306,-0.0507,0.0939,dop,aa,ww)[1]

    EF=-0.06
    dop=2*(0.5-dd(EF,-0.5908,0.0962,-0.1306,-0.0507,0.0939)+0.1)
    aa=AMAX(dop)[1]
    ww=AMAX(dop)[0]
    print ww
    gg3=CCT(EF,T,-0.5908,0.0962,-0.1306,-0.0507,0.0939,dop,aa,ww)[0]
    gg3t=CCT(EF,T,-0.5908,0.0962,-0.1306,-0.0507,0.0939,dop,aa,ww)[1]


    EF=-0.08
    dop=2*(0.5-dd(EF,-0.5908,0.0962,-0.1306,-0.0507,0.0939)+0.1)
    aa=AMAX(dop)[1]
    ww=AMAX(dop)[0]
    print ww
    gg4=CCT(EF,T,-0.5908,0.0962,-0.1306,-0.0507,0.0939,dop,aa,ww)[0]
    gg4t=CCT(EF,T,-0.5908,0.0962,-0.1306,-0.0507,0.0939,dop,aa,ww)[1]

    EF=-0.10
    dop=2*(0.5-dd(EF,-0.5908,0.0962,-0.1306,-0.0507,0.0939)+0.1)
    aa=AMAX(dop)[1]
    ww=AMAX(dop)[0]
    print ww
    gg5=CCT(EF,T,-0.5908,0.0962,-0.1306,-0.0507,0.0939,dop,aa,ww)[0]
    gg5t=CCT(EF,T,-0.5908,0.0962,-0.1306,-0.0507,0.0939,dop,aa,ww)[1]

    EF=-0.12
    dop=2*(0.5-dd(EF,-0.5908,0.0962,-0.1306,-0.0507,0.0939)+0.1)
    aa=AMAX(dop)[1]
    ww=AMAX(dop)[0]
    print ww
    gg6=CCT(EF,T,-0.5908,0.0962,-0.1306,-0.0507,0.0939,dop,aa,ww)[0]
    gg6t=CCT(EF,T,-0.5908,0.0962,-0.1306,-0.0507,0.0939,dop,aa,ww)[1]


    EF=-0.14
    dop=2*(0.5-dd(EF,-0.5908,0.0962,-0.1306,-0.0507,0.0939)+0.1)
    aa=AMAX(dop)[1]
    ww=AMAX(dop)[0]
    print ww
    gg7=CCT(EF,T,-0.5908,0.0962,-0.1306,-0.0507,0.0939,dop,aa,ww)[0]
    gg7t=CCT(EF,T,-0.5908,0.0962,-0.1306,-0.0507,0.0939,dop,aa,ww)[1]
    ff1.write(str(T)+'\t')
    ff1.write(str(gg0)+'\t')
    ff1.write(str(gg)+'\t')
    ff1.write(str(gg1)+'\t')
    ff1.write(str(gg2)+'\t')
    ff1.write(str(gg3)+'\t')
    ff1.write(str(gg4)+'\t')
    ff1.write(str(gg5)+'\t')
    ff1.write(str(gg6)+'\t')
    ff1.write(str(gg7)+'\n')
