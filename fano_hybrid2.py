from numpy import *
from math import *


NN=60
N1=60
NN=300


ff1=open('hybrid_mm_2.dat','w')

bogo=loadtxt('test_0.01_refmm_2.dat', dtype =float)   
GG=loadtxt('test_single_0.05_refmm.dat', dtype =float)

SS1=shape(GG)

SS2=shape(bogo)


S1=SS1[0]
d1=SS1[1]

print SS1
print SS2

s1=3
s2=100

dd=int(0.15/0.6*S1)
for i in range(dd+1,S1):
    omega1=bogo[i][0]
    GB1=bogo[i][1]
    N1=bogo[i][11]

    GB2=bogo[i][2]
    N2=bogo[i][12]


    GB3=bogo[i][3]
    N3=bogo[i][13]

    GB4=bogo[i][4]
    N4=bogo[i][14]

    GB5=bogo[i][5]
    N5=bogo[i][15]

    GB6=bogo[i][6]
    N6=bogo[i][16]

    GB7=bogo[i][7]
    N7=bogo[i][17]

    GB8=bogo[i][8]
    N8=bogo[i][18]

    GB9=bogo[i][9]
    N9=bogo[i][19]

    GB10=bogo[i][10]
    N10=bogo[i][20]


    GS1=GG[i-dd][1]
    NS1=GG[i-dd][11]

    GS2=GG[i-dd][2]
    NS2=GG[i-dd][12]

    GS3=GG[i-dd][3]
    NS3=GG[i-dd][13]

    GS4=GG[i-dd][4]
    NS4=GG[i-dd][14]

    GS5=GG[i-dd][5]
    NS5=GG[i-dd][15]

    GS6=GG[i-dd][6]
    NS6=GG[i-dd][16]

    GS7=GG[i-dd][7]
    NS7=GG[i-dd][17]

    GS8=GG[i-dd][8]
    NS8=GG[i-dd][18]
    
    GS9=GG[i-dd][9]
    NS9=GG[i-dd][19]
        
    GS10=GG[i-dd][10]
    NS10=GG[i-dd][20]
    
    GF1=(s1*GB1*N1-s2*GS1*NS1)/(N1+NS1)/400

    GF2=(s1*GB2*N2-s2*GS2*NS2)/(N2+NS2)/400

    GF3=(s1*GB3*N3-s2*GS3*NS3)/(N3+NS3)/400

    GF4=(s1*GB4*N4-s2*GS4*NS4)/(N4+NS4)/400

    GF5=(s1*GB5*N5-s2*GS5*NS5)/(N5+NS5)/400

    GF6=(s1*GB6*N6-s2*GS6*NS6)/(N6+NS6)/400

    GF7=(s1*GB7*N7-s2*GS7*NS7)/(N7+NS7)/400

    GF8=(s1*GB8*N8-s2*GS8*NS8)/(N8+NS8)/400

    GF9=(s1*GB9*N9-s2*GS9*NS9)/(N9+NS9)/400

    GF10=(s1*GB10*N10-s2*GS10*NS10)/(N10+NS10)/400

    ff1.write(str(omega1)+'\t')
    ff1.write(str(GF1)+'\t')
    ff1.write(str(GF2)+'\t')
    ff1.write(str(GF3)+'\t')
    ff1.write(str(GF4)+'\t')
    ff1.write(str(GF5)+'\t')
    ff1.write(str(GF6)+'\t')
    ff1.write(str(GF7)+'\t')
    ff1.write(str(GF8)+'\t')
    ff1.write(str(GF9)+'\t')
    ff1.write(str(GF10)+'\n')
