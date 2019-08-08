from numpy import *
from math import *


NN=60
N1=60
NN=300


ff1=open('hybrid_mm1.dat','w')

bogo=loadtxt('test.dat', dtype =float)   
GG=loadtxt('test_single_2.dat', dtype =float)

SS1=shape(bogo)

S1=SS1[0]
d1=SS1[1]

s1=10
s2=1
for i in range(NN):
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


    GS1=GG[i][1]
    NS1=GG[i][11]

    GS2=GG[i][2]
    NS2=GG[i][12]

    GS3=GG[i][3]
    NS3=GG[i][13]

    GS4=GG[i][4]
    NS4=GG[i][14]

    GS5=GG[i][5]
    NS5=GG[i][15]

    GS6=GG[i][6]
    NS6=GG[i][16]

    GS7=GG[i][7]
    NS7=GG[i][17]

    GS8=GG[i][8]
    NS8=GG[i][18]
    
    GS9=GG[i][9]
    NS9=GG[i][19]
        
    GS10=GG[i][10]
    NS10=GG[i][20]
    
    GF1=(s1*GB1*N1-s2*GS1*NS1*600)/(N1+NS1)/400

    GF2=(s1*GB2*N2-s2*GS2*NS2*600)/(N2+NS2)/400

    GF3=(s1*GB3*N3-s2*GS3*NS3*600)/(N3+NS3)/400

    GF4=(s1*GB4*N4-s2*GS4*NS4*600)/(N4+NS4)/400

    GF5=(s1*GB5*N5-s2*GS5*NS5*600)/(N5+NS5)/400

    GF6=(s1*GB6*N6-s2*GS6*NS6*600)/(N6+NS6)/400

    GF7=(s1*GB7*N7-s2*GS7*NS7*600)/(N7+NS7)/400

    GF8=(s1*GB8*N8-s2*GS8*NS8*600)/(N8+NS8)/400

    GF9=(s1*GB9*N9-s2*GS9*NS9*600)/(N9+NS9)/400

    GF10=(s1*GB10*N10-s2*GS10*NS10*600)/(N10+NS10)/400

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
