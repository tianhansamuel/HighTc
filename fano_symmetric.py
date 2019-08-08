from numpy import *
from math import *


NN=60
N1=60
NN=300


ff1=open('gg.dat','w')

hybrid=loadtxt('hybrid_mm_1.dat', dtype =float)   


SS1=shape(hybrid)

S1=SS1[0]
d1=SS1[1]

print SS1

for i in range(1,S1):
   o1=hybrid[i][0]
   g1=hybrid[i][1]/hybrid[S1-i][1]
   g2=hybrid[i][2]/hybrid[S1-i][2]
   g3=hybrid[i][3]/hybrid[S1-i][3]
   g4=hybrid[i][4]/hybrid[S1-i][4]
   g5=hybrid[i][5]/hybrid[S1-i][5]
   g6=hybrid[i][6]/hybrid[S1-i][6]
   g7=hybrid[i][7]/hybrid[S1-i][7]
   g8=hybrid[i][8]/hybrid[S1-i][8]
   g9=hybrid[i][9]/hybrid[S1-i][9]
   g10=hybrid[i][10]/hybrid[S1-i][10]   
   ff1.write(str(o1)+'\t')
   ff1.write(str(g1)+'\t')
   ff1.write(str(g2)+'\t')
   ff1.write(str(g3)+'\t')
   ff1.write(str(g4)+'\t')
   ff1.write(str(g5)+'\t')
   ff1.write(str(g6)+'\t')
   ff1.write(str(g7)+'\t')
   ff1.write(str(g8)+'\t')
   ff1.write(str(g9)+'\t')
   ff1.write(str(g10)+'\n')
