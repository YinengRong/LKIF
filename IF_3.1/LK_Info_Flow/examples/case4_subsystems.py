# case 4 (subsystems):
import numpy as np
from LK_Info_Flow import causal


#generate dataset
a11=-0.5;a21= 0.5;a31= 0.2; b1=1.0;
a12= 0.0;a22=-0.2;a32=-0.6; b2=1.0;
a13=-0.2;a23= 0.4;a33=-0.2; b3=1.0;

b11=-0.2;b21=-0.5;b31= 0.0; b4=1.0;
b12= 0.5;b22=-0.6;b32= 0.4; b5=1.0;
b13=-0.1;b23=-0.4;b33=-0.5; b6=1.0;

esp1=0.5
esp3=0.0


xx=np.zeros([20001,6])
xx[0,:]=np.array([0.5,0.5,0.5,0.5,0.5,0.5])
for i in range(20000):
    xx[i+1,0]=a11*xx[i,0]+a21*xx[i,1]+a31*xx[i,2]+b1*np.random.normal();
    xx[i+1,1]=a12*xx[i,0]+a22*xx[i,1]+a32*xx[i,2]+b2*np.random.normal();
    xx[i+1,2]=a13*xx[i,0]+a23*xx[i,1]+a33*xx[i,2]+b3*np.random.normal()+esp3*xx[i,5];
    xx[i+1,3]=b11*xx[i,3]+b21*xx[i,4]+b31*xx[i,5]+b1*np.random.normal()-esp1*xx[i,0];
    xx[i+1,4]=b12*xx[i,3]+b22*xx[i,4]+b32*xx[i,5]+b2*np.random.normal();
    xx[i+1,5]=b13*xx[i,3]+b23*xx[i,4]+b33*xx[i,5]+b3*np.random.normal();
#_________________________________________________________________________________


#calculate the causality
ind=[3,6]; #X0,X1,X2 subsystem A;X3,X4,X5 subsystem B
IF_g=causal.groups_est(xx=xx,ind=ind)
print('TA->B: %8.4f  TB->A: %8.4f'%(IF_g['TAB'],IF_g['TBA']))