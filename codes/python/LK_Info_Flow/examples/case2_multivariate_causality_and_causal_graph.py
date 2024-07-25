# case 2 (multivariate causality):
import numpy as np
np.random.seed(0)  


# generate dataset
a11=0;   a21=0;   a31=-0.6;  a41=0;   a51=-0.0;  a61=0;   b1=0.1;
a12=-0.5;a22=0;   a32=-0.0;  a42=0;   a52=0.0;   a62=0.8; b2=0.7;
a13=0;   a23=0.7; a33=-0.6;  a43=0;   a53=-0.0;  a63=0;   b3=0.5;
a14=0;   a24=0;   a34=-0.;   a44=0.7; a54=0.4;   a64=0;   b4=0.2;
a15=0;   a25=0;   a35=0;     a45=0.2; a55=0.0;   a65=0.7; b5=0.8;
a16=0;   a26=0;   a36=0;     a46=0;   a56=0.0;   a66=-0.5;b6=0.3;

M=6;
xx=np.zeros([100001,M])
xx[0]=np.array([0.4,0.5,0.6,0.7,0.6,0.7])
for i in range(100000):
    xx[i+1,0]=a11*xx[i,0]+a21*xx[i,1]+a31*xx[i,2]+a41*xx[i,3]+a51*xx[i,4]+a61*xx[i,5]+b1+np.random.normal();
    xx[i+1,1]=a12*xx[i,0]+a22*xx[i,1]+a32*xx[i,2]+a42*xx[i,3]+a52*xx[i,4]+a62*xx[i,5]+b2+np.random.normal();
    xx[i+1,2]=a13*xx[i,0]+a23*xx[i,1]+a33*xx[i,2]+a43*xx[i,3]+a53*xx[i,4]+a63*xx[i,5]+b3+np.random.normal();
    xx[i+1,3]=a14*xx[i,0]+a24*xx[i,1]+a34*xx[i,2]+a44*xx[i,3]+a54*xx[i,4]+a64*xx[i,5]+b4+np.random.normal();
    xx[i+1,4]=a15*xx[i,0]+a25*xx[i,1]+a35*xx[i,2]+a45*xx[i,3]+a55*xx[i,4]+a65*xx[i,5]+b5+np.random.normal();
    xx[i+1,5]=a16*xx[i,0]+a26*xx[i,1]+a36*xx[i,2]+a46*xx[i,3]+a56*xx[i,4]+a66*xx[i,5]+b6+np.random.normal();
#np.savetxt('data/case2_data.txt',xx)

##------calculate the causality---------------------------------
xx=np.loadtxt('data/case2_data.txt')
xx=xx[10000:].T;
from LK_Info_Flow import multi_causality_est
cau2=multi_causality_est(X=xx)#X [ndim,ntime_steps]


#information flow from column to raw
IF = np.squeeze(cau2['IF']) 
#normalized information flow
nIF= np.squeeze(cau2.get('nIF'));
#significant test: confidence levels of 90/95/99% err_e90/err_e95/err_e99;
err= np.squeeze(cau2.get('err_e99'));
#significant test (p-value)
p=np.squeeze(cau2.get('p'));

##---------causal graph---------------------------------------
from LK_Info_Flow.plot_causality import causal_graph 
causal_graph(causal_matrix=nIF,  name='causal_structure_case2')

