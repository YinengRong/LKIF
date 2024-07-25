# case 5 (time delay):
import numpy as np
from LK_Info_Flow import multi_causality_est
from LK_Info_Flow.plot_causality import causal_graph 

#generate data
xx=np.zeros([20000,5])
xx[:4,:]=np.random.normal(size=(4,5));
for i in range(4,20000-1):
    xx[i+1,0]=-0.95*np.sqrt(2)*xx[i,0]-0.9025*xx[i-1,0]+np.random.normal();
    xx[i+1,1]= 0.5*xx[i-1,0]+np.random.normal();
    xx[i+1,2]=-0.4*xx[i-2,0]+np.random.normal();
    xx[i+1,3]=-0.5*xx[i-1,0]+0.25*np.sqrt(2)*xx[i,3]+0.25*np.sqrt(2)*xx[i,4]+np.random.normal();
    xx[i+1,4]=-0.25*np.sqrt(2)*xx[i,3]+0.25*np.sqrt(2)*xx[i,4]+np.random.normal();
    
#np.savetxt('data/case5_data.txt',xx)

# Figure 2
xx=np.loadtxt('data/case5_data.txt')
xx1=np.zeros([20000-1000,7])
xx1[:,:5]=xx[1000:]
xx1[:,5]=xx[(1000-1):-1,0]#X5=X0[t-1]
xx1[:,6]=xx[(1000-2):-2,0]#X6=X0[t-2]

IF_0=multi_causality_est(X=xx1.T)
# exclude the fake causality from X(t+n) to X(t) (n>0)
a=IF_0['nIF']+0;
b=IF_0['p']+0;
a[5:,:,:]=0;
b[5:,:,:]=1;
#__________________________plot causal structure___________________________
causal_graph(causal_matrix=a,  significance=b,name='causal_structure_case5')


