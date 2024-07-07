# case 5 (time delay, latent confunders, cyclic causality):
# time delay and latent confounders

import numpy as np
from LK_Info_Flow import causal
import time 

#generate data
aii=-0.6;bii=0.3;aij=0.3;
xx=np.zeros([200000,13])

xx[:4,:]=np.random.normal(size=(4,13));
for i in range(4,200000-1):
    xx[i+1,0]=-aii*xx[i,0]    +bii*np.random.normal();
    xx[i+1,1]= aij*xx[i,0]    +aii*xx[i,1]  +aij*xx[i,2]  +aij*xx[i-2,4] \
              +aij*xx[i,5]    +aij*xx[i,6]  +aij*xx[i-1,8]+aij*xx[i,12]  \
              +bii*np.random.normal();
    xx[i+1,2]= aii*xx[i,2]    +aij*xx[i,3]  +bii*np.random.normal();
    xx[i+1,3]= aii*xx[i,3]    +bii*np.random.normal();
    xx[i+1,4]= aii*xx[i,4]    +bii*np.random.normal();
    xx[i+1,5]= aij*xx[i,1]    +aii*xx[i,5]  +bii*np.random.normal();
    xx[i+1,6]= aii*xx[i,6]    +bii*np.random.normal();
    xx[i+1,7]= aij*xx[i,6]    +aii*xx[i,7]  +bii*np.random.normal();
    xx[i+1,8]=-aii*xx[i,8]    +bii*np.random.normal();
    xx[i+1,9]= aij*xx[i,8]    +aii*xx[i,9]  +bii*np.random.normal();
    xx[i+1,10]=-aij*xx[i,1]   +aii*xx[i,10] +bii*np.random.normal();
    xx[i+1,11]= aij*xx[i,10]  +aii*xx[i,11] +bii*np.random.normal();
    xx[i+1,12]=-aij*xx[i,11]  +aii*xx[i,12] +bii*np.random.normal();

#(xx.T).shape->(13,200000)


#calculate the causality
ts=time.time()
IF_1=causal.multi_causality_est_OLS(X=xx.T,max_lag=3)
te=time.time()
print(te-ts)

#plot causal structure
## 1  causal graph
from LK_Info_Flow.plot_causality import causal_graph 
a=abs(IF_1['nIF'])
b=abs(IF_1['p'])
causal_graph(causal_matrix=IF_1['nIF'],  significance=IF_1['p'],  c_threshold=0.001,      s_threshold=0.01,   name='model')
#             the IF or nIF matrix            err_e90-99,p        threshold of causality  and significance    name of graph





## 2  heatmap
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
c_threshold=0.001
s_threshold=0.01
a=abs(IF_1['nIF'])
b=abs(IF_1['p'])
a[np.where(b>s_threshold)]=np.NaN
nM=a.shape[0]

a[np.where(a<c_threshold)]=np.NaN

c=np.reshape(a[:,:,::1],[13,13*3]).T
for i in range(nM):
    c[i,i]=np.NaN

data=pd.DataFrame(c.transpose())

dict_={'orientation':'vertical',"label":"normalized information flow (%)",\
       "ticklocation":"right", "alpha":0.8,"cmap":"cmap"}

plt.figure(figsize=(12, 3))
xticklabels=['0'    ,'1'    ,'2'    ,'3'    ,'4'    ,'5'    ,'6'    ,'7'    ,'8'    ,'9'    ,'10'    ,'11'    ,'12'    , \
             '0(-1)','1(-1)','2(-1)','3(-1)','4(-1)','5(-1)','6(-1)','7(-1)','8(-1)','9(-1)','10(-1)','11(-1)','12(-1)', \
             '0(-2)','1(-2)','2(-2)','3(-2)','4(-2)','5(-2)','6(-2)','7(-2)','8(-2)','9(-2)','10(-2)','11(-2)','12(-2)', ]
plot=sns.heatmap(data*100,linewidths=0.8,linecolor='black',cmap="jet",vmax=15,vmin=0.2,cbar_kws=dict_,xticklabels=xticklabels)

plt.xlabel('cause');
plt.ylabel('effect');