# case 3(panel data, discontinuous time series or ensemble data):
## 3.1________________________________________________
import numpy as np
from LK_Info_Flow import multi_causality_est
import time
from tqdm import tqdm 


# generate data
a11=0.3;a21=0;  a31=0;  a41=0;  b1=0.4;
a12=0.5;a22=0.7;a32=0.1;a42=0;  b2=0.5;
a13=0;  a23=0.3;a33=0.5;a43=0;  b3=0.6;
a14=0.2;a24=0.4;a34=0.3;a44=0.1;b4=0.3;
case_num=1000;
xx=np.zeros([case_num,1001,4])
xx[:,0,:]=np.array([0.4,0.5,0.6,0.7])
print('----------generating the data----------')
for icase in tqdm(range(case_num)):
    if icase!=0:
        xx[icase,0,:]=xx[icase,0,:]+np.random.normal(0,1,4)
    for i in range(1000):
        xx[icase,i+1,0]=a11*xx[icase,i,0]+a21*xx[icase,i,1]+a31*xx[icase,i,2]+a41*xx[icase,i,3]+b1*np.random.normal();
        xx[icase,i+1,1]=a12*xx[icase,i,0]+a22*xx[icase,i,1]+a32*xx[icase,i,2]+a42*xx[icase,i,3]+b2*np.random.normal();
        xx[icase,i+1,2]=a13*xx[icase,i,0]+a23*xx[icase,i,1]+a33*xx[icase,i,2]+a43*xx[icase,i,3]+b3*np.random.normal();
        xx[icase,i+1,3]=a14*xx[icase,i,0]+a24*xx[icase,i,1]+a34*xx[icase,i,2]+a44*xx[icase,i,3]+b4*np.random.normal();


# build the panel data
## Select a segment of 10 time length for each case and combine it into panel data
## and build the corresponding temperal index
X = np.zeros([10*case_num,4]);
t = np.zeros([10*case_num,]);
for j in range(case_num):
    i=int(np.floor(np.random.uniform()*case_num));
    X[10*j:10*(j+1),:]=xx[i,-10:,:];
    t[10*j:10*(j+1),]=np.linspace(0,9,10);


#X=np.loadtxt('E:\\BaiduSyncdisk\\papers\\author\\package\\R\\case3_data_X.txt')
#t=np.loadtxt('E:\\BaiduSyncdisk\\papers\\author\\package\\R\\case3_data_t.txt')
#calculate the causality (panel data)
print('start calculate causality:')
time_start=time.time()
IF_panel=multi_causality_est(X=X[:,:4].T,series_temporal_order=t)

T21=IF_panel.get('IF').squeeze()
err=IF_panel.get('err_e90').squeeze()
time_end=time.time()
print('est_panel: T1->2: %8.4f e90: %8.4f'%(T21[1,0], err[1,0]))
print('est_panel: T2->1: %8.4f e90: %8.4f'%(T21[0,1], err[0,1]))
print('time cost: %8.4f s'%(time_end-time_start))


#calculate the causality in one case 
time_start=time.time()
IF_time_series=multi_causality_est(X=xx[1,:,:].T)
T21=IF_time_series.get('IF').squeeze()
err=IF_time_series.get('err_e99').squeeze()
time_end=time.time()
print('est_time_series: T1->2: %8.4f e90: %8.4f'%(T21[1,0], err[1,0]))
print('est_time_series: T2->1: %8.4f e90: %8.4f'%(T21[0,1], err[0,1]))
print('time cost: %8.4f s'%(time_end-time_start))

## 3.2 ___________________________________________________________________________________________________________
#For time series data with missing measurements in time or a set of ensemble data, we need to provide an additional time corresponding to each data point (set the parameter "series_temporal_order")


# -*- coding: utf-8 -*-
"""
Created on Wed Sep 13 08:37:30 2023

@author: Yinen
"""

import numpy as np
from LK_Info_Flow import multi_causality_est
import scipy.io as sio 
from tqdm import tqdm
NT = 12000


# delay cyclic causal structure
AT = np.zeros((3, 3, 12))
for i in range(3):
    AT[i,i,:]=0.4
AT[1, 1, :] = -0.9
AT[0, 1, 0] = -0.1
AT[0, 1, 1] = -0.2
AT[0, 1, 2] = -0.1

AT[1, 0, 3] = -0.1
AT[1, 0, 4] = -0.2
AT[1, 0, 5] = -0.2

AT[1, 0, 6] = -0.3
AT[1, 0, 7] = -0.4
AT[1, 0, 8] = -0.3

AT[1, 0, 9] = -0.2
AT[1, 0, 10] = -0.2
AT[1, 0, 11] = -0.1

#plot the cyclic element A12 and A21
import matplotlib.pyplot as plt
x1=np.squeeze(AT[0,1,:])
x2=np.squeeze(AT[1,0,:])
plt.plot(x1, marker='o', markeredgecolor='r', markersize=10, linewidth=2, color='r')
plt.plot(x2, marker='o', markeredgecolor='b', markersize=10, linewidth=2, color='b')
plt.legend(['a12', 'a21'])
plt.xticks(range(0,12), [str(i) for i in range(1,12)] + ['0'], fontsize=16)
plt.xlim([-0.5, 11.5])
plt.show()



### repeat the experiment 100 time
B = np.eye(3) * 0.3 + np.ones((3, 3)) * 0.3
X = np.zeros((3, NT+1200))
NIF=[];P=[];SEIF=[];IFs=[];
NIF1=[];P1=[];SEIF1=[];IFs1=[];
for in_ in tqdm(range(1, 101)):

##Exclude the first 1200 time steps and select the stable time seires for information flow calculation.
    vt = 1200
    for it in range(1, NT+vt):
        X[:, it] = np.dot(AT[:, :, (it+2) % 12].T, X[:, it-1]) + np.dot(B, np.random.randn(3, 1)).T

    nn = X.shape
    
    N = X[:, 1201::12].shape
    
    xx = np.zeros((3, N[1]*2))
    t = np.zeros((1, N[1]*2))

# select Feb and Mar
    for i in range(3):
        xx[i, :] = np.reshape(np.array([X[i, 1200::12], X[i, 1201::12]]).T, (1, N[1]*2))
    t = np.reshape(np.array([np.arange(1200,nn[1],12), np.arange(1201,nn[1],12)]).T, (1, N[1]*2))

## information flow for panel data    
    IF=multi_causality_est(X=xx,series_temporal_order=t)
    NIF.append(IF['nIF'])
    P.append(IF['p'])
    SEIF.append(IF['SEIF'])
    IFs.append(IF['IF'])
## information flow for temperal data        
    IF=multi_causality_est(X=xx)
    NIF1.append(IF['nIF'])
    P1.append(IF['p'])
    SEIF1.append(IF['SEIF'])
    IFs1.append(IF['IF'])

# list-> numpy.array
NIF=np.array(NIF);NIF1=np.array(NIF1);
IFs=np.array(IFs);IFs1=np.array(IFs1);
P=np.array(P);P1=np.array(P1);
SEIF=np.array(SEIF);SEIF1=np.array(SEIF1);

print('Ground truth is: X2->X1')
## outputs
#average IF
print('IF(panel)')
mean_IF=np.squeeze(np.average(IFs,axis=0))
print(np.around(mean_IF,2))
print('IF(temperal)')
mean_IF1=np.squeeze(np.average(IFs1,axis=0))
print(np.around(mean_IF1,2))

#average nIF
print('normalized IF(panel)')
mean_NIF=np.squeeze(np.average(NIF,axis=0))
print(np.around(mean_NIF,2))
print('normalized IF(temperal)')
mean_NIF1=np.squeeze(np.average(NIF1,axis=0))
print(np.around(mean_NIF1,2))

#average p-value
print('p-value(panel)')
mean_P=np.squeeze(np.average(P,axis=0))
print(np.around(mean_P,2))
print('p-value(temperal)')
mean_P1=np.squeeze(np.average(P1,axis=0))
print(np.around(mean_P1,2))

#average effect size
print('effect size(panel)')
mean_E=np.average(np.squeeze(np.abs(IFs / SEIF / np.sqrt(1000-3))),axis=0)
print(np.around(mean_E,2))
key=(mean_E[0,1]+mean_E[1,0])/2
print('effect size(temperal)')
mean_E1=np.average(np.squeeze(np.abs(IFs1 / SEIF1 / np.sqrt(2000-3-1))),axis=0)
print(np.around(mean_E1,2))
key1=(mean_E1[0,1]+mean_E1[1,0])/2

A = np.eye(3)
A[0, 1] = 1
sum_val = np.sum(np.abs(A-np.squeeze(np.abs(IFs / SEIF / np.sqrt(1000-3)) >key)), axis=0)
print('structure Hanming Distance(panel)')
print(sum_val)

sum_val = np.sum(np.abs(A-np.squeeze(np.abs(IFs1 / SEIF1 / np.sqrt(2000-3-1)) >key1)), axis=0)
print('structure Hanming Distance(time_series)')
print(sum_val)
    
    
    