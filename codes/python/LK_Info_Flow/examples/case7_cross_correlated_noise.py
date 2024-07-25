# case 7 (data with cross-correlated noise):

# -*- coding: utf-8 -*-
"""
Created on Sun Dec 31 15:03:43 2023
To generate the cross-correlated noise, we use the functions (var_to_tsdata, var_specrad, var_decay, genvar) which are all rewritten based on the Matlab script of MVGC matlab toolbox by Barnett and Seth:
Barnett, L., and A. K. Seth, 2014: The MVGC multivariate Granger causality toolbox: A new approach to Granger-causal inference. J. Neurosci. Methods, 223, 50–68, https://doi.org/10.1016/j.jneumeth.2013.10.018
@author: Yinen
"""

import numpy as np
from LK_Info_Flow import multi_causality_est
from scipy.linalg import eig  
from scipy.linalg import cholesky  
import scipy.io as scio

def var_to_tsdata(A, SIG, m, N=1, mtrunc=None, decayfac=100):  
    if len(A.shape) == 1:  
        n = A.shape[0]  
        A = np.array([A])  
    else:  
        n = A.shape[1]  
  
    if mtrunc is None:  
        if decayfac is None:  
            decayfac = 100  
        rho = var_specrad(A)  
        assert rho < 1, 'unstable VAR'
        mtrunc = int((np.log(np.finfo(float).eps) - decayfac) / np.log(rho))  # enough time for autocovariance to decay to fp accuracy (and then some)  
    else:  
        assert np.isscalar(mtrunc) and np.isinteger(mtrunc) and mtrunc >= 0, 'mtrunc parameter must be a non-negative integer'
    try:
        C = cholesky(SIG, lower=True)  
    except:
        print( 'covariance matrix not positive-definite')  

    if N > 1:  # multi-trial  
        X = np.zeros((n, m, N))  
        E = np.zeros((n, m, N))  
        for r in range(N):  
            X[:, :, r], E[:, :, r] = genvar(A, C @ np.random.randn(n, m + mtrunc), mtrunc)  
        
    else:  # single trial  
        X, E = genvar(A, C @ np.random.randn(n, m + mtrunc), mtrunc)  
    return X, E, mtrunc

  
def var_specrad(A, newrho=None):  
    A_shape = A.shape  
    if len(A_shape)==2:
        A=A.reshape([A_shape[0],A_shape[1],1])
    n, n1, p = A.shape  
    assert n1 == n, 'VAR coefficients matrix has bad shape'  
    pn1 = (p-1)*n  
    if pn1!=0:
    # construct VAR coefficients for 1-lag problem  
        A1 = np.concatenate((A.reshape(n, p*n), np.eye(pn1), np.zeros((pn1, n))), axis=1)  
    else:
        A1 = A.reshape(n, p*n)
  
    # calculate spectral radius  
    eig_A,_ = eig(A1)
    rho = max(np.abs(eig_A))  
  
    if newrho is None or len(newrho) == 0:  
        #assert len(rho) <= 1, 'too many output parameters'  
        return rho  
    else:  
        return var_decay(A, newrho/rho), rho  # adjusted coefficients, previous value of spectral radius

def var_decay(A, dfac):  
    p = A.shape[2]  # Assuming A is a 3D numpy array  
    f = dfac  
    for k in range(p):  
        A[:, :, k] *= f  # Multiply each slice of A by f  
        f *= dfac  # Update f  
    return A  # Return modified A

def genvar(A, E, trunc=0):  
    assert np.isrealobj(A) and (A.ndim == 2 or A.ndim == 3), 'VAR coefficients must be a 2D or 3D matrix'  
    assert np.isrealobj(E) and E.ndim == 2, 'Residuals covariance must be a row vector or 2D matrix'  
    if A.ndim == 2:
        n1, n2 = A.shape
        A=A.reshape([n1,n2,1])
    n1, n2, p = A.shape  
    n, m = E.shape  
      
    assert n1 == n2, 'VAR coefficients blocks not square'  
    assert n1 == n,  'Residuals covariance matrix doesn''t match VAR coefficients matrix'  
    assert trunc >= 0 & trunc < m,'bad truncation'
      
    X = E.copy()  # Initialize X as a copy of E  
      
    for t in range(p+1, m):  
        for k in range(1, p+1):  
            X[:, t] = X[:, t] + np.dot(A[:, :, k-1], X[:, t-k])  # Update X using the VAR coefficients  
    if trunc >= 0 and trunc < m:  
        X = X[:, trunc+1:]  # Truncate X  
        if len(E.shape) == 2:  
            E = E[:, trunc+1:]  # Truncate E if it's a matrix  
            
    return X, E  # Return X and E (if there are more than one output argument)



if __name__=='__main__':
    NT = 100000  
    # Case 1  
    AT1 = [[0.5, -0.5], [-0.0, 0.6]]  
    
    # Case 2  
    AT2 = [[-0.5, 0.5], [-0.2, 0.5]]  
    
    # Case 3  
    AT3 = [[0.5, -0.2], [-0.5, 0.25]]  
    
    # Case 4  
    AT4 = [[0.25, -0.1], [-0.2, 0.1]]

    AT=np.array(AT1)



    j = 1  
    CE_3d=[];CX_3d=[];NIF_3d=[];IFs_3d=[];p_3d=[];SEIF_3d=[]
    for i in np.arange(-0.5, 0.5, 0.01):#tqdm(np.arange(-0.5, 0.5, 0.01)):  
        #print('%.2f'%i)  
        SIFT = np.array([[0.5, i], [i, 0.5]])    
        X, E, _ = var_to_tsdata(AT, SIFT, NT)
      
        CE_3d.append(np.corrcoef(E))
        CX_3d.append(np.corrcoef(X))
      
        cau = multi_causality_est(X)  
        NIF_3d.append(cau['nIF'])
        p_3d.append(cau['p'])
        SEIF_3d.append(cau['SEIF'])
        IFs_3d.append(cau['IF'])
      
        j += 1

    scio.savemat('data_1.mat',{'CX': CX_3d, 'CE': CE_3d,'NIF': NIF_3d})


## plot 
import matplotlib.pyplot as plt
fig=plt.figure(figsize=(12,6), dpi=100)

text_xy=[-0.95, 0.7]
texts=['a)','b)','c)','d)']
ylim=[[10**(-7),10**0.5],[10**(-4),10**0.2],[10**(-4),10**0.2],[10**(-4),10**0.2]]

# normalized IF
for i in range(4):
    data1=scio.loadmat('data_'+str(i+1)+'.mat')
    C=data1['CE']
    x=C[:,0,1]
    NIF=data1['NIF']
    y1=abs(NIF[:,0,1])
    y2=abs(NIF[:,1,0])
    plt.subplot(221+i)
    #plot
    plt.semilogy(x, y1, lw=1.5, c='b', alpha=0.7)
    plt.semilogy(x, y2, lw=1.5, c='r', alpha=0.7)
    plt.legend([r'$\tau_{21}$',r'$\tau_{12}$'],loc='lower left')
    plt.xlim([-1,1])
    plt.ylim(ylim[i])
    plt.text(text_xy[0],text_xy[1],texts[i],fontdict=dict(fontsize=14, color='k',weight='bold'))


# effect size
fig=plt.figure(figsize=(12,6), dpi=100)

text_y=[0.00003,0.015,0.015,0.015]
ylim=[[1e-5,1],[1e-2,1],[1e-2,1],[1e-2,1]]
for i in range(4):
    data1=scio.loadmat('data'+str(i+1)+'.mat')
    C=data1['CE']
    x=C[0,1]
    dIF=abs(data1['IF']/data1['SEIF']/np.sqrt(NT))  #effect size
    y1=abs(dIF[0,1])
    y2=abs(dIF[1,0])
    
    plt.subplot(221+i)
    plt.semilogy(x, y1, lw=1.5, c='b', alpha=0.7)
    plt.semilogy(x, y2, lw=1.5, c='r', alpha=0.7)
    plt.text(0.8,text_y[i],texts[i],         fontdict=dict(fontsize=16, color='k',
                           #family='monospace',#fontsize 'serif', 'sans-serif', 'cursive', 'fantasy', 'monospace'
                           weight='bold',#'light', 'normal', 'medium', 'semibold', 'bold', 'heavy', 'black'
                          
                          )
    )
    plt.legend([r'$X_2->X_1$',r'$X_1->X_2$'])
    plt.ylim(ylim[i])
    plt.xlim([-0.95,0.95])



# relative importance of IF(X2->X1)/IF(X1->X2)
fig=plt.figure(figsize=(9,3), dpi=100)
color=['b','y','g','r']
for i in range(1,4):
    data1=scio.loadmat('data'+str(i+1)+'.mat')
    C=data1['CE']
    x=C[0,1]
    NIF=data1['NIF']
    y1=abs(NIF[0,1])
    y2=abs(NIF[1,0])
    #plt.subplot(221+i)
    #绘图命令
    plt.plot(x, y1/y2, lw=2, c=color[i], alpha=0.7)


plt.legend([r'Case 2',r'Case 3',r'Case 4'])
plt.plot([-1,1], [1 ,1], lw=4, c='k', alpha=0.7)
plt.xlim([-1,1])
plt.ylim([-0.2,1.5])
