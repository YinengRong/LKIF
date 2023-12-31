# -*- coding: utf-8 -*-
"""
Created on Wed Jul 12 10:45:37 2023

@author: Yinen

Rong Yineng (yinengrong@foxmail.com)
see https://github.com/YinengRong/LKIF for details and examples

On input:
   X: matrix storing the M time series (each as Nx1 column vectors)
   max_lag: time order of lags (default 1)
     >=1 lag time step
     -1 Determine the lag order of the data based on the AIC criterion.
     -2 Determine the lag order of the data based on the BIC criterion.
   np: integer >=1, time advance in performing Euler forward 
     differencing, e.g., 1, 2. Unless the series are generated
     with a highly chaotic deterministic system, np=1 should be used. 
     (default 1)
   dt: frequency of sampling 
     (default 1)
   series_teporal_order: Nx1 column vector that records the timestamps 
     of each sample, with a minimum sampling interval of dt. This option
     is used for panel data or datasets with missing measurements..
     (default [])
   significance_test:  1  do the significance test (default)
                       0  not (to save the computation time)
On output:
   IF:               information flow
   nIF:              normalized information flow
   max_lag:          time order of lags (in IF)
   SEIF:             standard error of information flow
   err_e90/e95/e99: standard error at 90/95/99# significance level
   dnoise:           dH1_noise/dt with correlated noise
   p:                p-value of information flow

Citations: 
   X.S. Liang, 2016: Information flow and causality as rigorous notions ab initio. Phys. Rev. E, 94, 052201.
   X.S. Liang, 2014: Unraveling the cause-effect relation between time series. Phys. Rev. E 90, 052150.
   X.S. Liang, 2015: Normalizing the causality between time series. Phys. Rev. E 92, 022126.
   X.S. Liang, 2021: Normalized Multivariate Time Series Causality Analysis and Causal Graph Reconstruction. Entropy. 23. 679.
   X.S. Liang, Chen, D. K., Zhang, R. H. 2023: Quantitative causality, causality-aided discovery, and causal machine learning. Ocean-Land-Atmos Res.
"""

import numpy as npy
from scipy.stats import norm

def pre_panel_data(XX, t, q):
    n, _, m = XX.shape
    if len(t) == 0:
        t = npy.arange(m)
    tt = npy.full((q+2, m+q), -npy.pi**2)
    for k in range(q+1):
        tt[k, k:k+m] = t
    tt[q+1, :] = tt[q, :] - 1
    ind = npy.where(npy.max(npy.abs(npy.diff(tt,axis=0)) - 1,axis=0)==0)[0]
    q1 = q + 1
    M = len(ind)
    nq = n * q
    X0 = XX[:, 0, ind].reshape(n, M)
    XL = XX[:, 1:q1, ind].transpose(1,0,2).reshape(nq, M)
    return X0, XL

def infocrit(L, k, m):
    if m - k - 1 <= 0:
        aic = npy.nan
    else:
        aic = -2 * L + 2 * k * (m / (m - k - 1))
    bic = -2 * L + k * npy.log(m)
    return aic, bic

def lag_order(X, t, option):
    #  calculates the lag order of the data
    n, m = X.shape
    X = X - npy.mean(X, axis=1, keepdims=True)
    morder = npy.arange(1, min([max([npy.floor(m/(n**2+n)), 1]), 20]) + 1)
    nummo = len(morder)
    aic = npy.full(nummo, npy.nan)
    bic = npy.full(nummo, npy.nan)
    q = max(morder)
    q1 =q+1
    XX = npy.zeros((n, q1, m+q))
    for k in range(q1):
        XX[:, k, k:k+m] = X
    for i in range(nummo):
        q = int(morder[i])
        if q >= m:
            print('WARNING: model order too large (must be < %d)' % m)
            continue
        XX = XX[:, :, :m]
        X0, XL = pre_panel_data(XX, t, q)
        M = X0.shape[1]
        A = npy.linalg.lstsq(XL.T, X0.T, rcond=None)[0].T
        E = X0 - A @ XL
        DSIG = npy.linalg.det(E @ E.T / (M-1))
        if DSIG <= 0:
            print('WARNING: residuals covariance not positive definite')
            continue
        aic[i], bic[i] = infocrit(-(M/2) * npy.log(DSIG), q*n*n, M)
    if option == -1:
        max_lag = npy.where(aic == npy.min(aic))[0][0]
    else:
        max_lag = npy.where(bic == npy.min(bic))[0][0]
    return max_lag

def temporal_lag(X, m1, max_lag):
    # Removing causality from future to past when max_lag > 1
    return X[:m1].reshape([m1,m1,max_lag])

def multi_causality_est_OLS(X, max_lag=1, np=1, dt=1, series_temporal_order=None, significance_test=1):
    
    
    ## data pre-processing
    m, n = X.shape
    if series_temporal_order is None:
        series_temporal_order = npy.arange(1, dt * n + 1, dt)
    
    # max_lag
    if max_lag < 0:
        max_lag = lag_order(X, series_temporal_order/dt, max_lag)
    # panel data
    q1 = max_lag + 1
    XX = npy.zeros((m, q1, n+q1-1))
    for k in range(q1):
        XX[:, k, k:k+n] = X
    XX = XX[:, :, :n]
    XL, XR = pre_panel_data(XX, series_temporal_order/dt, max_lag)
    m2, n2 = XL.shape
    X1 = npy.vstack((XR, npy.ones((1, n2))))
    X2 = npy.vstack((XL, XR[:m2*(max_lag-1), :], npy.ones((1, n2))))
    
    
    ## calculate information flow
    # information flow(IF)
    A0 = X2*npy.mat(X1).I;#X2 @ npy.linalg.lstsq(X1.T, X1.T, rcond=None)[0]
    A = (A0 - npy.eye(m2*max_lag+1)) / dt
    At = A[:-1, :-1]
    C = npy.cov(X1[:-1, :])
    IF = ((C / npy.diag(C)).T)*npy.array(At)
    
    # normalized IF
    E = X2 - A0 @ X1
    SIG = E @ E.T / (n - np-m2-1)
    B = npy.sqrt(npy.abs(npy.diag(SIG[:-1, :-1])) / dt)
    dH1_noise = npy.sum(npy.array(B)**2, axis=0) / (2 * npy.diag(C))
    Z=(npy.sum(npy.abs(IF), axis=1) + npy.abs(dH1_noise));
    Z = npy.repeat(Z, m2*max_lag).reshape(m2*max_lag,m2*max_lag)
    nIF = IF / Z;
    
    ## output 
    IF_result = {}
    IF_result['IF'] = temporal_lag(IF, m2, max_lag) 
    IF_result['nIF'] = temporal_lag(nIF, m2, max_lag)
    IF_result['dnoise'] = dH1_noise[-m2:]
    IF_result['max_lag'] = max_lag
    ##significance test
    if significance_test == 1:
        se_a = npy.sqrt(npy.diag(SIG[:-1, :-1]).reshape(1,m2*max_lag)* 
                        npy.diag(npy.linalg.inv(X1[:-1, :] @ X1[:-1, :].T)).reshape(m2*max_lag,1)).T
        SE_IF = se_a * npy.abs(C / npy.diag(C)).T / dt # standarized error of IF
        p = 1 - (1 - norm.cdf(npy.abs(IF / SE_IF))) * 2 # p-value
        z99 = 2.56;z95 = 1.96;z90 = 1.65
        e90_IF = SE_IF * z90
        e95_IF = SE_IF * z95
        e99_IF = SE_IF * z99
        IF_result['SEIF'] = temporal_lag(SE_IF, m2, max_lag)
        IF_result['err_e90'] = temporal_lag(e90_IF, m2, max_lag)
        IF_result['err_e95'] = temporal_lag(e95_IF, m2, max_lag)
        IF_result['err_e99'] = temporal_lag(e99_IF, m2, max_lag)
        IF_result['p'] = temporal_lag(p, m2, max_lag)
        
    return IF_result



def groups_est(xx,ind,dt=1,np=1):
    if  npy.array(npy.array(xx.shape).shape)!=2:
        # print('multi_est need inputs like: if with x1,x2:  [ x1[:,] x2[:,m] ], [ x1[:,m] x2[:,] ], [ x1[:,m] x2[:,n] ], [ x1[:,] x2[:,] ]')
        # print(' if with x only: x[:,0]->x1 x[:,1]->x2 x[:,1:]->others ')
        return 0
    N=npy.array(xx.shape)
    nm = N[0]; M = N[1];
    
    k=0;
    
    
    dx = npy.zeros([nm-np,M]);
    for k in range(M):
        dx[:,k] = (xx[np:,k]- xx[:nm-np,k] )/(np*dt)
    x = xx[:nm-np, :];
    del xx
    
    NL,_=npy.array(x.shape)
    
    
    
    C=npy.cov(x.transpose());  #
    

        
    dC=npy.zeros([M,M])
    for kk in range(M):
        for k in range(M):
            dC[k,kk] = sum((x[:,k] - npy.mean(x[:,k])) * (dx[:,kk] - npy.mean(dx[:,kk]))); 
    dC=dC/(NL-1)

    try: 
        ann = npy.dot(npy.linalg.inv(C),dC)
    except:
        ann = npy.ones(shape=npy.shape(dC))*float('inf')
        # print ('x1==x2, please check the data')
    
    Cr=npy.zeros([M,M])
    Crr=Cr+0;
    Crr[:ind[0],:ind[0]]=C[:ind[0],:ind[0]];
    for i in range(ind[0],ind[1]):
        Crr[i,i]=1
    
    Cr[ind[0]:ind[1],ind[0]:ind[1]]=C[ind[0]:ind[1],ind[0]:ind[1]];
    for i in range(ind[0]):
        Cr[i,i]=1
    invCr=npy.linalg.inv(Cr);
    invCrr=npy.linalg.inv(Crr);
    AC=npy.dot(ann[:ind[1],:ind[1]],C[:ind[1],:ind[1]])
    TBA=npy.trace(npy.dot(invCr[ind[0]:ind[1],ind[0]:ind[1]],AC[ind[0]:ind[1],ind[0]:ind[1]].transpose()))-npy.trace(ann[ind[0]:ind[1],ind[0]:ind[1]])
    
    TAB=npy.trace(npy.dot(invCrr[:ind[0],:ind[0]],AC[:ind[0],:ind[0]].transpose()))-npy.trace(ann[:ind[0],:ind[0]])
    
    
    

    IF_result={'TAB':TAB,'TBA':TBA}
    return IF_result





