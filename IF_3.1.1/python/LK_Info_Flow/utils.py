# -*- coding: utf-8 -*-
"""
Created on Tue Jul 23 17:01:27 2024

@author: Yinen
"""

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
        max_lag = npy.where(aic == npy.min(aic))[0][0]+1
    else:
        max_lag = npy.where(bic == npy.min(bic))[0][0]+1
    return max_lag

def temporal_lag(X, m1, max_lag):
    # Removing causality from future to past when max_lag > 1
    return X[:m1].reshape([m1,m1,max_lag])


