# -*- coding: utf-8 -*-
"""
Created on Tue Jul 23 17:00:15 2024

@author: Yinen


multi_causality_est(X, max_lag=1, np=1, dt=1, series_temporal_order=None, significance_test=1)

% By Yineng Rong (yinengrong@foxmail.com)
% 
% On input:
%    X: matrix storing the M time series (each as Nx1 column vectors)
%    max_lag: time order of lags (default 1)
%    >=1 lag time step
%    -1 Determine the lag order of the data based on the AIC criterion.
%   - 2 Determine the lag order of the data based on the BIC criterion.
%    np: integer >=1, time advance in performing Euler forward 
%	 differencing, e.g., 1, 2. Unless the series are generated
%	 with a highly chaotic deterministic system, np=1 should be
%	 used. 
%    (default 1)
%    dt: frequency of sampling (default 1)
%
%    series_teporal_order: Nx1 column vectors, records the timestamp of 
%    each sample, with a minimum sampling interval of dt (used for 
%    panel data, or missing measurement data). =>case3
%    (default [])
%    significance_test:  1  do the significance test (default)
%           =0  not (to save the computation time)
% On output:
%    a structure value IF_result with sub
%    IF:  information flow
%    nIF: normalized information flow
%    SEIF: standard error of information flow
%    errr.e90/e95/e99:standard error at 90/95/99% confidence level
%    dnoise: dH1_noise/dt with correlated noise
%    dnoise_old: dH1_noise/dt without correlated noise (##  been commented out.)
%    nIF: normalized information flow without correlated noise
%    p: p-value of information flow
%
% Citations: 
%    X.S. Liang, 2016: Information flow and causality as rigorous notions
%    ab initio. Phys. Rev. E, 94, 052201.
%    X.S. Liang, 2014: Unraveling the cause-effect relation between time series. Phys. Rev. E 90, 052150.
%    X.S. Liang, 2015: Normalizing the causality between time series. Phys. Rev. E 92, 022126.
%    X.S. Liang, 2021: Normalized Multivariate Time Series Causality Analysis and Causal Graph Reconstruction. Entropy. 23. 679.

% Note: This is an alternative and simplified version of
%   multi_causality_est.m, which was originally coded by X.S. Liang
%   Here all the causal relations are inferred once for all.
% 

"""

from LK_Info_Flow.utils import *
import numpy as npy
def multi_causality_est(X, max_lag=1, np=1, dt=1, series_temporal_order=None, significance_test=1):
    
    
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
    #B = npy.sqrt(npy.abs(npy.diag(SIG[:-1, :-1])) / dt)  #Liang,2014 2021
    B = npy.sqrt(npy.abs(SIG[:-1, :-1]) / dt)
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
        p = (1 - norm.cdf(npy.abs(IF / SE_IF))) * 2 # p-value
        z99 = 2.56;z95 = 1.96;z90 = 1.65 # confidence level
        e90_IF = SE_IF * z90
        e95_IF = SE_IF * z95
        e99_IF = SE_IF * z99
        IF_result['SEIF'] = temporal_lag(SE_IF, m2, max_lag)
        IF_result['err_e90'] = temporal_lag(e90_IF, m2, max_lag)
        IF_result['err_e95'] = temporal_lag(e95_IF, m2, max_lag)
        IF_result['err_e99'] = temporal_lag(e99_IF, m2, max_lag)
        IF_result['p'] = temporal_lag(p, m2, max_lag)
        
    return IF_result