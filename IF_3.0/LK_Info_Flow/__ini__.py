"""
Created on Mon Sep  7 16:26:22 2020

@author: easycan @   yinengrong@foxmail.com

The LK_Info_flow module includes implementation of the Liang information flow algorithms.

updates: 
    v2.0  20230524
    According to the data type, Liang information flow(LIF) causal inference
    is divided into three categories: time series, panel data and subsystems
    The first two categories include univariate causality, multivariate 
    causality, normalized causality, and corresponding significance tests

On input:
    xx1,xx2: the series (n by 1 colum vectors)
    xx: the series (n by m colum vectors)  m<=n
    dt: sample time interval
        default: 1 unit time
    np: integer >=1, time advance in performing Euler forward
    differencing, e.g., 1, 2. Unless the series are generated
    with a highly chaotic deterministic system, np=1 should be
    used. 
        default: 1
    st: significance test 
        0 false
        1(default) true
    relative_IF: normalization of IF
        0 false
        1(default) true
    
    
    panel_data_est:
        t: time step order 
    
    
    groups_est
        ind: a 2 colum vector, ind[0] < ind[1] <= m. 
            ind[0]: the index that separates A from the system: 
            (1...ind[0]) forms A,  
    
            ind[1]: the index that separates B from the system: 
            (ind[0]+1...ind[1]) forms B.

On output:
    IF_result(time_series_est and panel_data_est):
        T21: info flow from column to rows, m X m matrix
        tau: normalized info flow
        p:  siginificance(p-value)
        e90:standard error at 90% level
        e95:standard error at 95% level
        e99:standard error at 99% level
    
    IF_result(groups_est):
        TAB: info flow from subspace A to subspace B
        TBA: info flow from subspace B to subspace A
        

    
    
    v1.2  20210709
    add multi_est_Xn in causality_calculate
    add est_panel in causality_calculate
    add est_panel_Xn in causality_calculate
    v1.1  20201125
    add more input cases
    add a version with less warning sentences (causality_calculate)
    v1.0  20200907
    rewrite from causality_est.m tau_est.m multi_causality_est.m 

********************************************************************************************************************
function [T21, err90, err95, err99] = causality_est(x1, x2, np)

Estimate T21, the information transfer from series X2 to series X1 
dt is taken to be 1.

On input:
     X1, X2: the series (n by 1 colum vectors)
     np: integer >=1, time advance in performing Euler forward 
	 differencing, e.g., 1, 2. Unless the series are generated
	 with a highly chaotic deterministic system, np=1 should be
	 used. 

On output:
   T21:  info flow from X2 to X1	(Note: Not X1 -> X2!)
   err90: standard error at 90% significance level
   err95: standard error at 95% significance level
   err99: standard error at 99% significance level

********************************************************************************************************************
[tau21, dH1_star, dH1_noise] = tau_est(xx1, xx2, np)

Estimate tau21, the normalized information transfer from X2 to X1 
On input:
   X1, X2: the series
   np: integer >=1, time advance in performing Euler forward 
   differencing, e.g., 1, 2. Unless the series are generated
   with a highly chaotic deterministic system, np=1 should be
   used.
On output:
   tau21:  relative info flow from X2 to X1	(Note: Not X1 -> X2!)
   dH1_star:  relative contribution from X1 to itself
		(Lyapunov exponent)
   dH1_noise: relative noise contribution
		(noise-to-signal ratio)
********************************************************************************************************************
function [T21, err90, err95, err99] = multi_causality_est(x, np)

Estimate T21, the information transfer from series X2 to series X1, among
M time series (stored as column vectors in X).
dt is taken to be 1.
On input:
    XX: matrix storing the M time series (each as Nx1 column vectors)
       X1 and X2 stored as the first two column vectors
    np: integer >=1, time advance in performing Euler forward 
	 differencing, e.g., 1, 2. Unless the series are generated
	 with a highly chaotic deterministic system, np=1 should be
	 used. 

 On output:
    T21:  info flow from X2 to X1	(Note: Not X1 -> X2!)
    err90: standard error at 90% significance level
    err95: standard error at 95% significance level
    err99: standard error at 99% significance level
*******************************************************************************************************************
Citations: 
     X. San Liang, 2014: Unraveling the cause-effect relation between time series. Phys. Rev. E 90, 052150.
     X. San Liang, 2015: Normalizing the causality between time series. Phys. Rev. E 92, 022126.
     X. San Liang, 2016: Information flow and causality as rigorous notions ab initio. Phys. Rev. E, 94, 052201.
     X. San Liang, 2021: Normalized Multivariate Time Series Causality Analysis and Causal Graph Reconstruction. Entropy. 23. 679.

  for panel data:
     Yineng Rong and X. San Liang, 2021: Panel Data Causal Inference Using a Rigorous Information Flow Analysis for Homogeneous, Independent and Identically Distributed Datasets. IEEE Access. 9, 47266-47274. 
  for subsystems:
     X. San Liang, 2022: The Causal Interaction between Complex Subsystems. Entropy. 24, 3.
    

"""
import scipy
import numpy
from .causal import multi_causality_est_OLS


__version__ = "3.0.0"