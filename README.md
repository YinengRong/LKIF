# LKIF
Created on Mon Sep  7 16:26:22 2020 

@Author: Yineng Rong @ yinengrong@foxmail.com
 (based on the original MATLAB codes by X. San Liang, which are available at www.ncoads.org)
        
@Maintainer: Yineng Rong @ yinengrong@foxmail.com

## About the module
Causal analysis, a fundamental problem in various disciplines, has recently been recognized as a promising approach to developing an explainable deep learning. However, incorporation of causality into artificial intelligence algorithms is faced with challenges such as ambiguity, non-quantifiability, and computational inefficiency in traditional formalisms. Over the past two decades, these challenges have been essentially fixed, with the development of a rigorous and **quantitative causality analysis - the Liang-Kleeman information flow theory** (Liang and Kleeman, 2005; Liang, 2008; 2014; 2016; 2021). This theory, which is based on a rigorous footing and derived from first principles, has resulted in notable scientific discoveries in fields ranging from finance, neuroscience, quantum mechanics, artificial intelligence to oceanography, meteorology and climate science. This module provides a practical implementation of the theory, complete with core codes and selected examples. All the codes are translated from the MATLAB scripts orginally written by X.S. Liang and can be downloaded from http://www.ncoads.org/.


**Citations:**

[1] X.S. Liang, 2016: Information flow and causality as rigorous notions ab initio. Phys. Rev. E, 94, 052202.

[2] X.S. Liang, 2014: Unraveling the cause-effect relation between time series. Phys. Rev. E, 90, 052150.

[3] X.S. Liang, 2021: Normalized multivariate time series causality analysis and causal graph reconstruction. Entropy, 23, 679.


[![License](https://img.shields.io/badge/license-MIT-blue.svg)](https://github.com/YinengRong/LKIF/blob/main/LICENSE)


## Python
### Requirements
* numpy
* scipy
* graphiz (plot causal graph)

### Installation
To install **LK_Info_flow** package, use `pip` as follows:
```sh
  pip install LK_Info_flow
```
or
```sh
pip install .\IF_3.1\dist\LK_Info_flow-py3-none-any.whl
```

### Standard Call for time series or panel data
```sh
from LK_Info_Flow import causal
IF_result=causal.multi_causality_est_OLS(X, max_lag=1, np=1, dt=1, series_temporal_order=None, significance_test=1):
IF_result['IF']
```

**Inputs**:
```sh
   X: matrix storing the M time series (each as Nx1 column vectors)
   
   max_lag: time order of lags (default 1)
   
     >=1 lag time step
   
     -1 Determine the lag order of the data based on the AIC criterion.
   
     -2 Determine the lag order of the data based on the BIC criterion.
     
   np: integer >=1, time advance in performing Euler forward differencing, e.g., 1, 2. Unless the series are generated with a highly chaotic deterministic system, np=1 should be used. 
   
   (default 1)
   
   dt: frequency of sampling 
   
   (default 1)
   
   series_teporal_order: Nx1 column vector that records the timestamps of each sample, with a minimum sampling interval of dt. This option is used for panel data or datasets with missing measurements.
   
   (default [])
   
   significance_test:  1  do the significance test (default)
   
                       0  not (to save the computation time)
```

**Outputs**:
```sh
   IF:               information flow
   
   nIF:              normalized information flow
   
   max_lag:          time order of lags (in IF)
   
   SEIF:             standard error of information flow
   
   err_e90/e95/e99: standard error at 90/95/99# significance level
   
   p:                p-value of information flow
```


### Standard Call for subsystems
```sh
from LK_Info_Flow import causal
IF_result=causal.group_est(X, ind, np=1, dt=1):
IF_result['TAB']
```

**Inputs**:
```sh
   X: matrix storing the M time series (each as Nx1 column vectors)

   ind: 2X1 vector (0<=ind[0]<ind[1]<=M); The series of the components of subsystem A are stored in column [0:ind[0]], and the components of subsystem B are in column [ind[0]:ind[1]]

   np(default 1): integer >=1, time advance in performing Euler forward differencing, e.g., 1, 2. Unless the series are generated with a highly chaotic deterministic system, np=1 should be used. 
   
   dt(default 1): frequency of sampling 
```

   
**outputs**:
```sh
    TAB: info flow from subspace A to subspace B
    TBA: info flow from subspace B to subspace A
```


**More details are in the example file (https://github.com/YinengRong/LKIF/blob/main/IF_3.0/LK_Info_Flow/examples/example.ipynb)**

There are 8 cases in the file:

1 bivariate causality

2 multivariable causality

3 causality analysis with panel data, discontinuous time series or ensemble data

4 causality between subsystems

5 time delay, confounder, cyclic causality

6 computational cost for large-scale Liang information flow analysis

7 causality analysis with data in the presencee of cross-correlated noise

8 time varying causality


## R
### Requirements
The library is developed in version 4.4.1 of R
* stats
* MASS
* igraph (plot causal graph)

### Standard Call for time series or panel data
```sh
source(LK_Info_Flow.R)
IF_result=multi_causality_est_OLS(X, max_lag=1, np=1, dt=1, series_temporal_order=NULL, significance_test=1):
IF_result$IF
```
=>inputs and outputs can reffer to python package
### Standard Call for subsystems
```sh
source(LK_Info_Flow.R)
IF_result=group_est(X, ind, np=1, dt=1):
IF_result$TAB
```
=>inputs and outputs can reffer to python package

## Matlab
### Requirements
The library is developed in version 2020a of matlab

### Standard Call for time series or panel data
```sh
addpath ...
IF_result=multi_causality_est_OLS(X, max_lag, np, dt, series_temporal_order, significance_test):
IF_result.IF
```
=>inputs and outputs can reffer to python package
### Standard Call for subsystems
```sh
addpath ...
IF_result=group_est(X, ind, np, dt):
IF_result.TAB
```
=>inputs and outputs can reffer to python package

## Citations:
* **X. San Liang**, 2014: Unraveling the cause-effect relation between time series. Phys. Rev. E 90, 052150.
* **X. San Liang**, 2015: Normalizing the causality between time series. Phys. Rev. E 92, 022126.
* **X. San Liang**, 2016: Information flow and causality as rigorous notions ab initio. Phys. Rev. E, 94, 052201.
* **X. San Liang**, 2021: Normalized Multivariate Time Series Causality Analysis and Causal Graph Reconstruction. Entropy. 23. 679.


**for panel data:**

* **Yineng Rong** and **X. San Liang**, 2021: Panel Data Causal Inference Using a Rigorous Information Flow Analysis for Homogeneous, Independent and Identically Distributed Datasets. IEEE Access. 9, 47266-47274.
  
  
**for subsystems:**
* **X. San Liang, 2022**: The Causal Interaction between Complex Subsystems. Entropy. 24, 3.

## Chagelogs in V3.1
Updated the formula for p-value
Added toolkit for R language


The codes are rewritten from the MATLAB scripts (causality_est.m, tau_est.m, multi_causality_est.m, etc.) on Prof. X. San Liang's lab website (http://www.ncoads.org/article/show/67.aspx)
