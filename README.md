# LKIF
Created on Mon Sep  7 16:26:22 2020 

@Author: Yineng Rong @ yinengrong@foxmail.com
 (based on the original MATLAB codes by X. San Liang, which are available at www.ncoads.org)
        
@Maintainer: Yineng Rong @ yinengrong@foxmail.com

## About the module
Causal analysis, a fundamental problem in various disciplines, has recently been recognized as a promising approach to developing an explainable deep learning. However, incorporation of causality into artificial intelligence algorithms is faced with challenges such as ambiguity, non-quantifiability, and computational inefficiency in traditional formalisms. Over the past two decades, these challenges have been essentially fixed, with the development of a rigorous and **quantitative causality analysis - the Liang-Kleeman information flow theory** (Liang and Kleeman, 2005; Liang, 2008; 2014; 2016; 2021). This theory, which is based on a rigorous footing and derived from first principles, has resulted in notable scientific discoveries in fields ranging from finance, neuroscience, quantum mechanics, artificial intelligence to oceanography, meteorology and climate science. This module provides a practical implementation of the theory, complete with core codes and selected examples. All the codes are translated from the MATLAB scripts orginally written by X.S. Liang and can be downloaded from http://www.ncoads.org/.




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
pip install .\codes\python\dist\LK_Info_flow-py3-none-any.whl
```

### Standard call for time series or panel data
```sh
from LK_Info_Flow import multi_causality_est
IF_result=multi_causality_est(X, max_lag=1, np=1, dt=1, series_temporal_order=None, significance_test=1):
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
   a structure value IF_result with sub

   IF:               information flow
   
   nIF:              normalized information flow
   
   SEIF:             standard error of information flow
   
   err_e90/e95/e99:  standard error at 90/95/99% confidence level
   
   p:                p-value of information flow
```


### Standard call for subsystems
```sh
from LK_Info_Flow import causality_subspace
IF_result=causality_subspace(X, ind, np=1, dt=1):
IF_result['TAB']
```

**Inputs**:
```sh
   X: matrix storing the M time series (each as Nx1 column vectors)

   ind: 2X1 vector (0<=ind[0]<ind[1]<=M); the index ind[0] that separates A from the system: [0:ind[0]] forms A, and the index ind[0] together with ind[1] separates B from the system: [ind[0]:ind[1]] forms B.

   np(default 1): integer >=1, time advance in performing Euler forward differencing, e.g., 1, 2. Unless the series are generated with a highly chaotic deterministic system, np=1 should be used. 
   
   dt(default 1): frequency of sampling 
```

   
**outputs**:
```sh
    TAB: info flow from subspace A to subspace B
    TBA: info flow from subspace B to subspace A
```


**More details are in the example file ([example.ipynb](./codes/python/LK_Info_Flow/examples/example.ipynb))**

There are 8 cases in the file:

1. Bivariate causality analysis (Liang, 2014);

2. Multivariable causality analysis (Liang, 2021);

3. Causality analysis with panel data, discontinuous time series or ensemble data (Rong and Liang, 2021)

4. Causal inference between different subsystems (Liang, 2022);

5. Takens' theorem

6. Computational cost for large-scale Liang information flow analysis

7. Causality analysis with data in the presencee of cross-correlated noise

8. Time varying causality analysis;


## R
### Requirements
The library is developed in version 4.4.1 of R
* stats
* MASS
* igraph (Optional, plot causal graph)
* R.matlab (Optional, to achieve data exchange between different software)

### Standard call for time series or panel data
```sh
source(LK_Info_Flow.R)
IF_result=multi_causality_est(X, max_lag=1, np=1, dt=1, series_temporal_order=NULL, significance_test=1):
IF_result$IF
```
=>inputs and outputs can reffer to python package
### Standard call for subsystems
```sh
source(LK_Info_Flow.R)
IF_result=causality_subspace(X, ind, np=1, dt=1):
IF_result$TAB
```
=>inputs and outputs can reffer to python package

## Matlab
### Requirements
The library is developed in version 2020a of matlab

### Standard call for time series or panel data
```sh
addpath ...
IF_result=multi_causality_est(X, np, dt,  max_lag, series_temporal_order, significance_test):
IF_result.IF
```
=>inputs and outputs can reffer to python package
### Standard call for subsystems
```sh
addpath ...
IF_result=causality_subspace(X, ind, np, dt):
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


The codes are rewritten from the MATLAB scripts (causality_est.m, tau_est.m, multi_causality_est.m, etc.) on Prof. X. San Liang's lab website (http://www.ncoads.org/article/show/67.aspx)
