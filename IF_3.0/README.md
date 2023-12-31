Created on Mon Sep  7 16:26:22 2020

@author: Yinengrong @   yinengrong@foxmail.com

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
   err_e90/e95/e99:  standard error at 90/95/99# significance level
   p:                p-value of information flow

Citations: 
   X.S. Liang, 2016: Information flow and causality as rigorous notions ab initio. Phys. Rev. E, 94, 052201.
   X.S. Liang, 2014: Unraveling the cause-effect relation between time series. Phys. Rev. E 90, 052150.
   X.S. Liang, 2015: Normalizing the causality between time series. Phys. Rev. E 92, 022126.
   X.S. Liang, 2021: Normalized multivariate time series causality analysis and causal graph reconstruction. Entropy. 23. 679.
   X.S. Liang, Chen, D. K., Zhang, R. H. 2023: Quantitative causality, causality-aided discovery, and causal machine learning. Ocean-Land-Atmos Res.


