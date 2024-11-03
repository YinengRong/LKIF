# -*- coding: utf-8 -*-
"""
Created on Tue Jul 23 17:01:48 2024

@author: Yineng Rong (yinengrong@foxmail.com)

% 
% function IF = causality_subspace(x, r, s, np)
%
% Infer the IF between two subspaces A and B
% (M time series stored as column vectors in X).
% Units: Nats per unit time 
% (dt is taken to be 1 here; the final result should be divided by dt).
%
% On input:
%    X: matrix storing the M time series (each as NLx1 column vectors),
%       the first r series forming subspace A, 
%	the next s-r series forming subspace B.
%
%    r: the index that separates A from the system: 
%	(1...r) forms A,  r < s <= M.
%
%    s: the index together with r separates B from the system:
%	 (r+1,...,s) forms B,  s>r, s<=M.
%
%    np: integer >=1, time advance in performing Euler forward 
%	 differencing, e.g., 1, 2. Unless the series are generated
%	 with a highly chaotic deterministic system, np=1 should be
%	 used. 
%
% On output IF[XXXX]:
%    TAB:  info flow from subspace A to subspace B
%    TBA:  info flow from subspace B to subspace A
%
% Citations: 
% Liang, X.S. The causal interaction between complex subsystems. 
%             Entropy, 2022, 24, 3.
%
% Note: Here the symbols M and NL correspond respectively 
%         to the symbols N and K in the paper.
%
% Orignal matlab codes by X. San Liang in 2022
% Authorized by X.S. Liang to release here.




"""
from LK_Info_Flow.utils import *
import numpy as npy
def causality_subspace(xx,ind,dt=1,np=1):
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
    Cr=Cr[:ind[1],:ind[1]]
    Crr=Crr[:ind[1],:ind[1]]
    invCr=npy.linalg.inv(Cr);
    invCrr=npy.linalg.inv(Crr);
    AC=npy.dot(ann[:ind[1],:ind[1]].T,C[:ind[1],:ind[1]])
    TBA=npy.trace(npy.dot(invCr[ind[0]:ind[1],ind[0]:ind[1]],AC[ind[0]:ind[1],ind[0]:ind[1]].transpose()))-npy.trace(ann[ind[0]:ind[1],ind[0]:ind[1]])
    
    TAB=npy.trace(npy.dot(invCrr[:ind[0],:ind[0]],AC[:ind[0],:ind[0]].transpose()))-npy.trace(ann[:ind[0],:ind[0]])
    
    
    

    IF_result={'TAB':TAB,'TBA':TBA}
    return IF_result




