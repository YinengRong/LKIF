function [tau21, dH1_star, dH1_noise] = multi_tau_est(xx, np)
% 
% function [tau21, dH1_star, dH1_noise] = multi_tau_est(x, np)
%
% Estimate tau21, normalized info flow from series X2 to series X1.
% dt is taken to be 1.
%
% On input:
%    XX: matrix storing the M time series (each as Nx1 column vectors)
%       X1 and X2 stored as the first two column vectors
%    np: integer >=1, time advance in performing Euler forward 
%	 differencing, e.g., 1, 2. Unless the series are generated
%	 with a highly chaotic deterministic system, np=1 should be
%	 used. 
%
% On output:
%    tau21:  relative info flow from X2 to X1	(Note: Not X1 -> X2!)
%    dH1_star:  normalized Lyapunov exponent of X1
%    dH1_noise: relative noise contribution (noise-to-signal ratio)
%
% Citations: 
%    X.S. Liang, 2021: Normalized multivariate time series causality
%    analysis and causal graph reconstruction. Entropy, 23, 679.
%    X.S. Liang, 2016: Information flow and causality as rigorous notions
%    ab initio. Phys. Rev. E, 94, 052201.
%    X.S. Liang, 2014: Unraveling the cause-effect relation between time series. Phys. Rev. E 90, 052150.
%    X.S. Liang, 2015: Normalizing the causality between time series. Phys. Rev. E 92, 022126.

% Note: This is the original matlab program for the normalization of 
%       multivariate causalities coded by X. San Liang in 2021. 
% Authorized by X.S. Liang to release here.


dt = 1;	


[nm, M] = size(xx);		% nm: length of series
				% M:  dim of system (or # of series)


dx1(:,1) = (xx(1+np:nm, 1) - xx(1:nm-np, 1)) / (np*dt);
x(:,1:M) = xx(1:nm-np, 1:M);


clear xx;

N = nm-np;




C = cov(x);		% MxM matrix



for k = 1 : M,
%dC(1,1) = sum((x(:,1) - mean(x(:,1))) .* (dx1 - mean(dx1))); 
dC(k,1) = sum((x(:,k) - mean(x(:,k))) .* (dx1 - mean(dx1))); 
end
dC = dC / (N-1);


a1n = inv(C) * dC;

%% The vector a1n stores coefficients a11, a12, a13,...,a1M
%% particularly, a12 = a1n(2,1), and hence T21 = C12/C11 * a12
%% but here they are not needed
% a12 = a1n(2,1);
% T21 = C(1,2)/C(1,1) * a12;
% T11 = a1n(1,1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% For convenience, let T1 stores the info flow from Xj --> X1
%     j = 1, ..., M, (where T1(1)=T11=dH1*/dt,  T1(2) = T21, ...
%
  for j = 1:M,
      T1(j,1) = C(1,j) / C(1,1)  *  a1n(j,1);
  end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

f1 = mean(dx1);
for k = 1 : M,
    f1 = f1 - a1n(k,1) * mean(x(:,k));
end

R1 = dx1 - f1;
for k = 1 : M,
    R1 = R1 - a1n(k,1) * x(:,k);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% noise-to-signal ratio %%%
   Q1 = sum(R1 .* R1);
   g11 = Q1 * dt / N;
dH1_noise = g11 / (2*C(1,1));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


   Z = sum(abs(T1)) + abs(dH1_noise);
tau21 = T1(2) / Z;
dH1_star = T1(1) / Z;
dH1_noise = dH1_noise / Z;



