function [T21, err90, err95, err99] = multi_causality_est(xx, np)
% function [T21, err90, err95, err99] = causality_est(x, np)
% Estimate T21, the information transfer from series X2 to series X1, among
% M time series (stored as column vectors in X).
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
%    T21:  info flow from X2 to X1	(Note: Not X1 -> X2!)
%    err90: standard error at 90% significance level
%    err95: standard error at 95% significance level
%    err99: standard error at 99% significance level
%
% Citations: 
%    X.S. Liang, 2016: Information flow and causality as rigorous notions
%    ab initio. Phys. Rev. E, 94, 052201.
%    X.S. Liang, 2014: Unraveling the cause-effect relation between time series. Phys. Rev. E 90, 052150.
%    X.S. Liang, 2015: Normalizing the causality between time series. Phys. Rev. E 92, 022126.
%    X.S. Liang, 2021: Normalized Multivariate Time Series Causality Analysis and Causal 
%                        Graph Reconstruction. Entropy. 23. 679.
dt = 1;	

[nm, M] = size(xx);

dx1(:,1) = (xx(1+np:nm, 1) - xx(1:nm-np, 1)) / (np*dt);
x(:,1:M) = xx(1:nm-np, 1:M);


clear xx;

N = nm-np;





C = cov(x);		% MxM matrix



for k = 1 : M
%dC(1,1) = sum((x(:,1) - mean(x(:,1))) .* (dx1 - mean(dx1))); 
dC(k,1) = sum((x(:,k) - mean(x(:,k))) .* (dx1 - mean(dx1))); 
end
dC = dC / (N-1);


a1n = inv(C) * dC;

a12 = a1n(2,1);


%
% Information transfer: T21 = C12/C11 * a12
%

 T21 = C(1,2)/C(1,1) * a12;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

f1 = mean(dx1);
for k = 1 : M
    f1 = f1 - a1n(k,1) * mean(x(:,k));
end

R1 = dx1 - f1;
for k = 1 : M
    R1 = R1 - a1n(k,1) * x(:,k);
end

Q1 = sum(R1 .* R1);

b1 = sqrt(Q1 * dt / N);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% covariance matrix of the estimator of (f1, a11, a12,..., a1M, b1)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



NI(1,1) = N * dt /b1/b1;
NI(M+2,M+2) = 3*dt/b1^4 *sum(R1 .* R1) - N/b1/b1;

for k = 1 : M
    NI(1,k+1) = dt/b1/b1 * sum(x(:,k));
end
NI(1,M+2) = 2*dt/b1^3 *sum(R1);

for k = 1 : M
  for j = 1 : M
     NI(j+1,k+1) = dt/b1/b1 * sum(x(:,j) .* x(:,k));
  end
end

for k = 1 : M
    NI(k+1,M+2) = 2*dt/b1^3 * sum(R1 .* x(:,k));
end

for j = 1 : M+2
    for k = 1 : j
    NI(j,k) = NI(k,j);
    end
end

invNI = inv(NI);
var_a12 = invNI(3,3);

%
% approx. variance of a12, corr. to variance of T21



%

var_T21 = (C(1,2) / C(1,1))^2 * var_a12;

% disp([ abs(T21)  abs(a1n(1,1)) abs(b1^2 / (2. * C(1,1))) abs(T21) + abs(a1n(1,1)) + abs(b1^2 / (2. * C(1,1)))+sum(abs(C(1,3:end)/C(1,1) .* a1n(3:end)'))  ])
% disp([(T21) + (a1n(1,1)) + (b1^2 / (2. * C(1,1)))+sum((C(1,3:end)/C(1,1) .* a1n(3:end)'))  ])
%
% From the standard normal distribution table, 
% at level alpha=95%, z=1.96
%		 99%, z=2.56
%		 90%, zx=1.65
%

	z99 = 2.56;
	z95 = 1.96;
	z90 = 1.65;

	err90 = sqrt(var_T21) * z90;
	err95 = sqrt(var_T21) * z95;
	err99 = sqrt(var_T21) * z99;


