function [t21,T21,dH1_noise,err90,err95,err99] = multi_causality_est_all_panel_tau(xx, np,st,t)
% recoded by Rong Yineng (yinengrong@foxmail.com)
% function [t21,T21, err90,err95,err99] = multi_causality_est_all(xx, np,st)
% Estimate T21, the information transfer from series Xk to series X1, among
% M(1<k<=M) time series (stored as column vectors in X).
% dt is taken to be 1.
%
% On input:
%    XX: matrix storing the M time series (each as Nx1 column vectors)
%       X1 and X2 stored as the first two column vectors
%    np: integer >=1, time advance in performing Euler forward 
%	 differencing, e.g., 1, 2. Unless the series are generated
%	 with a highly chaotic deterministic system, np=1 should be
%	 used. 
%    st:  st=1  do the significance test
%           =0  not (to save the computation time)
% On output:
%    t21:  normalized information flow
%    T21:  information flow
%    err90: standard error at 90% significance level
%    err95: standard error at 95% significance level
%    err99: standard error at 99% significance level
%
% Citations: 
%    X.S. Liang, 2016: Information flow and causality as rigorous notions
%    ab initio. Phys. Rev. E, 94, 052201.
%    X.S. Liang, 2014: Unraveling the cause-effect relation between time series. Phys. Rev. E 90, 052150.
%    X.S. Liang, 2015: Normalizing the causality between time series. Phys. Rev. E 92, 022126.
%    X.S. Liang, 2021: Normalized Multivariate Time Series Causality Analysis and Causal Graph Reconstruction. Entropy. 23. 679.

dt = 1;	
[nm,~] = size(xx);
j=1;
for i=2:nm
    if t(i)-t(i-1)==dt
        ind(j)=i;j=j+1;
    end
end

[nm, M] = size(xx);

for k = 1 : M
    dx(:,k) = (xx(ind, k) - xx(ind-np, k)) / (np*dt);
end
x(:,1:M) = xx(ind-np, 1:M);

clear xx;
[nm, M] = size(x);
N = nm;



C = cov(x);		% MxM matrix


for kk=1 : M
    for k = 1 : M
    %dC(1,1) = sum((x(:,1) - mean(x(:,1))) .* (dx1 - mean(dx1))); 
    dC(k,kk) = sum((x(:,k) - mean(x(:,k))) .* (dx(:,kk) - mean(dx(:,kk)))); 
    end
end
dC = dC / (N-1);


ann = inv(C) * dC;

% a12 = a1n(2,1);


%
% Information transfer: T21 = C12/C11 * a12
%

% T21 = C(1,2)/C(1,1) * a12;
T21 =zeros(M,M);
 for k=1:M
    T21(k,:) = C(k,:)./C(k,k) .* ann(:,k)';
%     T21(k,k)=NaN;
 end

 


f = mean(dx);
f = f - sum((ann'.*repmat(mean(x),[M,1])),2)';


R = dx - repmat(f,[N,1]);
R = R - x * ann;

Q = sum(R .* R);

b = sqrt(Q * dt / N);
 if st==1
 
%%%%%%%%%%%%significance test%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% covariance matrix of the estimator of (f1, a11, a12,..., a1M, b1)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



NI(1,1,:) = N * dt ./b.^2;
NI(M+2,M+2,:) = 3*dt./b.^4 .*Q - N./b.^2;

for k = 1 : M
    NI(1,k+1,:) = dt./b.^2 * sum(x(:,k));
end
NI(1,M+2,:) = 2*dt./b.^3 .*sum(R);

for k = 1 : M
  for j = 1 : M
     NI(j+1,k+1,:) = dt./b.^2 .* sum(x(:,j) .* x(:,k));
  end
end

for k = 1 : M
    NI(k+1,M+2,:) = 2*dt./b.^3 .* sum(R .* repmat(x(:,k),[1,M]));
end

for j = 1 : M+2
    for k = 1 : j
    NI(j,k,:) = NI(k,j,:);
    end
end

%
% approx. variance of a12, corr. to variance of T21



for k=1:M
tmp=diag(inv(NI(:,:,k)));
var_a12(k,:)=tmp(2:M+1);
end


% var_T21 = (C(1,2) / C(1,1))^2 * var_a12;

%var_T21 = (C(1,:) / C(1,1)).^2 .* var_a12;
for  k=1:M
    var_T21(k,:) = (C(k,:)./C(k,k)).^2 .* var_a12(k,:);
end
% p=1-(1-normcdf(abs(T21./sqrt(var_T21))))*2;

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
 else
    err90 = [];
	err95 = [];
	err99 = [];
 end
dH1_noise = b.^2 ./ (2. * diag(C)');
t21=T21./(sum(abs(T21)')+abs(dH1_noise))';

