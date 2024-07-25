function [T21, err90, err95, err99] = causality_est(xx1, xx2, np)
% 
% function [T21, err90, err95, err99] = causality_est(x1, x2, np)
%
% Estimate T21, the information transfer from series X2 to series X1 
% dt is taken to be 1.
%
% On input:
%    X1, X2: the series (n by 1 colum vectors)
%    np: integer >=1, time advance in performing Euler forward 
%	 differencing, e.g., 1, 2. Unless the series are generated
%	 with a highly chaotic deterministic system, np=1 should be
%	 used. 
%
% On output:
%    T21:  info flow from X2 to X1	(Note: Not X1 -> X2!)
%    err90: standard error at 90% confidence level
%    err95: standard error at 95% confidence level
%    err99: standard error at 99% confidence level
%
% Citations: 
%    X. San Liang, 2014: Unraveling the cause-effect relation between time series. Phys. Rev. E 90, 052150.
%    X. San Liang, 2015: Normalizing the causality between time series. Phys. Rev. E 92, 022126.

% Note: This is the original matlab program for causality estimation 
% for 2 time series by X.S. Liang in 2014. 
% Authorized by X.S. Liang to release here.


dt = 1;	


[nm, one] = size(xx1);


dx1(:,1) = (xx1(1+np:nm, 1) - xx1(1:nm-np, 1)) / (np*dt);
 x1(:,1) = xx1(1:nm-np, 1);

dx2(:,1) = (xx2(1+np:nm, 1) - xx2(1:nm-np, 1)) / (np*dt);
 x2(:,1) = xx2(1:nm-np, 1);

clear xx1 xx2;

N = nm-np;



C = cov(x1, x2);

dC(1,1) = sum((x1 - mean(x1)) .* (dx1 - mean(dx1))); 
dC(1,2) = sum((x1 - mean(x1)) .* (dx2 - mean(dx2))); 
dC(2,1) = sum((x2 - mean(x2)) .* (dx1 - mean(dx1))); 
dC(2,2) = sum((x2 - mean(x2)) .* (dx2 - mean(dx2))); 
dC = dC / (N-1);


% C_infty = cov(x1, x2);
C_infty = C;

detc = det(C);

a11 = C(2,2) * dC(1,1) - C(1,2) * dC(2,1);
a12 = -C(1,2) * dC(1,1) + C(1,1) * dC(2,1);
% a21 = -C(1,2) * dC(1,2) + C(1,1) * dC(2,2);
% a22 = C(2,2) * dC(1,2) - C(1,2) * dC(2,2);

a11 = a11 / detc;
a12 = a12 / detc;
% a21 = a21 / detc;
% a22 = a22 / detc;

f1 = mean(dx1) - a11 * mean(x1) - a12 * mean(x2);
% f2 = mean(dx2) - a21 * mean(x1) - a22 * mean(x2);

R1 = dx1 - (f1 + a11*x1 + a12*x2);
% R2 = dx2 - (f2 + a21*x1 + a22*x2);

Q1 = sum(R1 .* R1);
% Q2 = sum(R2 .* R2);

b1 = sqrt(Q1 * dt / N);
% b2 = sqrt(Q2 * dt / N);

%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% covariance matrix of the estimation of (f1, a11, a12, b1)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
NI(1,1) = N * dt / b1/b1;
NI(2,2) = dt/b1/b1 * sum(x1 .* x1);
NI(3,3) = dt/b1/b1 * sum(x2 .* x2);
NI(4,4) = 3*dt/b1^4 * sum(R1 .* R1) - N/b1/b1;
NI(1,2) = dt/b1/b1 * sum(x1);
NI(1,3) = dt/b1/b1 * sum(x2);
NI(1,4) = 2*dt/b1^3 * sum(R1);
NI(2,3) = dt/b1/b1 * sum(x1 .* x2);
NI(2,4) = 2*dt/b1^3 * sum(R1 .* x1);
NI(3,4) = 2*dt/b1^3 * sum(R1 .* x2);

NI(2,1) = NI(1,2);
NI(3,1) = NI(1,3);    NI(3,2) = NI(2,3);
NI(4,1) = NI(1,4);    NI(4,2) = NI(2,4);   NI(4,3) = NI(3,4);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

invNI = inv(NI);
var_a12 = invNI(3,3);		
%
% approx. variance of a12, corr. to variance of T21
%


%
% Information transfer: T21 = C12/C11 * a12
%			T12 = C21/C22 * a21
%
 T21 = C_infty(1,2)/C_infty(1,1) * (-C(2,1)*dC(1,1) + C(1,1)*dC(2,1)) / detc;
% T12 = C_infty(2,1)/C_infty(2,2) * (-C(1,2)*dC(2,2) + C(2,2)*dC(1,2)) / detc;
%

var_T21 = (C_infty(1,2)/C_infty(1,1))^2 * var_a12;


%
% From the standard normal distribution table, 
% significance level alpha=95%, z=1.96
%		           99%, z=2.56
%			   90%, z=1.65
%
	z99 = 2.56;
	z95 = 1.96;
	z90 = 1.65;
    z80 = 1.28;
    
	err90 = sqrt(var_T21) * z90;
	err95 = sqrt(var_T21) * z95;
	err99 = sqrt(var_T21) * z99;

