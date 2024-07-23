function [TAB, TBA] = IF_subspace_v1(xx, r, s, np)
% 
% function [TAB, TBA] = IF_subspace(x, r, s, np)
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
% On output:
%    TAB:  info flow from subspace A to subspace B
%    TBA:  info flow from subspace B to subspace A
%    err90: standard error at 90% confidence level
%    err95: standard error at 95% confidence level
%    err99: standard error at 99% confidence level
%
% Citations: 
% Liang, X.S. The causal interaction between complex subsystems. 
%             Entropy, 2022, 24, 3.
%
% Note: Here the symbols M and NL correspond respectively 
%         to the symbols N and K in the paper.


dt = 1;		% dt is set to 1 here; 
		% the final result should be divided by dt.


[nm, M] = size(xx);

for i=1:M,  dx1(:,i) = (xx(1+np:nm, i) - xx(1:nm-np, i)) / (np*dt);  end
x(:,1:M) = xx(1:nm-np, 1:M);


clear xx;

NL = nm-np;



C = cov(x);		% MxM covariance matrix

Cr=C(1:s,1:s);  Cr(1:s,1:r)=0;  Cr(1:r,1:s)=0;  for k=1:r, Cr(k,k) = 1; end
	% C(1:s,1:s) with first r rows and first r columns removed,
	% replaced by an identity submatrix

Crr=C(1:s,1:s);  Crr(1:s,r+1:s)=0;  Crr(r+1:s,1:s)=0;  for k=r+1:s, Crr(k,k) = 1; end
	% C(1:s,1:s) with last s-r rows and last s-r columns removed,
	% replaced by an identity submatrix.


for i = 1 : M,
for k = 1 : M,
dC(k,i) = sum((x(:,k) - mean(x(:,k))) .* (dx1(:,i) - mean(dx1(:,i) ))); 
end
end
dC = dC / (NL-1);

invC = inv(C);
invCr = inv(Cr);


for i = 1 : M,
A(:,i) = invC * dC(:,i);
end

A = A';		% a(i,k)	ith row in the dyn sys eqn


%
% Rate of information flow from subspace A to subspace B; 
% see Eq. (28) in Liang (2022).
%
% TAB = \sum_{i=r+1}^s  [\sum_{j=r+1}^s  invCr(i,j) * 
%	(\sum_{k=1}^s A(i,k)*C(k,j))  -   A(i,i)]
%



   AC = A(1:s,1:s) * C(1:s,1:s); 		%  A(i,:) * C(:,j)


TAB = trace(invCr(r+1:s,r+1:s) * AC(r+1:s,r+1:s)') - trace(A(r+1:s,r+1:s));
 


invCrr = inv(Crr);
TBA = trace(invCrr(1:r,1:r) * AC(1:r,1:r)') - trace(A(1:r,1:r));




%%%%%%%%% Significance test %%%%%%%%%%
%%%%%%%%% To be implemented %%%%%%%%%%
	err90 = 0;
	err95 = 0;
	err99 = 0;





