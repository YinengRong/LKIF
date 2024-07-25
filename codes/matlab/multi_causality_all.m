function IF_result=multi_causality_all(X,np,dt,max_lag,series_temporal_order,significance_test)
% On input:
%    X: matrix storing the M time series (each as Nx1 column vectors)
%    max_lag: time order of lags (default 1)
%    >=1 lag time step
%    -1 Determine the lag order of the data based on the AIC criterion.
%   - 2 Determine the lag order of the data based on the BIC criterion.
%    np: integer >=1, time advance in performing Euler forward 
%	 differencing, e.g., 1, 2. Unless the series are generated
%	 with a highly chaotic deterministic system, np=1 should be
%	 used. 
%    (default 1)
%    dt: frequency of sampling (default 1)
%
%    series_teporal_order: Nx1 column vectors, records the timestamp of 
%    each sample, with a minimum sampling interval of dt (used for 
%    panel data, or missing measurement data).
%    (default [])
%    significance_test:  1  do the significance test (default)
%           =0  not (to save the computation time)
% On output:
%    a structure value IF_result with sub
%    IF:  information flow
%    nIF: normalized information flow
%    SEIF: standard error of information flow
%    errr.e90/e95/e99:standard error at 90/95/99% confidence level
%    dnoise: dH1_noise/dt with correlated noise
%    dnoise_old: dH1_noise/dt without correlated noise (##  been commented out.)
%    nIF: normalized information flow without correlated noise
%    p: p-value of information flow
%
% Citations: 
%    X.S. Liang, 2016: Information flow and causality as rigorous notions
%    ab initio. Phys. Rev. E, 94, 052201.
%    X.S. Liang, 2014: Unraveling the cause-effect relation between time series. Phys. Rev. E 90, 052150.
%    X.S. Liang, 2015: Normalizing the causality between time series. Phys. Rev. E 92, 022126.
%    X.S. Liang, 2021: Normalized Multivariate Time Series Causality Analysis and Causal Graph Reconstruction. Entropy. 23. 679.

% Note: This is an alternative and simplified version of
%   multi_causality_est.m, which was originally coded by X.S. Liang
%   Here all the causal relations are inferred once for all.
% 
% Author: Yineng Rong (yinengrong@foxmail.com)

switch nargin
    case 1
        np=1;dt=1;max_lag=1;series_temporal_order=[];significance_test=1;
    case 2
        dt=1;np=1;series_temporal_order=[];significance_test=1;
    case 3
        dt=1;series_temporal_order=[];significance_test=1;
    case 4
        series_temporal_order=[];significance_test=1;
    case 5
        significance_test=1;
end



%% data pre-processing
[m,n]=size(X);
if isempty(series_temporal_order)
    series_temporal_order=1:dt:dt*n;
end    
    
%max lag
if max_lag<0
    max_lag=lag_order(X,series_temporal_order/dt,max_lag);
end
% panel data
q1=max_lag+1;
XX = zeros(m,q1,n+q1-1);
for k = 0:q1-1
    XX(:,k+1,k+1:k+n) = X; % k-lagged observations
end
XX=XX(:,:,1:n);
[XL,XR]=pre_panel_data(XX,series_temporal_order/dt,max_lag);

[m2,n2]=size(XL);
X1=[XR;ones(1,n2)];% [X and 1], where the 1 is a const
X2=[XL;XR(1:m2*(max_lag-1),:);ones(1,n2)];% [X and 1], where the 1 is a const


%% calculate information flow
%information flow(IF)
A0=X2/X1;
A=(A0-eye(m2*max_lag+1))/dt;

At=A(1:end-1,1:end-1);
C=cov(X1(1:end-1,:)');
IF=C./diag(C).*At;           

%standarized error of IF
E= X2(1:end,:)-A0(1:end,1:end)*X1(1:end,:);
SIG = (E*E')/(n-np-m2-1);                                                 % in Liang (2014,2021) it is SIG = diag(diag((E*E')/(n0-np)));
se_a=sqrt(diag(SIG(1:end-1,1:end-1))*diag(inv(X1(1:end-1,:)*X1(1:end-1,:)'))');
SE_IF=se_a.*abs(C./diag(C))/dt;
p=(1-normcdf(abs(IF./(SE_IF))))*2;

%significance_test
if significance_test==1
    z99 = 2.56;z95 = 1.96;z90 = 1.65;
    e90_IF=SE_IF*z90;
    e95_IF=SE_IF*z95;
    e99_IF=SE_IF*z99;
end

% normalized IF
B=sqrt(abs(SIG(1:end-1,1:end-1)/dt));
%dH1noise/dt
dH1_noise = sum(B.^2) ./ (2. * diag(C)');  
dH1_noise_old = (diag(B.^2) ./ (2. * diag(C)))';                           %in Liang(2015, 2021)
nIF=IF./(sum(abs(IF)')+abs(dH1_noise))';
nIF_old=IF./(sum(abs(IF)')+abs(dH1_noise_old))';


%% output 

% Removing causality from future to past when max_lag > 1
IF_result.IF      =     temporal_lag(IF,m2,max_lag);
IF_result.nIF     =    temporal_lag(nIF,m2,max_lag);
IF_result.SEIF    =  temporal_lag(SE_IF,m2,max_lag);
IF_result.err.e90 = temporal_lag(e90_IF,m2,max_lag);
IF_result.err.e95 = temporal_lag(e95_IF,m2,max_lag);
IF_result.err.e99 = temporal_lag(e99_IF,m2,max_lag);
IF_result.p       =      temporal_lag(p,m2,max_lag);
IF_result.dnoise  =  dH1_noise(end-m2+1:end);
IF_result.dnoise_old = dH1_noise_old(end-m2+1:end)';
IF_result.nIF_old     =    temporal_lag(nIF_old,m2,max_lag);
IF_result.B       =                               B;
end

function X=temporal_lag(X,m1,max_lag)
    %X=flip(reshape(X(end-m1+1:end,:),[m1,m1,max_lag]),3);
    X=reshape(X(1:m1,:),[m1,m1,max_lag]);
end

function max_lag=lag_order(X,t,option)

[n,m] = size(X);
X = X-mean(X,2); % no constant term, normalise

% store lags
morder=1:min([max([floor(m/(n^2+n)),1]),20]);
q = max(morder);
q1 = q+1;
XX = zeros(n,q1,m+q);
for k = 0:q
    XX(:,k+1,k+1:k+m) = X; % k-lagged observations
end

nummo = length(morder);

aic = nan(nummo,1);
bic = nan(nummo,1);

% loop through model orders

for i = 1:nummo

    q = morder(i);



    if q >= m
        fprintf(2,'  WARNING: model order too large (must be < %d)\n',m);
        continue
    end

    %panel data
    XX=XX(:,:,1:m);
    [X0,XL]=pre_panel_data(XX,t,q);
    
    M = size(X0,2);
    
    % stack lags

    owstate = warn_supp;
    A = X0/XL;                                         % OLS using QR decomposition
    [iw,~,wid]= warn_test(owstate,[],false);
    if isbad(A)                                        % something went badly wrong
        fprintf(2,'  WARNING: VAR estimation failed\n');
        continue % show-stopper
    end
    if iw % rank-deficient?
        fprintf(2,'  WARNING: VAR estimation may be problematic (%s)',wid);
        % not necessarily a show-stopper - carry on
    end

    E   = X0-A*XL;            % residuals
    DSIG = det((E*E')/(M-1)); % residuals covariance matrix determinant
    if DSIG <= 0
        if ~iw, if ~verb, fprintf('model order = %d',q); end; end
        fprintf(2,'  WARNING: residuals covariance not positive definite\n');
        continue % show-stopper
    end

    [aic(i),bic(i)] = infocrit(-(M/2)*log(DSIG),q*n*n,M); % -(M/2)*log(DSIG) is max log-likelihood
end
if option==-1
    max_lag=find(aic==min(aic));
else
    max_lag=find(bic==min(bic));
end

end

function [aic,bic] = infocrit(L,k,m)

if m-k-1 <= 0
    aic = NaN;
else
    aic = -2*L + 2*k*(m/(m-k-1)); % Note AIC without correction = -2*L + 2*k
end
bic = -2*L + k*log(m);


end


function [X0,XL]=pre_panel_data(XX,t,q)
    [n,~,m]=size(XX);
    if ~isempty(t)
        assert(length(t)==m, 'length of series order must be same as the data X');
        tt = zeros(q+2,m+q)-pi*pi;
        for k = 0:q
            tt(k+1,k+1:k+m) = t; % k-lagged observations
        end
        tt(q+2,:)=tt(q+1,:)-1;
        ind=find(max(abs(diff(tt)+1))==0);
    else 
        ind=q+1:m;
    end
    
    
    q1 = q+1;
    M = length(ind);
    nq = n*q;

    % stack lags

    
    
    X0 = reshape(XX(:,1,ind,:),n,M);
    XL = reshape(XX(:,2:q1,ind,:),nq,M);
end
