'
multi_causality_est <- function(X, max_lag=1, np=1, dt=1, series_temporal_order=NULL, significance_test=1) 
source("LK_Info_Flow.R")
IF_result=multi_causality_est(X, np=1, dt=1, max_lag=1, series_temporal_order=NULL, significance_test=1):

Estimate information flow, the information transfer from columns to raws among
M time series (stored as a MXM matrix).

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

IF_result$XXX:

   max_lag:          time order of lags (in IF)

bellowing outputs are Estimate information flow, the information transfer from columns to raws among
M time series (stored as a M X M X max_lag matrix).

   IF:               information flow 
   nIF:              normalized information flow
   SEIF:             standard error of information flow
   err_e90/e95/e99:  standard error at 90/95/99# significance level
   dnoise:           dH1_noise/dt with correlated noise
   p:                p-value of information flow
'  

library(stats)
library(MASS)
library(igraph)


multi_causality_est <- function(X, max_lag=1, np=1, dt=1, series_temporal_order=NULL, significance_test=1) {
    res <- list()
    dm <- dim(X)
#    m <- dm[1]
#    n <- dm[2]
  m = nrow(X)
  n = ncol(X)


  if (is.null(series_temporal_order)) {
    series_temporal_order <- seq(from = 1, to = (dt * n), by = dt)
  }

  if (max_lag < 0) {
    max_lag = lag_order(X, series_temporal_order/dt, max_lag)
  }

  q1 = max_lag + 1
  XX = array(0, dim = c(m, q1, n+q1-1))

  for (k in 1:q1) {
    XX[, k, k:(k+n-1)] <- X
  }

  XX <- XX[, , 1:n]


  panel_data = pre_panel_data(XX, series_temporal_order/dt, max_lag)
  XL = panel_data$X0
  XR = panel_data$XL

  m2 = nrow(XL)
  n2 = ncol(XL)


  X1 = rbind(XR, matrix(1, nrow = 1, ncol = n2))
  if (max_lag>1){
    X2 = rbind(XL, XR[1:(m2 * (max_lag-1)), ], matrix(1, nrow = 1, ncol = n2))
    }
    else {
      X2 = rbind(XL, matrix(1, nrow = 1, ncol = n2))
    }




  A0 = X2 %*% ginv(X1)
  A = (A0 - diag(m2*max_lag+1)) / dt
  At = A[1:(nrow(A)-1), 1:(ncol(A)-1)]
  C = cov(t(X1[1:(nrow(X1)-1), ]))#cov(X1[1:(nrow(X1)-1), ])


  IF = ((C / diag(C))) * At

  E = X2 - A0 %*% X1
  SIG = E %*% t(E) / (n - np - m2 - 1)

  B = sqrt(abs(SIG[1:(nrow(SIG)-1), 1:(ncol(SIG)-1)]) / dt)
  dH1_noise = colSums(B^2) / (2 * diag(C))
  Z = rowSums(abs(IF)) + abs(dH1_noise)
  Z = matrix(rep(Z, m2*max_lag), nrow = m2*max_lag, ncol = m2*max_lag)
  nIF = IF / Z




  if (significance_test == 1) {
    diag_SIG <- diag(SIG[-nrow(SIG), -ncol(SIG)])
    reshaped_sig <- matrix(diag_SIG, ncol = 1, nrow = m2*max_lag)
    diag_x <- diag(solve(X1[-nrow(X1), ] %*% t(X1[-nrow(X1), ])))
    reshaped_x <- matrix(diag_x,ncol = m2*max_lag, nrow = 1)
    se_a <- sqrt(reshaped_sig %*% reshaped_x)


    SE_IF <- se_a * abs(C / diag(C)) / dt  # standarized error of IF
    p <- (1 - pnorm(abs(IF / SE_IF))) * 2  # p-value

    z90 = 1.65
    z95 = 1.96
    z99 = 2.56

    e90_IF = SE_IF * z90
    e95_IF = SE_IF * z95
    e99_IF = SE_IF * z99
    res=list('IF' = temporal_lag(t(IF), m2,max_lag), 
         'nIF' = temporal_lag(t(nIF), m2,max_lag), 
         'dnoise' = dH1_noise[1:(length(dH1_noise)-m2)], 
         'max_lag' = max_lag, 
         'SEIF' = temporal_lag(t(SE_IF), m2,max_lag), 
         'err_e90' = temporal_lag(t(e90_IF), m2,max_lag), 
         'err_e95' = temporal_lag(t(e95_IF), m2,max_lag), 
         'err_e99' = temporal_lag(t(e99_IF), m2,max_lag),  
         'p' = temporal_lag(t(p), m2,max_lag)
    )
  } else {
    res=list('IF' = temporal_lag(IF, m2,max_lag), 
         'nIF' = temporal_lag(nIF, m2,max_lag), 
         'dnoise' = temporal_lag(dH1_noise, m2,max_lag), 
         'max_lag' = max_lag 
    )
  }

  return(res)
}


groups_est <- function(xx, ind, dt=1, np=1){
  if (length(dim(as.matrix(xx))) != 2){
     stop("multi_est need inputs be like: if with x1,x2:  [ x1[:,] x2[:,m] ], [ x1[:,m] x2[:,] ], [ x1[:,m] x2[:,n] ], [ x1[:,] x2[:,] ], if with x only: x[:,0]->x1 x[:,1]->x2 x[:,1:]->others")
  }
  N <- dim(xx)
  nm <- N[1]; M <- N[2];

  dx <- matrix(0, nm-np, M)
  for (k in 1:M){
    dx[,k] <- (xx[(np+1):length(xx[,k]),k] - xx[1:(nm-np),k])/np*dt
  }

  x <- xx[1:(nm-np),]
  rm(xx)

  NL <- length(x[,1])

  C <- cov(x)

  dC <- matrix(0, M, M)

  for (kk in 1:M){
    for (k in 1:M){
      dC[k,kk] <- sum((x[,k] - mean(x[,k])) * (dx[,kk] - mean(dx[,kk])))
    }
  }

  dC <- dC/(NL-1)

  # Try-catch equivalent in R for linalg.inv()
  ann <- tryCatch(solve(C) %*% dC, 
                  error = function(e) matrix(Inf, M, M))


  Cr <- matrix(0, M, M)
  Crr <- Cr

  Crr[1:ind[1], 1:ind[1]] <- C[1:ind[1], 1:ind[1]]

  for (i in (ind[1]+1):ind[2]){
    Crr[i,i] <- 1
  }

  Cr[(ind[1]+1):ind[2], (ind[1]+1):ind[2]] <- C[(ind[1]+1):ind[2], (ind[1]+1):ind[2]]

  for (i in 1:ind[1]){
    Cr[i,i] <- 1
  }

  invCr <- solve(Cr)
  invCrr <- solve(Crr)
  AC <- ann[1:ind[2],1:ind[2]] %*% C[1:ind[2],1:ind[2]]
  TBA <- sum(diag(invCr[(ind[1]+1):ind[2], (ind[1]+1):ind[2]] %*% AC[(ind[1]+1):ind[2], (ind[1]+1):ind[2]])) - sum(diag(ann[(ind[1]+1):ind[2], (ind[1]+1):ind[2]]))

  TAB <- sum(diag(invCrr[1:ind[1],1:ind[1]] %*% AC[1:ind[1],1:ind[1]])) - sum(diag(ann[1:ind[1], 1:ind[1]]))


  IF_result <- list('TAB' = TAB, 'TBA' = TBA)
  return(IF_result)
}


plot_causal_graph <- function(A, max_lag, f_name = "model",plot_style='fr') {  
    # A=>matrix
    A <- as.matrix(A)  
      
    N <- dim(A)  
    edge_w<-v_t<-v_f<-NULL
    for (i in 1:N[1]){
      for (j in 1:N[2]) {
        if (is.na(A[i,j])){}
        else {
           v_f=c(v_f,i)
           v_t=c(v_t,j)
           edge_w=c(edge_w,round(A[i,j],2))
        }
      }
    }
    v_all<-unique(c(v_f,v_t))



    v_name<-NULL
    for (i in 1:length(v_all)){
      if (v_all[i]<=N[1]){
        v_name=c(v_name,paste('',as.character(v_all[i]-1),'',seq=''))
      }
      else {
        lag_num=floor((v_all[i]-1)/N[1])
        v_name=c(v_name,paste('',as.character(v_all[i]-lag_num*N[1]-1),'(-',as.character(max_lag-lag_num),')',seq=''))
      }
    }
    
    edges <- data.frame(from = v_f, to = v_t)
    g <- graph_from_data_frame(edges, directed = TRUE)
    V(g)$name <- v_name
    E(g)$weight <-edge_w
    png(filename = paste(f_name,".png",sep=''), width = 1200, height = 800, units = "px", res = 100)
    if (plot_style=='fr'){
      plot(g, edge.label = E(g)$weight, vertex.label = V(g)$name, vertex.label.color = "black", vertex.size = 10, edge.arrow.size = 0.5,layout= layout_with_fr(g))  # layout_with_fr，layout_with_kk
    }
    else {
      plot(g, edge.label = E(g)$weight, vertex.label = V(g)$name, vertex.label.color = "black", vertex.size = 10, edge.arrow.size = 0.5,layout= layout_with_kk(g))  # layout_with_fr，layout_with_kk
    }
    dev.off()
}

causal_graph <- function(causal_matrix, significance = NULL, c_threshold = 0.001, s_threshold = 0.01, f_name = 'model',plot_style='fr') {  
    a <- abs(causal_matrix)  
    if (!is.null(significance)) {  
        b <- abs(significance)  
        a[b > s_threshold] <- NA  
    }  
      
    a[a < c_threshold] <- NA
    N=dim(a)
    cc <- matrix(NA, nrow = N[1], ncol = N[2] * N[3])
    for (i in 1:N[2]) {
      j<-(i-1)*N[3]+1
      cc[,j:(j+N[3]-1)]<-a[,i,]
    }
    for (i in 1:N[1]) {
      cc[i,i]<- NA
    }
    plot_causal_graph(cc, N[3], f_name = f_name,plot_style=plot_style )
}

pre_panel_data <- function(XX, t, q){
  dims = dim(XX)
  n = dims[1]
  m = dims[3]

  if(length(t) == 0){
    t = 1:m
  }

  tt = matrix(-pi^2, nrow = q+2, ncol = m+q)

  for(k in 1:(q+1)){
    tt[k, k:(k+m-1)] = t
  }

  tt[q+2, ] = tt[q+1, ] - 1
  dt = abs(diff(tt) + 1)
  max_dt = apply(dt, 2, max)
  dims <- c(1, length(max_dt))
  max_dt <- array(max_dt,dim=dims)
  ind = which(max_dt <0.0000001)
    dims <- c(1, length(ind))
  ind <- array(ind,dim=dims)
  
  
  M = length(ind)

  nq = n * q

  X0 = matrix(nrow = n, ncol = M)

  for(i in 1:n){
    X0[i, ] = XX[i, 1, ind]
  }
  XL = array(dim = c(nq, M))

  sub_XX <- XX[, 2:(q+1), ind]
  #reordered_XX <- aperm(sub_XX, c(2, 1, 3)) 
  XL <- array(data = as.vector(sub_XX), dim = c(nq, M))
  
  return(list(X0=X0, XL=XL))
}



infocrit <- function(L, k, m){
    if(m - k - 1 <= 0){
        aic <- NaN
    } else {
        aic <- -2 * L + 2 * k * (m/(m - k - 1))
    }
    bic <- -2 * L + k * log(m)
    return(list(aic=aic, bic=bic))
}



lag_order <- function(X, t, option) {

  n <- nrow(X)
  m <- ncol(X)

  X <- sweep(X, 1, rowMeans(X))

  morder <- 1:min(max(floor(m/(n^2 + n)), 1), 20)
  nummo <- length(morder)

  aic <- rep(NaN, nummo)
  bic <- rep(NaN, nummo)

  q <- max(morder)
  q1 <- q + 1

  XX <- array(0, c(n, q1, m + q))
  for (k in 0:q) {
    XX[, k + 1, (k + 1):(k + m)] <- X
  }

  for (i in 1:nummo) {
    q <- morder[i]

    if (q >= m) {
      print(paste('WARNING: model order too large (must be <', m, ')'))
      next
    }

    XX <- XX[, , 1:m]

    pre_panel_data_output <- pre_panel_data(XX, t, q)
    X0 <- pre_panel_data_output$X0
    XL <- pre_panel_data_output$XL

    M <- ncol(X0)

    A <- ginv(t(XL)) %*% t(X0)
    E <- X0 - (XL %*% A)

    DSIG <- det(t(E) %*% E / (M - 1))

    if (DSIG <= 0) {
      print('WARNING: residuals covariance not positive definite')
      next
    }

    infocrit_output <- infocrit(-(M/2) * log(DSIG), q*n*n, M)
    aic[i] <- infocrit_output$aic
    bic[i] <- infocrit_output$bic
  }

  if (option == -1) {
    max_lag <- min(which(aic == min(aic, na.rm = TRUE)))
  } else {
    max_lag <- min(which(bic == min(bic, na.rm = TRUE)))
  }

  return(max_lag)
}




temporal_lag <- function(X, m1, max_lag) {
  # Removing causality from future to past when max_lag > 1
  #return(array(X[1:m1,], dim = c(m1, m1, max_lag)))
  XX <- aperm((array(X[1:m1,], dim = c(m1, max_lag, m1))),c(1,3,2))
  return (XX)
}
