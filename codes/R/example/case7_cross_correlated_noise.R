# case 7 (data with cross-correlated noise):
source("../LK_Info_Flow.R")
library(R.matlab)
# To generate the cross-correlated noise, we use functions to generate noise based on the MVGC MATLAB toolbox by Barnett and Seth.

# Define the function for generating cross-correlated noise
var_to_tsdata <- function(A, SIG, m, N=1, mtrunc=NULL, decayfac=100) {
    if (length(dim(A)) == 1) {
        n <- length(A)
        A <- matrix(A, nrow=1)
    } else {
        n <- ncol(A)
    }

    if (is.null(mtrunc)) {
        if (is.null(decayfac)) {
            decayfac <- 100
        }
        rho <- var_specrad(A)
        if (rho >= 1) {
            stop("unstable VAR")
        }
        mtrunc <- as.integer((log(.Machine$double.eps) - decayfac) / log(rho))
    } else {
        stopifnot(is.numeric(mtrunc), is.integer(mtrunc), mtrunc >= 0)
    }

    tryCatch({
        C <- chol(SIG, lower=TRUE)
    }, error = function(e) {
        stop("covariance matrix not positive-definite")
    })

    if (N > 1) {  # multi-trial
        X <- array(0, dim=c(n, m, N))
        E <- array(0, dim=c(n, m, N))
        for (r in 1:N) {
            temp <- C %*% matrix(rnorm(n * (m + mtrunc)), nrow=n)
            tmp <- genvar(A, temp, mtrunc)
            E[,,r] <- tmp$E
            X[,,r] <- tmp$X
        }
    } else {  # single trial
        temp <- t(C) %*% matrix(rnorm(n * (m + mtrunc)), nrow=n)
        

        tmp <- genvar(A, temp, mtrunc)
        E <- tmp$E
        X <- tmp$X
    }

    return(list(X=X, E=E, mtrunc=mtrunc))
}

var_specrad <- function(A, newrho=NULL) {
    A_shape <- dim(A)
    if (length(A_shape) == 2) {
        A <- array(A, dim=c(A_shape[1], A_shape[2], 1))
    }
    A_shape <- dim(A)
    n <- A_shape[1]
    n1 <- A_shape[2]
    p <- A_shape[3]

    
    stopifnot(n1 == n)


    pn1 <- (p-1)*n

    if (pn1 != 0) {
        A1 <- cbind(c(A, nrow=n, ncol=p*n), diag(pn1), matrix(0, nrow=pn1, ncol=n))
    } else {
        A1 <- matrix(c(A), nrow=n, ncol=p*n)
    }

    # calculate spectral radius
    eig_A <- eigen(A1)$values
    rho <- max(abs(Re(eig_A)))

    if (is.null(newrho) | length(newrho) == 0) {
        return(rho)
    } else {
        return(list(var_decay(A, newrho/rho), rho))
    }
}

var_decay <- function(A, dfac) {
    p <- dim(A)[3]  # Assuming A is a 3D array
    f <- dfac
    for (k in 1:p) {
        A[,,k] <- A[,,k] * f
        f <- f * dfac
    }
    return(A)
}


genvar <- function(A, E, trunc=0) {
    stopifnot(is.numeric(A) & is.numeric(E) & (length(dim(A)) == 2 | length(dim(A)) == 3))
    stopifnot(is.numeric(E) & length(dim(E)) == 2)

    if (length(dim(A)) == 2) {
        A <- array(A, dim=c(dim(A)[1], dim(A)[2], 1))
    }

    dim_A <- dim(A)
    n1 <- dim_A[1]
    n2 <- dim_A[2]
    p <- dim_A[3]

    dim_E <- dim(E)
    n <- dim_E[1]
    m <- dim_E[2]

    stopifnot(n1 == n2)
    stopifnot(n1 == n)
    stopifnot(trunc >= 0 & trunc < m)

    X <- E  # Initialize X as a copy of E

    for (t in (p+1):m) {
        for (k in 1:p) {
            X[,t] <- X[,t] + A[,,k] %*% X[,t-k]  # Update X using the VAR coefficients
        }
    }

    if (trunc >= 0 & trunc < m) {
        X <- X[, (trunc + 1):m]  # Truncate X
        if (length(dim(E)) == 2) {
            E <- E[, (trunc + 1):m]  # Truncate E if it's a matrix
        }
    }

    return(list(X=X, E=E))
}


gendata <- function(AT,NT,k) {
    j <- 1
    CE_3d <- CX_3d <- NIF_3d <- IFs_3d <- p_3d <- SEIF_3d <- array(NA, dim=c(201, 2, 2))

    pb <- txtProgressBar(min = -0.5, max = 0.5, style = 3)  
    for (i in seq(-0.5, 0.5, by=0.01)) {
        setTxtProgressBar(pb, i)  
        SIFT <- matrix(c(0.5, i, i, 0.5), nrow=2, byrow=TRUE)
        result <- var_to_tsdata(AT, SIFT, NT)

        CE_3d[j,,] <- cor(t(result$E))
        CX_3d[j,,] <- cor(t(result$X))

        cau <- multi_causality_est(result$X)

        NIF_3d[j,,] <- cau$nIF
        p_3d[j,,] <- cau$p
        SEIF_3d[j,,] <- cau$SEIF
        IFs_3d[j,,] <- cau$IF

        j <- j + 1
}
close(pb)

# save to .mat
writeMat(paste("case7_data_", as.character(k), ".mat",sep=""),CX = CX_3d, CE = CE_3d, NIF = NIF_3d, SEIF = SEIF_3d, IF = IFs_3d, P=p_3d)

}


NT <- 100000

# Case 1
AT1 <- matrix(c(0.5, -0.5, -0.0, 0.6), nrow=2, byrow=TRUE)

# Case 2
AT2 <- matrix(c(-0.5, 0.5, -0.2, 0.5), nrow=2, byrow=TRUE)

# Case 3
AT3 <- matrix(c(0.5, -0.2, -0.5, 0.25), nrow=2, byrow=TRUE)

# Case 4
AT4 <- matrix(c(0.25, -0.1, -0.2, 0.1), nrow=2, byrow=TRUE)

AT <- AT1
gendata(AT1,NT,1)
gendata(AT2,NT,2)
gendata(AT3,NT,3)
gendata(AT4,NT,4)



# 创建绘图窗口
par(mfrow=c(2,2))
dev.new(width=12, height=6, unit='in', res=100) 

text_xy <- c(-0.95, 0.7)
texts <- c('a)', 'b)', 'c)', 'd)')
ylim <- list(c(10**(-7), 10**0.5), c(10**(-4), 10**0.2), c(10**(-4), 10**0.2), c(10**(-4), 10**0.2))

# normalized IF
for (i in 1:4) {
  data1 <- readMat(paste('data', i, '.mat', sep=''))
  C <- data1$CE
  x <- C[1, ]
  NIF <- data1$NIF
  y1 <- abs(NIF[1, ])
  y2 <- abs(NIF[, 1])

  plot_ind <- i*10 + 1 

  # 绘图
  plot(x, y1, type='l', col='blue', lwd=1.5, ylim=ylim[[i]], xlab='', ylab='', main='', log='y')
  lines(x, y2, col='red', lwd=1.5)
  legend('bottomleft', legend=c(expression(paste(tau[21])), expression(paste(tau[12]))), col=c('blue', 'red'), lwd=2)
  axis(1, pos=axex(origin=ylim[i][1], lengths=ylim[i][2]))
  axis(2, las=2)
  text(text_xy[1], text_xy[2], texts[i], col='black', cex=1.4, font=2)
}

library(ggplot2)
# 创建一个函数来绘制每个子图  
plot_function <- function(data, text) {  
  ggplot(data, aes(x = x)) +  
    geom_line(aes(y = y1), color = "blue", alpha = 0.7, size = 1.5) +  
    geom_line(aes(y = y2), color = "red", alpha = 0.7, size = 1.5) +  
    scale_y_log10() +  
    labs(title = text) +  
    xlim(c(-1, 1)) +  
    #ylim(c(data$ylim_lower[1], data$ylim_upper[1])) +  
    theme_minimal() +  
    theme(legend.position = "bottomleft") +  
    guides(color = guide_legend(title = NULL)) #+  
    #annotate("text", x = -0.95, y = 0.7, label = text, size = 6, fontface = "bold", color = "black")  
}  
  

# 创建绘图窗口
par(mfrow=c(2,2))
dev.new(width=12, height=6, unit='in', res=100) 

text_xy <- c(-0.95, 0.7)
texts <- c('a)', 'b)', 'c)', 'd)')
ylim <- list(c(10**(-7), 10**0.5), c(10**(-4), 10**0.2), c(10**(-4), 10**0.2), c(10**(-4), 10**0.2))




# normalized IF
par(mfrow=c(2, 2))
for (i in 1:4) {
  data1 <- readMat(paste("case7_data_", as.character(i), ".mat",sep=""))
  C <- data1$CE
  x <- C[1:101,1,2 ]
  NIF <- data1$NIF
  y1 <- abs(NIF[1:101,1,2])
  y2 <- abs(NIF[1:101,2,1])
  data <- data.frame(x = x, y1 = y1, y2 = y2, ylim_lower = 10^(-7), ylim_upper = 10^(0.5))
p <-  plot_function(data, texts[i])

ggsave(paste("case7_nIF_", as.character(i),".png",sep=''), plot = p, width = 11.11, height = 8.33, dpi = 600)
 
}