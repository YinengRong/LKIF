# case 3(panel data, discontinuous time series or ensemble data):
## 3.1________________________________________________
library(progress)
source("../LK_Info_Flow.R")
# Define parameters
a <- matrix(c(0.3, 0, 0, 0.4,
              0.5, 0.7, 0.1, 0.5,
              0, 0.3, 0.5, 0.6,
              0.2, 0.4, 0.3, 0.1), nrow = 4, byrow = TRUE)
b <- c(0.4, 0.5, 0.6, 0.3)
case_num <- 1000

# Generate data
set.seed(123)  # for reproducibility
xx <- array(0, dim = c(case_num, 1001, 4))
xx[,1,] <- matrix(rep(c(0.4, 0.5, 0.6, 0.7), case_num), nrow = case_num, byrow = TRUE)

pb <- txtProgressBar(min = 0, max = case_num, style = 2)  
for (icase in 1:case_num) {
  setTxtProgressBar(pb, icase)  
  if (icase > 1) {
    xx[icase,1,] <- xx[icase,1,] * rnorm(4, mean = 0, sd = 0)
  }
  for (i in 1:1000) {
    xx[icase,i + 1,1] <- sum(a[1,] * xx[icase,i,] + b[1] * rnorm(1))
    xx[icase,i + 1,2] <- sum(a[2,] * xx[icase,i,] + b[2] * rnorm(1))
    xx[icase,i + 1,3] <- sum(a[3,] * xx[icase,i,] + b[3] * rnorm(1))
    xx[icase,i + 1,4] <- sum(a[4,] * xx[icase,i,] + b[4] * rnorm(1))
  }
}
close(pb)
print('Finished generating the data.')

## Select a segment of 10 time length for each case and combine it into panel data
## and build the corresponding temporal index
X = array(0, dim=c(10*case_num,4))
t = array(0, dim=c(10*case_num))
for(j in 1:case_num){
    i=floor(runif(1)*case_num)+1
    X[(10*j-9):(10*j),]=xx[i,(991:1000),]
    t[(10*j-9):(10*j)]=0:9
}



cat('start calculate causality:\n')
# Record start time
time_start <- Sys.time()
#write.table(X, file = "case3_data_X.txt", sep = "\t", col.names = F, row.names = F)
X <- as.matrix(read.table("case3_data_X.txt", header = FALSE))

#write.table(t, file = "case3_data_t.txt", sep = "\t", col.names = F, row.names = F)
t <- as.matrix(read.table("case3_data_t.txt", header = FALSE))

# Calculate causality with temperal index (panel data) 
IF_panel <- multi_causality_est(X = t(X[,1:4]), series_temporal_order = t)
time_end <- Sys.time()

T21=IF_panel$IF[,,1]
err=IF_panel$err_e90[,,1]


# Print the estimated values and errors
cat(sprintf('est_panel: T1->2: %8.4f e90: %8.4f\n', T21[2,1], err[2,1]))
cat(sprintf('est_panel: T2->1: %8.4f e90: %8.4f\n', T21[1,2], err[1,2]))

# Calculate and print the time cost
time_end <- Sys.time()
cat(sprintf('time cost: %8.4f s\n', as.numeric(difftime(time_end, time_start, units = "secs"))))



# Calculate causality with temperal index (panel data) 
IF_panel <- multi_causality_est(X = t(X[,1:4]))
time_end <- Sys.time()

T21=IF_panel$IF[,,1]
err=IF_panel$err_e90[,,1]


# Print the estimated values and errors
cat(sprintf('est: T1->2: %8.4f e90: %8.4f\n', T21[2,1], err[2,1]))
cat(sprintf('est: T2->1: %8.4f e90: %8.4f\n', T21[1,2], err[1,2]))

# Calculate and print the time cost
time_end <- Sys.time()
cat(sprintf('time cost: %8.4f s\n', as.numeric(difftime(time_end, time_start, units = "secs"))))




## 3.2_________________________________________________________________________________________________
NT <- 12000

AT <- array(0, dim = c(3, 3, 12))
for (i in 1:3) {
  AT[i, i, ] <- 0.4
}
AT[2, 2, ] <- -0.9
AT[1, 2, 1] <- -0.1
AT[1, 2, 2] <- -0.2
AT[1, 2, 3] <- -0.1
AT[2, 1, 4] <- -0.1
AT[2, 1, 5] <- -0.2
AT[2, 1, 6] <- -0.2
AT[2, 1, 7] <- -0.3
AT[2, 1, 8] <- -0.4
AT[2, 1, 9] <- -0.3
AT[2, 1, 10] <- -0.2
AT[2, 1, 11] <- -0.2
AT[2, 1, 12] <- -0.1

#plot the cyclic element A12 and A21
x1 <- c(AT[1, 2, ])
x2 <- c(AT[2, 1, ])

plot(1:12, x1, type = "o", col = "red", pch = 19, cex = 2, lwd = 2, xlab = "Time", ylab = "Value", ylim = range(c(x1, x2)))
points(1:12, x2, type = "o", col = "blue", pch = 19, cex = 2, lwd = 2)
legend("topright", legend = c("a12", "a21"), col = c("red", "blue"), lwd = c(2, 2), pch = c(19, 19))
axis(1, at = 1:12, labels = c(1:12), cex.axis = 1.5)


# Repeat the experiment 100 times
B <- diag(3) * 0.3 + matrix(0.3, nrow = 3, ncol = 3)
X <- matrix(0, nrow = 3, ncol = NT + 1200)

NIF <- P <- SEIF <- IF <- array(NA, dim=c(100, 3, 3))
NIF1 <- P1 <- SEIF1 <- IF1 <- array(NA, dim=c(100, 3, 3))

pb <- txtProgressBar(min = 0, max = 100, style = 3)  
for (in_ in 1:100) {
  setTxtProgressBar(pb, in_)  
  ## Exclude the first 1200 time steps and select the stable time series for information flow calculation.
  vt <- 1200
  for (it in 1:(NT+vt-1)) {
    X[, it+1] <- t(AT[, , (it+2) %% 12 + 1]) %*% X[, it]+(B %*% matrix(rnorm(3), nrow = 3))
  }

  nn <- dim(X)
  N <- dim(X[, seq(from = 1201, to = ncol(X), by = 12)])
  #N <- dim(X[, (1201 + 1):ncol(X)])

  xx <- matrix(0, nrow = 3, ncol = N[2]*2)
  t <- matrix(0, nrow = 1, ncol = N[2]*2)

  # Select Feb and Mar
  for (i in 1:3) {
    part1 <- X[i, seq(1201, ncol(X), by = 12)]  
    part2 <- X[i, seq(1202, ncol(X), by = 12)] 
    
    xx[i, seq(from = 1, to = N[2]*2, by = 2)] <- part1
    xx[i, seq(from = 2, to = N[2]*2, by = 2)] <- part2
  }
  t[seq(from = 1, to = N[2]*2, by = 2)] <-seq(from = 1, to = N[2]*2*6, by = 12)
  t[seq(from = 2, to = N[2]*2, by = 2)] <-seq(from = 2, to = N[2]*2*6, by = 12)

  ## Information flow for panel data
  
  cau2 <- multi_causality_est(X = xx, series_temporal_order = t)
  NIF[in_,,] <- cau2$nIF[,,1]
   IF[in_,,] <- cau2$IF[,,1]
    P[in_,,] <- cau2$p[,,1]
  SEIF[in_,,]<- cau2$SEIF[,,1]
  ## Information flow for temporal data
  cau2 <- multi_causality_est(X = xx)
  NIF1[in_,,] <- cau2$nIF[,,1]
   IF1[in_,,] <- cau2$IF[,,1]
    P1[in_,,] <- cau2$p[,,1]
  SEIF1[in_,,]<- cau2$SEIF[,,1]
  
}
close(pb)





cat('Ground truth is: X2->X1\n')
# Average IF
cat('IF(panel)\n')
mean_IF <- apply(IF, c(2, 3), mean)
print(round(mean_IF, 2))
cat('IF(temperal)\n')
mean_IF1 <- apply(IF1, c(2, 3), mean)
print(round(mean_IF1, 2))


# average nIF
cat("normalized IF(panel)\n")
mean_NIF <- apply(NIF, c(2, 3), mean)
mean_NIF <- round(mean_NIF, 2)
print(mean_NIF)

cat("normalized IF(temporal)\n")
mean_NIF1 <- apply(NIF1, c(2, 3), mean)
mean_NIF1 <- round(mean_NIF1, 2)
print(mean_NIF1)

# average p-value
cat("p-value(panel)\n")
mean_P <- apply(P, c(2, 3), mean)
mean_P <- round(mean_P, 2)
print(mean_P)

cat("p-value(temporal)\n")
mean_P1 <- apply(P1, c(2, 3), mean)
mean_P1 <- round(mean_P1, 2)
print(mean_P1)

# average effect size
cat("effect size(panel)\n")
mean_E <- apply(abs(IF / SEIF / sqrt(1000-3)), c(2,3), mean)
mean_E <- round(mean_E,2)
print(mean_E)
key <- (mean_E[1, 2] + mean_E[2, 1]) / 2

cat("effect size(temporal)\n")
mean_E1 <- apply(abs(IF1 / SEIF1 / sqrt(2000-3-1)), c(2,3), mean)
mean_E1 <- round(mean_E1, 2)
print(mean_E1)
key1 <- (mean_E1[1, 2] + mean_E1[2, 1]) / 2

A <- array(0, dim=c(100, 3, 3))
for (i in 1:100) {
  A[i, , ] <- diag(3)
  A[i, 1, 2] <- 1
}
mean_distance_panel <- apply(abs(A - (abs(IF / SEIF / sqrt(1000-3)) > key)), c(2,3),sum)
cat("structure Hanming Distance(panel)\n")
print(mean_distance_panel)

mean_distance_temporal <- apply(abs(A - (abs(IF1 / SEIF1 / sqrt(2000-3-1)) > key1)), c(2,3),sum)
cat("structure Hanming Distance(time_series)\n")
print(mean_distance_temporal)