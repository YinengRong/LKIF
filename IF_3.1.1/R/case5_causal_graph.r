source("LK_Info_Flow.R")
library(MASS)

#____________________________________generate data_________________________
aii <- -0.6
bii <- 0.3
aij <- 0.3

# Generate data matrix
xx <- matrix(0, nrow = 200000, ncol = 13)

xx[1:4, ] <- matrix(rnorm(4*13), nrow = 4)
for (i in 5:200000-1) {
    xx[i+1, 1] <- -aii * xx[i, 1] + bii * rnorm(1)
    xx[i+1, 2] <- aij * xx[i, 1] + aii * xx[i, 2] + aij * xx[i, 3] + aij * xx[i-2, 5] +
                  aij * xx[i, 6] + aij * xx[i, 7] + aij * xx[i-1, 9] + aij * xx[i, 13] + bii * rnorm(1)
    xx[i+1, 3] <- aii * xx[i, 3] + aij * xx[i, 4] + bii * rnorm(1)
    xx[i+1, 4] <- aii * xx[i, 4] + bii * rnorm(1)
    xx[i+1, 5] <- aii * xx[i, 5] + bii * rnorm(1)
    xx[i+1, 6] <- aij * xx[i, 2] + aii * xx[i, 6] + bii * rnorm(1)
    xx[i+1, 7] <- aii * xx[i, 7] + bii * rnorm(1)
    xx[i+1, 8] <- aij * xx[i, 7] + aii * xx[i, 8] + bii * rnorm(1)
    xx[i+1, 9] <- -aii * xx[i, 9] + bii * rnorm(1)
    xx[i+1, 10] <- aij * xx[i, 9] + aii * xx[i, 10] + bii * rnorm(1)
    xx[i+1, 11] <- -aij * xx[i, 2] + aii * xx[i, 11] + bii * rnorm(1)
    xx[i+1, 12] <- aij * xx[i, 11] + aii * xx[i, 12] + bii * rnorm(1)
    xx[i+1, 13] <- -aij * xx[i, 12] + aii * xx[i, 13] + bii * rnorm(1)
}
#write.table(xx, file = "case5_data.txt", sep = "\t", col.names = F, row.names = F)
xx <- as.matrix(read.table("case5_data.txt", header = FALSE))
xx=t(xx)
ts <- Sys.time()
IF_1 <- multi_causality_est_OLS(X = xx, max_lag = 3)
te <- Sys.time()
print(te-ts)

causal_graph(causal_matrix = IF_1$nIF, significance = IF_1$p, c_threshold = 0.001, s_threshold = 0.01, name = 'case5')