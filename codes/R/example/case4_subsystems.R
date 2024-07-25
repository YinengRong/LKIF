# Replace LK_Info_Flow with your R equivalent package
source("../LK_Info_Flow.R")

#generate dataset
a11 <- -0.5; a21 <- 0.5; a31 <- 0.2; b1 <- 1.0
a12 <- 0.0; a22 <- -0.2; a32 <- -0.6; b2 <- 1.0
a13 <- -0.2; a23 <- 0.4; a33 <- -0.2; b3 <- 1.0

b11 <- -0.2; b21 <- -0.5; b31 <- 0.0; b4 <- 1.0
b12 <- 0.5; b22 <- -0.6; b32 <- 0.4; b5 <- 1.0
b13 <- -0.1; b23 <- -0.4; b33 <- -0.5; b6 <- 1.0

esp1 <- 0.5
esp3 <- 0.0

xx <- matrix(0, 20001, 6)
xx[1,] <- rep(0.5, 6)
for (i in 1:20000) {
  xx[i+1,1] <- a11*xx[i,1] + a21*xx[i,2] + a31*xx[i,3] + b1*rnorm(1)
  xx[i+1,2] <- a12*xx[i,1] + a22*xx[i,2] + a32*xx[i,3] + b2*rnorm(1)
  xx[i+1,3] <- a13*xx[i,1] + a23*xx[i,2] + a33*xx[i,3] + b3*rnorm(1) + esp3*xx[i,6]
  xx[i+1,4] <- b11*xx[i,4] + b21*xx[i,5] + b31*xx[i,6] + b1*rnorm(1) - esp1*xx[i,1]
  xx[i+1,5] <- b12*xx[i,4] + b22*xx[i,5] + b32*xx[i,6] + b2*rnorm(1)
  xx[i+1,6] <- b13*xx[i,4] + b23*xx[i,5] + b33*xx[i,6] + b3*rnorm(1)
}

#calculate the causality
ind <- c(3,6) #X0,X1,X2 subsystem A;X3,X4,X5 subsystem B
IF_g <- causality_subspace(xx=xx, ind=ind)
cat(sprintf('  TA->B: %8.4f  TB->A: %8.4f', IF_g$TAB, IF_g$TBA))