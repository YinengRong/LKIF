# case 1 (bivariate causality):
set.seed(0)
source("LK_Info_Flow.R") 

#generate dataset

a11 <- 0.3
a12 <- -0.4
a22 <- 0.7
a21 <- 0
b1 <- 0.5
b2 <- 0.5

x <- rep(0, 10000)
y <- rep(0, 10000)
x[1] <- 0.4
y[1] <- 0.3

#simulate the time series
for(i in 1:(length(x)-1)){
 x[i+1] <- a11 * x[i] + a21 * y[i] + b1 * rnorm(1)
 y[i+1] <- a12 * x[i] + a22 * y[i] + b2 * rnorm(1)
}

#print the structure of system
cat(sprintf("x(i+1)=%.2f * x(i) + %.2f * y(i) + %.2f W\n", a11,a21,b1))
cat(sprintf("y(i+1)=%.2f * x(i) + %.2f * y(i) + %.2f W\n", a12,a22,b2))

#initialization
#X <- matrix(c(x, y), nrow = 2)  
#
X <- rbind(x,y)

## write.table(X, file = "case1_data.txt", sep = "\t", col.names = F, row.names = F)
X <- as.matrix(read.table("case1_data.txt", header = FALSE))

X1=X[1:2,201:10000] 
#calculate the causality
causal_graph <- multi_causality_est(X1)

#information flow
IF <- causal_graph$IF[,,1]

#normalized information flow
nIF <- causal_graph$nIF[,,1]

#significant test
e99 <- causal_graph$err_e99[,,1]

if (abs(IF[1,2]) > e99[1,2]) {
  cat(paste('y -> x percent:', format(nIF[1,2]*100, nsmall = 2), '%\n'))
} else {
  cat('y not -> x\n')
}

if (abs(IF[2,1]) > e99[2,1]) {
  cat(paste('x -> y percent:', format(nIF[2,1]*100, nsmall = 2), '%\n'))
} else {
  cat('x not -> y\n')
}