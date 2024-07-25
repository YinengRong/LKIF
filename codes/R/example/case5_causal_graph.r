source("../LK_Info_Flow.R")
library(MASS)

#____________________________________generate data_________________________
# 生成数据  
n <- 20000  
p <- 5  
xx <- matrix(0, nrow = n, ncol = p)  
xx[1:4, ] <- rnorm(4 * p)  
  
for (i in 5:(n-1)) {  
  xx[i+1, 1] <- -0.95 * sqrt(2) * xx[i, 1] - 0.9025 * xx[i-1, 1] + rnorm(1)  
  xx[i+1, 2] <- 0.5 * xx[i-1, 1] + rnorm(1)  
  xx[i+1, 3] <- -0.4 * xx[i-2, 1] + rnorm(1)  
  xx[i+1, 4] <- -0.5 * xx[i-1, 1] + 0.25 * sqrt(2) * xx[i, 4] +   
                 0.25 * sqrt(2) * xx[i, 5] + rnorm(1)  
  xx[i+1, 5] <- -0.25 * sqrt(2) * xx[i, 4] + 0.25 * sqrt(2) * xx[i, 5] + rnorm(1)  
}  
#write.table(xx, file = "case5_data.txt", sep = "\t", col.names = F, row.names = F)
#xx <- as.matrix(read.table("case5_data.txt", header = FALSE))

xx1 <- xx[1000:n, ]  
xx1 <- cbind(xx1, xx[999:(n-1), 1], xx[998:(n-2), 1])  
IF_0 <- multi_causality_est_OLS(t(xx1))  
a <- IF_0$nIF  
b <- IF_0$p  
a[,6:nrow(a) , ] <- 0  
b[,6:nrow(b) ,] <- 1  

causal_graph(causal_matrix = a, significance = b, c_threshold = 0.001, s_threshold = 0.01, f_name = 'case5', plot_style='fr') 

