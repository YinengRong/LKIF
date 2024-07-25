source("LK_Info_Flow.R")
# Set the random seed
#set.seed(0)

# Generate dataset
  

A <- matrix(c(0, -0.5, 0, 0, 0, 0,  
             0, 0, 0.7, 0, 0, 0,  
            -0.6, -0.0, -0.6, -0.0, 0, 0,  
             0, 0, 0, 0.7, 0.2, 0,  
             0.0, 0.0, -0.0, 0.4, 0.0, 0.0,  
             0, 0.8, 0, 0, 0.7, -0.5), nrow=6, byrow=TRUE)  
b <- c(0.1, 0.7, 0.5, 0.2, 0.8, 0.3)  
  
M <- 6  
xx <- matrix(0, nrow=100001, ncol=M)  
xx[1,] <- c(0.4, 0.5, 0.6, 0.7, 0.6, 0.7)  
  
for (i in 2:100001) {  
  xx[i,] <- A %*% xx[i-1,] + b + rnorm(M)  
}  


# Initialization
Nxx <- dim(xx)
xx <- t(xx[10001:nrow(xx), ])


# Calculate the causality
#cau2 <- multi_causality_est(X = xx)
# write.table(xx, file = "case2_data.txt", sep = "\t", col.names = F, row.names = F)
X <- as.matrix(read.table("case2_data.txt", header = FALSE))
cau2 <- multi_causality_est(X = t(X))


causal_graph(causal_matrix = cau2$nIF, f_name = 'case2', plot_style='kk')

nIF=cau2$nIF[,,1]
IF=cau2$IF[,,1]
err=cau2$err_e99[,,1]
p=cau2$p[,,1]


# Qualitative results
for (i in 1:Nxx[2]) {
  for (j in 1:Nxx[2]) {
    if (abs(IF[i,j]) > err[i,j]) {
      IF[i,j] <- 1
    }
    else {
        IF[i,j] <- 0
    }
  }
}

# Print results xj -> xi
cat('xj -> xi:\n')
f <- '    i'
for (i in 1:Nxx[2]) {
  f <- paste0('    ', f)
}
cat(f, '\n')
f <- '         1'
for (i in 1:(Nxx[2]-1)) {
  f <- paste0(f, sprintf('%5d', i+1))
}
cat(f, '\n')
for (j in 1:Nxx[2]) {
  if (j == floor(Nxx[2] / 2)) {
    f <- ' j'
  } else {
    f <- '  '
  }
  f <- paste0(f, sprintf('%5d', j))
  for (i in 1:Nxx[2]) {
    f <- paste0(f, sprintf('  %d  ', IF[j,i]))
  }
  cat(f, '\n')
}

# quantitative causality

# Print results xj -> xi
cat('Tj -> i:\n')
f <- '    i'
for (i in 1:Nxx[2]) {
  f <- paste0('    ', f)
}
cat(f, '\n')
f <- '            1'
for (i in 1:(Nxx[2]-1)) {
  f <- paste0(f, sprintf('%10d', i+1))
}
cat(f, '\n')
for (j in 1:Nxx[2]) {
  if (j == floor(Nxx[2] / 2)) {
    f <- ' j'
  } else {
    f <- '  '
  }
  f <- paste0(f, sprintf('%4d', j))
  for (i in 1:Nxx[2]) {
    f <- paste0(f, sprintf(' %8.4f ', nIF[j,i]))
  }
  cat(f, '\n')
}

# significant test (err99)

cat('e99:\n')
f <- '    i'
for (i in 1:Nxx[2]) {
  f <- paste0('    ', f)
}
cat(f, '\n')
f <- '            1'
for (i in 1:(Nxx[2]-1)) {
  f <- paste0(f, sprintf('%10d', i+1))
}
cat(f, '\n')
for (j in 1:Nxx[2]) {
  if (j == floor(Nxx[2] / 2)) {
    f <- ' j'
  } else {
    f <- '  '
  }
  f <- paste0(f, sprintf('%4d', j))
  for (i in 1:Nxx[2]) {
    f <- paste0(f, sprintf(' %8.4f ', err[j,i]))
  }
  cat(f, '\n')
}

#significant test (p-value)

cat('p:\n')
f <- '    i'
for (i in 1:Nxx[2]) {
  f <- paste0('    ', f)
}
cat(f, '\n')
f <- '            1'
for (i in 1:(Nxx[2]-1)) {
  f <- paste0(f, sprintf('%10d', i+1))
}
cat(f, '\n')
for (j in 1:Nxx[2]) {
  if (j == floor(Nxx[2] / 2)) {
    f <- ' j'
  } else {
    f <- '  '
  }
  f <- paste0(f, sprintf('%4d', j))
  for (i in 1:Nxx[2]) {
    f <- paste0(f, sprintf(' %8.4f ', p[j,i]))
  }
  cat(f, '\n')
}
