source("LK_Info_Flow.R")
# Set the random seed
set.seed(0)

# Generate dataset
a <- matrix(c(0.3, 0, 0, 0,
              0.5, 0.7, 0.1, 0,
              0, 0.3, 0.5, 0,
              0.2, 0.4, 0.3, 0.1), nrow = 4, byrow = TRUE)

b <- c(0.4, 0.5, 0.6, 0.7)

xx <- matrix(0, nrow = 1001, ncol = 4)
xx[1, ] <- c(0.4, 0.5, 0.6, 0.7)

for (i in 1:1000) {
    xx[i+1, 1] <- a[1,1]*xx[i,1] + a[1,2]*xx[i,2] + a[1,3]*xx[i,3] + a[1,4]*xx[i,4] + b[1]*rnorm(1)
    xx[i+1, 2] <- a[2,1]*xx[i,1] + a[2,2]*xx[i,2] + a[2,3]*xx[i,3] + a[2,4]*xx[i,4] + b[2]*rnorm(1)
    xx[i+1, 3] <- a[3,1]*xx[i,1] + a[3,2]*xx[i,2] + a[3,3]*xx[i,3] + a[3,4]*xx[i,4] + b[3]*rnorm(1)
    xx[i+1, 4] <- a[4,1]*xx[i,1] + a[4,2]*xx[i,2] + a[4,3]*xx[i,3] + a[4,4]*xx[i,4] + b[4]*rnorm(1)
}

#print the structure of the system
cat(sprintf("x1(i+1)=%.2f * x1(i) + %.2f * x2(i) + %.2f * x3(i) + %.2f * x4(i) + %.2f W\n", a[1,1], a[1,2], a[1,3], a[1,4], b[1]))
cat(sprintf("x2(i+1)=%.2f * x1(i) + %.2f * x2(i) + %.2f * x3(i) + %.2f * x4(i) + %.2f W\n", a[2,1], a[2,2], a[2,3], a[2,4], b[2]))
cat(sprintf("x3(i+1)=%.2f * x1(i) + %.2f * x2(i) + %.2f * x3(i) + %.2f * x4(i) + %.2f W\n", a[3,1], a[3,2], a[3,3], a[3,4], b[3]))
cat(sprintf("x4(i+1)=%.2f * x1(i) + %.2f * x2(i) + %.2f * x3(i) + %.2f * x4(i) + %.2f W\n", a[4,1], a[4,2], a[4,3], a[4,4], b[4]))

# Initialization
Nxx <- dim(xx)
t <- seq(0, 1000, by = 1) # time series

# Calculate the causality
#cau2 <- multi_causality_est_OLS(X = t(xx), series_temporal_order = t)
#write.table(t(xx), file = "case2_data.txt", sep = "\t", col.names = F, row.names = F)
X <- as.matrix(read.table("case2_data.txt", header = FALSE))
cau2 <- multi_causality_est_OLS(X = X, series_temporal_order = t)

nIF=cau2$nIF[,,1]
IF=cau2$IF[,,1]
err=cau2$err_e99[,,1]
p=cau2$p[,,1]


# Qualitative results
for (i in 1:4) {
  for (j in 1:4) {
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
f <- '    j'
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
    f <- ' i'
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
f <- '    j'
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
    f <- ' i'
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
f <- '    j'
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
    f <- ' i'
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
f <- '    j'
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
    f <- ' i'
  } else {
    f <- '  '
  }
  f <- paste0(f, sprintf('%4d', j))
  for (i in 1:Nxx[2]) {
    f <- paste0(f, sprintf(' %8.4f ', p[j,i]))
  }
  cat(f, '\n')
}