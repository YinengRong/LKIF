source("../LK_Info_Flow.R")

a11 <- 0.3; a12 <- -0.4; a22 <- 0.7; a21 <- 0; b1 <- 0.4; b2 <- 0.5
x <- numeric(100000)
y <- numeric(100000)
x[1] <- 0.4
y[1] <- 0.3

for (i in 1:(50000-1)) {
  x[i+1] <- a11*x[i] + a21*y[i] + b1*runif(1, -1, 1)
  y[i+1] <- a12*x[i] + a22*y[i] + b2*runif(1, -1, 1)
}

cat(paste("x(i+1)=", format(a11, nsmall = 2), "* x(i) +", format(a21, nsmall = 2), "* y(i) +", format(b1, nsmall = 2), "W\n"))
cat(paste("y(i+1)=", format(a12, nsmall = 2), "* x(i) +", format(a22, nsmall = 2), "* y(i) +", format(b2, nsmall = 2), "W\n"))

a11 <- 0.3; a12 <- 0; a22 <- 0.7; a21 <- 0.3; b1 <- 0.4; b2 <- 0.5
for (i in (50000-1):(100000-1)) {
  x[i+1] <- a11*x[i] + a21*y[i] + b1*runif(1, -1, 1)
  y[i+1] <- a12*x[i] + a22*y[i] + b2*runif(1, -1, 1)
}

cat(paste("x(i+1)=", format(a11, nsmall = 2), "* x(i) +", format(a21, nsmall = 2), "* y(i) +", format(b1, nsmall = 2), "W\n"))
cat(paste("y(i+1)=", format(a12, nsmall = 2), "* x(i) +", format(a22, nsmall = 2), "* y(i) +", format(b2, nsmall = 2), "W\n"))

cat("generating data:\n")
window_size <- 1000
T <- numeric(100000)
E99 <- numeric(100000)
T1 <- numeric(100000)
E991 <- numeric(100000)

pb <- txtProgressBar(min = 10000, max = 90000, style = 3)  
for (i in 10000:90000) {
  setTxtProgressBar(pb, i)  
  tmp <- matrix(0, nrow=2, ncol=window_size*2+1)
  tmp[1, ] <- x[(i-window_size):(i+window_size)]
  tmp[2, ] <- y[(i-window_size):(i+window_size)]
  cau <- multi_causality_est(X=tmp)
  T21 <- cau$IF[,,1]
  e99 <- cau$err_e90[,,1]
  T[i] <- T21[1, 2]
  E99[i] <- e99[1, 2]
  T1[i] <- T21[2, 1]
  E991[i] <- e99[2, 1]
}
close(pb)

library(ggplot2)
data <- data.frame(x=seq(20000, 80000), T=T1[20000:80000], E99=E991[20000:80000],y=rep(0, 60001))
p <- ggplot(data, aes(x=x, y=T)) + 
     geom_line(aes(x=x, y=T),color="red") +
     geom_line(aes(x=x, y=y), color="black") +
     geom_ribbon(aes(ymin=T-E99, ymax=T+E99), fill="blue", alpha=0.2) +
     xlim(c(20000, 80000)) +
     labs(title="information flow from X2 to X1 with 90 percent significant test")

p + theme_minimal() + theme(plot.title = element_text(hjust = 0.5, size = 16), plot.background = element_rect(fill = "white"), 
                            plot.margin = margin(1, 1, 1, 1, "cm"))

ggsave("case8_time_vary.png", plot = p, width = 11.11, height = 8.33, dpi = 72)
