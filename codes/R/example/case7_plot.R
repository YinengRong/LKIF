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

  
 