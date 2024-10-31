stnwep22 <- function(x,B, lam0, lam1, alpha, los)  #observed x, B = Bootstrap size, rest are MLE parameters
{
  samp1 <- list()
  for(i in 1:B)
  {
    samp1[[i]] <- sort(weibull(length(x),lam0,lam1,alpha))
  }
  
  est1 <- matrix(0 , nrow = B , ncol = length(x))
  for(i in 1:B)
  {
    for(j in 1:length(x))
    {
      est1[i,j] <- samp1[[i]][j]
    }
  }
  est1_new <- numeric()
  for(k in 1:length(x))
  {
    est1_new[k] <- mean(est1[,k])
  }
  teststat1 <-numeric()
  for(l in 1:B)
  {
    teststat1[l] <- max(abs(samp1[[l]]- est1_new))
  }
  
  ord_test <- sort(teststat1)
  cn <- ord_test[floor((1-los)*B)]
  
  #teststat_org <- max(abs(sort(x)) - est1_new)
  
  
  teststat_org <- max(abs(sort(weibull(length(x), lam0, lam1, alpha))) - est1_new) ################# warninggggggggg 
  p1 <- length(which(teststat1 > teststat_org))/length(teststat1)
  out <- c(teststat_org , p1)
  
  
  # if(teststat_org > cn)
  #  {
  # print("Data is coming from WEP where lamda_0 != lambda_1")
  #}else
  #{
  #  print("Data is coming from WEP where lamda_0!= lambda_1")
  #}
  
  
  # ord_test2 <- sort(teststat2)
  # cn2 <- ord_test2[990]
  # teststat_org2 <- max(abs(samp2[[1]] - est2_new))
  # p2 <- length(which(teststat2 > teststat_org2))/length(teststat2)
  # out2 <- c(teststat_org2 , p2)
  # print(paste("Test Statistic for Model WEP(alpha,lambda,lambda) =", out[1] , " and P-value =",out[2])  )
  # print(paste("Test Statistic for Model WEP(alpha,lambda_0,lambda_1) =", out2[1] , " and P-value =",out2[2])  )
  # par(mfrow = c(1,2))
  # density(teststat1 , col = "red" , xlab = "Test Statistic" , main = expression(paste("WEP (", alpha,",",lambda,",",lambda ,") Model")) )
  # density(teststat2 , col = "blue" , xlab = "Test Statistic" , main = expression(paste("WEP (" ,alpha,",",lambda[0],",",lambda[1] ,") Model")))
  return(list(boot_statistic = ord_test, critical = cn, output = out))
}


# This returns ordered W
# Bootstrapped Cn(beta) values
# original test statistic value
# p value




weibull <- function(n,lambda_not,lambda_one,alpha)
{
  U <- runif(n+1 , 0 , 1)
  d <- numeric()
  for(i in 1:n)
  {
    d[i] <- min((((-1/lambda_not)*(log(U[i+1])))^(1/alpha)) , (((-1/lambda_one)*(log(U[i])))^(1/alpha)))
  }
  
  return(d)
}

v <- weibull(n = 75, lambda_not = 0.98, lambda_one = 0.98, alpha = 2)
sss <- stnwep22(v, B = 1000, lam0 = 0.994, lam1 = 0.994, alpha = 2.127, los = 0.05)




w <- as.data.frame(sss$boot_statistic)

ggplot(w, aes(x = w[,1])) +
  geom_histogram(aes(y = ..density..), binwidth = 0.005, fill = "#556B2F", color = "#556B2F",
                 alpha = 0.4) +
  geom_segment(aes(x = sss$critical, xend = sss$critical, y = 0, yend = 9.75), 
               color = "#8A2E2E", linetype = "solid", size = 1, alpha = 0.3) +
  geom_segment(aes(x = sss$output[1], xend = sss$output[1], y = 0, yend = 9.75), 
               color = "#4A5D99", linetype = "solid", size = 1, alpha = 0.3 ) +
  labs(
    x = expression(paste("W")),  # Example with alpha
    y = "Density"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 16, face = "bold"), 
    axis.title.x = element_text(size = 14),
    axis.title.y = element_text(size = 14),
    axis.text = element_text(size = 12),
    panel.grid.minor = element_blank()
  )















