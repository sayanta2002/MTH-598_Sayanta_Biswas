inv <- function(lam, alpha, u)
{
  g <- (-1/lam)*log(u)
  c0 <- g^(1/alpha)
  return(c0)
}



synthetic_data_generator <- function(lam0, lam1, alpha, n) #uniform(0,1) of size n
{
  u <- runif(n+2, 0,1)
  w <- numeric()
  for (i in 1:(length(u)-1))
  {
    
    w[i] <- min(inv(lam = lam0, alpha = alpha, u = u[i]), ###change needed
                inv(lam = lam1, alpha = alpha, u = u[i-1]))
    
  }
  vec <- as.vector(w)
  vec1 <- vec[2:length(vec)]
  #max(vec1)
  #density(vec1)
  #vec1_round <- round(vec1, digits = 2)
  time <- seq(from =1, to = length(vec1), by = 1)
  df <- as.data.frame(cbind(time, vec1))
  return(df)
} #output is a dataframe


# a vector containing WEP(lamda0, lambda1, alpha)

s <- synthetic_data_generator(lam0 = 0.15, lam1 = 0.04, alpha = 3.0, n = 75)
plot.ts(s$vec1)






### checking if the sample mean matches

meaner <- function(vector,lam0, lam1, alpha)
{
  ff <- vector
  lambda <- (lam0 + lam1)
  mean <- (lambda^(-1/alpha)) * gamma((1/alpha) + 1)
  samp <- mean(ff)
  
  return(c(samp, mean, length(vector)))
}

meaner(s$vec1, 0.15,0.04,3)


































