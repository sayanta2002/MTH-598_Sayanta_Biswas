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
