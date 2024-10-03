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










### The MLE function for every xi 's or others

MLE_start <- function(vector, lam0, lam1, alpha) {
  gamma <- lam0/lam1
  beta <- (gamma)^(1/alpha)
  l <- length(vector)
  
  index1 <- numeric()
  index2 <- numeric()
  index1_plus_one <- numeric()
  index2_plus_one <- numeric()
  i0_beta <- numeric()
  
  for (i in 1:(l-1)) {
    if (vector[i] * beta < vector[i+1]) {
      index1[i] <- i
      index1_plus_one[i] <- i+1
    } else if (vector[i] * beta > vector[i+1]) {
      index2[i] <- i
      index2_plus_one[i] <- i+1
    } else {
      i0_beta[i] <- vector[i]
    }
  }
  
  # Remove NA entries
  i_1 <- index1[!is.na(index1)]
  i_2 <- index2[!is.na(index2)]
  i1_plus1 <- index1_plus_one[!is.na(index1_plus_one)]
  i2_plus1 <- index2_plus_one[!is.na(index2_plus_one)]
  i0_beta <- i0_beta[!is.na(i0_beta)]
  
  union <- sort(union(i1_plus1, i2_plus1))
  
  x1 <- vector[i_1]
  x2 <- vector[i_2]
  x1_p <- vector[i1_plus1]
  x2_p <- vector[i2_plus1]
  x0 <- i0_beta
  
  un_x1_p_x2_p <- vector[union]
  
  n0 <- length(x0)
  n1 <- length(x1)
  n2 <- length(x2)
  
  ####### here starts g2 function
  new <- vector[2:(length(vector)-1)]
  t1 <- sum(x1^alpha) + (1+gamma)*((sum(x1_p^alpha) + sum(x2^alpha))) +
    gamma * sum(x2_p^alpha) + (gamma^2 + gamma + 1) * sum(x0^alpha) 
  t2 <- (1 + gamma) * sum(new^alpha)
  g2_al_gam <- t1 - t2
  
  ######## lam1_hat
  lam1hat <- (n1 + n2 + 1)/g2_al_gam
  
  ######## h2
  
  h2 <- log(vector[1]) + sum(log(un_x1_p_x2_p))
  
  ####### Here starts likelihood
  l1 <- (n1 + n2 + 1) * (log(lam1hat) + log(alpha)) - (n1 + n2 + 1) -
    (n0 - 1) * log(1 + gamma) 
  
  l2 <- ((n0 - (1/alpha)) + n2) * log(gamma) + (alpha - 1) * h2
  
  l_profiled <- l1 + l2 + (alpha - 1) * h2
  
  
  
  return(list(lamda0 = lam0, lamda1 = lam1, gamma = gamma, 
              beta = beta, 
              x1 = x1, x2 = x2, x0 = x0, 
              x1_p = x1_p, x2_p = x2_p, 
              union_plus = un_x1_p_x2_p,
              n0 = n0, n1 = n1, n2 = n2, g2 = g2_al_gam,
              lam1hat = lam1hat,
              likelihood = l_profiled))
}






# The iterate_MLE_dataframe function remains the same
iterate_MLE_dataframe <- function(vector, dataframe) {
  mega_list <- list()
  
  for (i in 1:nrow(dataframe)) {
    row <- dataframe[i, ]
    result <- MLE_start(vector, lam0 = row$lambda, lam1 = row$lam1, alpha = row$alpha)
    mega_list[[i]] <- result
  }
  
  return(mega_list)
}



mega_results <- iterate_MLE_dataframe(data_vector, grid_df)

# Updated function to extract key results into a data frame, now including lam1hat
extract_key_results <- function(mega_list) {
  result <- lapply(seq_along(mega_list), function(i) {
    list(
      row = i,
      lambda0 = mega_list[[i]]$lamda0,
      lambda1 = mega_list[[i]]$lamda1,
      alpha = grid_df$alpha[i],  # We take alpha from the original dataframe
      gamma = mega_list[[i]]$gamma,
      beta = mega_list[[i]]$beta,
      n0 = mega_list[[i]]$n0,
      n1 = mega_list[[i]]$n1,
      n2 = mega_list[[i]]$n2,
      g2 = mega_list[[i]]$g2,
      lam1hat = mega_list[[i]]$lam1hat,
      likelihood = mega_list[[i]]$likelihood
    )
  })
  
  return(do.call(rbind, lapply(result, data.frame)))
}































