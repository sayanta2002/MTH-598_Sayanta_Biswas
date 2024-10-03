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


MLE_start_1 <- function(vector,alpha) 
{
  beta <- 1
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
  
  ####### here starts g1 function
  new <- vector[2:(length(vector)-1)]
  p1 <- sum(x1^alpha) + 2* sum(x1_p^alpha) + 2 * sum(x2^alpha) + sum(x2_p^alpha) 
  p2 <- 3 * sum(x0^alpha) - 2 * sum(new^alpha)
  g1 <- p1 + p2
  
  ######## lambda_alpha
  
  lam_alpha <- (n1 + n2 +1) / g1
  
  ######## h_alpha
  
  h <- (n1 + n2 + 1) * log(alpha) + (n1 + n2 +1) * lam_alpha + 
    alpha * (log(vector[1]) + sum(log(un_x1_p_x2_p)))
  
  return(list(x1 = x1, x2 = x2, x0 = x0, x1_p = x1_p, 
              x2_p = x2_p, union_plus = un_x1_p_x2_p,
              n0 = n0, n1 = n1, n2 = n2, lam_hat_alpha = lam_alpha,
              g1 = g1, alpha = alpha, h_alpha = h))

}




alpha <- seq(1,4, 0.001)
expand <- expand.grid(alpha = alpha)
grid_alpha <- as.data.frame(expand)


iterate_MLE_dataframe <- function(vector, dataframe) {
  mega_list <- list()
  
  for (i in 1:nrow(dataframe)) {
    row <- dataframe[i, ]
    result <- MLE_start_1(vector, alpha = row)
    mega_list[[i]] <- result
  }
  
  return(mega_list)
}


# This is just an example. Replace with your actual data.
s <- synthetic_data_generator(lam0 = 1.2, lam1 = 0.9, alpha = 2 , n = 62)


mega_results <- iterate_MLE_dataframe(s$vec1, grid_alpha)

# Updated function to extract key results into a data frame, now including lam1hat
extract_key_results <- function(mega_list) {
  result <- lapply(seq_along(mega_list), function(i) {
    list(
      row = i,
      alpha = mega_list[[i]]$alpha,  # We take alpha from the original dataframe grid_alpha
      n0 = mega_list[[i]]$n0,
      n1 = mega_list[[i]]$n1,
      n2 = mega_list[[i]]$n2,
      g2 = mega_list[[i]]$g1,
      lam_profile = mega_list[[i]]$lam_hat_alpha,
      h = mega_list[[i]]$h_alpha
    )
  })
  
  return(do.call(rbind, lapply(result, data.frame)))
}

results_df <- extract_key_results(mega_results)
plot(x = results_df$alpha, y = results_df$h)

# Find the row where lambda is maximum
max_likelihood_row <- results_df[which.max(results_df$h), ]
max_likelihood_row


