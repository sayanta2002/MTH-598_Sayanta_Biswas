library(forecast)
library(gridExtra)
library(ggplot2)


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

v <- weibull(n = 50, lambda_not = 0.98, lambda_one = 0.98, alpha = 2)


#mle_alpha_new , mle_lambda , mle_alpha_new_1 , mle_lamda_0 , mle_lamda_


plot_acf_pacf <- function(data_vector) {
  # Calculate ACF and PACF values
  acf_vals <- acf(data_vector, plot = FALSE)
  pacf_vals <- pacf(data_vector, plot = FALSE)
  
  # Convert ACF and PACF to data frames for ggplot
  acf_df <- data.frame(Lag = acf_vals$lag, ACF = acf_vals$acf)
  pacf_df <- data.frame(Lag = pacf_vals$lag, PACF = pacf_vals$acf)
  
  # PACF plot
  pacf_plot <- ggplot(acf_df, aes(x = Lag, y = ACF)) +
    geom_bar(stat = "identity", fill = "salmon", color = "salmon") +
    geom_hline(yintercept = 0, color = "black") +
    geom_hline(yintercept = c(-1.96/sqrt(length(data_vector)), 1.96/sqrt(length(data_vector))),
               color = "#CB6D51", linetype = "dashed") +
    labs(x = "Lag", y = "PACF") +
    theme_minimal(base_size = 15)
  
  # Arrange plots side by side
  grid.arrange(pacf_plot)
}

plot_acf_pacf(v)


