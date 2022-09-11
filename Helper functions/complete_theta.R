complete_theta <- function(theta, iota_dist, S0) {

  if(iota_dist == "exponential") theta[c("shape", "lambda")] <- NULL
  if(iota_dist == "weibull")     theta <- add_gamma(theta)
  theta <- add_beta(theta, S0) # theta[["beta"]] will be c(b1, b2)

  return(theta)

}
