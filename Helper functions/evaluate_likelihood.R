#' Log likelihood (contributions of the infected).
#'
#'
#'
#'
#' @param SS sufficient statistics of the current configuration of the latent data
#'
#' @return log likelihood
#' @export
#'
loglik_infec <- function(beta_i, SS, gener, b, iota_dist = "exponential") {

  # Setup
  # beta   <- theta[["beta"  ]]
  # gamma  <- theta[["gamma" ]]
  #shape  <- theta[["shape" ]]
  #lambda <- theta[["lambda"]]

  # print(list("SS" = SS)) # SS[["compatible"]] : subscript out of bounds

  n_T             <- SS[["n_T"            ]]
  I_tau_T         <- SS[["I_tau_T"        ]]
  S_tau_T         <- SS[["S_tau_T"        ]]
  integral_SI     <- SS[["integral_SI"    ]]
  iota_removed    <- SS[["iota_removed"   ]]
  iota_infectious <- SS[["iota_infectious"]]


  # contribution of infections
  loglik_infec <- if(gener){
    n_T * log(beta_i) + sum(log(I_tau_T) - b * log(S_tau_T)) - beta_i  * integral_SI
  } else {
    n_T * log(beta_i) + sum(log(I_tau_T)) -beta_i  * integral_SI # could log(I_tau_T) be log(0)? Did this occur in the original code? Why would I(t) be 0 at the time of an infection?
  }

  ### adjustment for small segments.

  # if(is.nan(loglik_infec)){ # loglik_infec NaN
  #   loglik_infec <- 0.00001
  # }

  # ###


  return(loglik_infec)

}


#' Log likelihood (contributions of the removal) which does not account for iota_infect.
#'
#'
#'
#'
#'
#' @param SS sufficient statistics of the current configuration of the latent data
#'
#' @return log likelihood removal
#' @export
#'
loglik_remov <- function(gamma, SS, gener, b, iota_dist = "exponential") {

  n_T             <- SS[["n_T"            ]]
  I_tau_T         <- SS[["I_tau_T"        ]]
  S_tau_T         <- SS[["S_tau_T"        ]]
  integral_SI     <- SS[["integral_SI"    ]]
  iota_removed    <- SS[["iota_removed"   ]]
  iota_infectious <- SS[["iota_infectious"]]

  # contribution of removals
  loglik_remov <-
    if(iota_dist == "exponential") {
      sum(stats::dexp(iota_removed   , gamma        , log   = TRUE                    )) # removals before t_end (observed)
      # sum(stats::pexp(iota_infectious, gamma        , log.p = TRUE, lower.tail = FALSE))   # removals after  t_end (not observed) P[T > t]
    } else if(iota_dist == "weibull") {
      sum(dweibull2  (iota_removed   , shape, lambda, log   = TRUE                    ))
      # sum(pweibull2  (iota_infectious, shape, lambda, log.p = TRUE, lower.tail = FALSE))
    }

  return(loglik_remov)

}

#' Log likelihood (contributions of the removal) which does account for iota_infect.
#'
#'
#'
#'
#'
#' @param SS sufficient statistics of the current configuration of the latent data
#'
#' @return log likelihood removal
#' @export
#'
loglik_remov_infectious <- function(gamma, SS, gener, b, iota_dist = "exponential") {

  n_T             <- SS[["n_T"            ]]
  I_tau_T         <- SS[["I_tau_T"        ]]
  S_tau_T         <- SS[["S_tau_T"        ]]
  integral_SI     <- SS[["integral_SI"    ]]
  iota_removed    <- SS[["iota_removed"   ]]
  iota_infectious <- SS[["iota_infectious"]]

  # contribution of removals
  loglik_remov <-
    if(iota_dist == "exponential") {
      sum(stats::dexp(iota_removed   , gamma        , log   = TRUE                    )) + # removals before t_end (observed)
        sum(stats::pexp(iota_infectious, gamma        , log.p = TRUE, lower.tail = FALSE))   # removals after  t_end (not observed) P[T > t]
    } else if(iota_dist == "weibull") {
      sum(dweibull2  (iota_removed   , shape, lambda, log   = TRUE                    )) +
        sum(pweibull2  (iota_infectious, shape, lambda, log.p = TRUE, lower.tail = FALSE))
    }

  return(loglik_remov)

}
