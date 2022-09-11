#' Compute the sufficient statistics from the latent space
#'
#' These sufficient statistics are used to conduct inference on the parameters
#'
#' @inheritParams rprop_x
#'
#' @param return_SI logical; whether to return the trajectories of S and I
#'
#' @return sufficient statistics
#' @export
#'
#'
suff_stat <- function (x, Y, gener, b, return_SI = FALSE) {
  
  stopifnot(is.logical(gener), is.numeric(b), b > 0, is.logical(return_SI))
  
  if (!x[["compatible"]]) 
    return(list(compatible = FALSE))
  
  tau_T <- x[["tau_T"]]
  tau_J <- x[["tau_J"]]
  t_end <- Y[["t_end"]]
  I0 <- Y[["I0"]]
  S0 <- Y[["S0"]]
  
  infected <- is.finite(tau_T)
  infected_during <- infected & tau_T > 0
  recovered <- is.finite(tau_J)
  infectious <- infected & (!recovered)
  n_T <- sum(infected_during)
  n_J <- sum(recovered)
  iota_removed <- tau_J[recovered] - tau_T[recovered]
  iota_infectious <- t_end - tau_T[infectious]
  
  if (any(iota_removed < 1e-12)) 
    return(list(compatible = FALSE))
  
  tau_T <- tau_T[infected_during]
  tau_J <- tau_J[recovered]
  tau <- c(tau_T, tau_J)
  order_tau <- order(tau)
  chi <- c(rep(TRUE, n_T), rep(FALSE, n_J))[order_tau]
  
  delta_I <- ifelse(chi, 1, -1)
  delta_S <- ifelse(chi, -1, 0)
  I_tau <- c(I0, I0 + cumsum(delta_I))
  S_tau <- c(S0, S0 + cumsum(delta_S))
  I_tau_T <- I_tau[c(chi, FALSE)]
  S_tau_T <- S_tau[c(chi, FALSE)]
  
  if (any(I_tau_T == 0)) 
    return(list(compatible = FALSE))
  
  tau <- tau[order_tau]
  dtau <- diff(c(0, tau, t_end))
  integral_SI <- if (gener) 
    sum(dtau * I_tau * S_tau^(1 - b))
  else sum(dtau * I_tau * S_tau)
  integral_I <- sum(dtau * I_tau)
  
  if (return_SI) 
    return(list(n_T = n_T, n_J = n_J, S = S_tau, I = I_tau, t = c(0, tau)))
  
  SS <- list(compatible = TRUE, n_T = n_T, n_J = n_J, iota_removed = iota_removed, 
             iota_infectious = iota_infectious, integral_SI = integral_SI, 
             integral_I = integral_I, I_tau_T = I_tau_T, S_tau_T = S_tau_T)
  return(SS)
}
