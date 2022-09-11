#' Log proposal density for the latent data for time-varying beta
#' (the change-point version keeps track of values of beta at each observation time T_k).
#'
#' @inheritParams rprop_x
#'
#' @param beta_vector vector of beta where each component, beta_i, corresponds to beta for a time unit (such as a week).
#' @param theta initial value for the parameters
#' @param Y observed data
#' @param x current configuration of the latent data
#' @param i_update index set of particles whose values are updated
#'
#' @return log density
#' @export
#'
b_dprop_x <- function(beta_vector,
                      theta, Y, x, i_update,
                      gener, b, iota_dist, approx
) {

  # Setup

  T_k    <- Y[["T_k"  ]] # does not count the observation at time = 0.
  I0     <- Y[["I0"   ]]
  ts     <- Y[["ts"   ]]
  t_end  <- Y[["t_end"]]
  K      <- length(T_k)
  T_sum  <- sum(T_k) + I0


  ###### for time-varying beta: calculate values of beta at each observation time T_k
  beta <- numeric(K) # vector of beta_i's

  for(k in 1 : K){
    wk = floor( k / (length(T_k) / length(beta_vector)) ) + 1 # compute the beta_i in which we observed T_k[k] (+1 for indexing of vectors in R.)
    if(k == K){
      wk = floor(  k / (length(T_k) / length(beta_vector)) ) # adjustment for the last ts value. ([k+1] adjusts ts to line up with T_k)
    }
    # compute beta values at the observation times.
    beta[k] <- beta_vector[wk]
  }
  ######


  gamma  <- theta[["gamma" ]]
  lambda <- theta[["lambda"]]
  shape  <- theta[["shape" ]]
  # beta   <- theta[["beta"  ]]

  tau_T  <- x    [["tau_T"]]
  tau_J  <- x    [["tau_J"]]
  I_k    <- x    [["I_k"  ]]
  S_k    <- x    [["S_k"  ]]


  # Contribution of infections
  contribution_infection <- if(approx == "poisson") {
    0
  } else if(approx == "ldp") { # linear death process

    # exclude initially infectious particles
    i_update_tau_T <- setdiff(i_update, 1 : I0)

    # mu_k is a vector. beta is now a vector corresponding to the the beta values at each T_k.
    mu_k  <- if(gener)  beta * I_k * S_k^(-b)  else  beta * I_k
    low   <- ts[1 : K      ] # TODO: compute only once and attach to Y
    upp   <- ts[2 : (K + 1)]

    mu_i  <- rep(mu_k, T_k) # replicates mu_1 (the first element) T_1 times and so on...
    low_i <- rep(low , T_k)
    upp_i <- rep(upp , T_k)

    sum(dexp_trunc_log( # inf error in truncate log?
      tau_T[i_update_tau_T     ], mu_i [i_update_tau_T - I0],
      low_i[i_update_tau_T - I0], upp_i[i_update_tau_T - I0]
    ))
  } # end-if(approximation)

  ### INF error
  if(contribution_infection == Inf){
    print(list("tau_T" = tau_T[i_update_tau_T     ], "low_i" = low_i[i_update_tau_T - I0], "upp_i" = upp_i[i_update_tau_T - I0], "mu" = mu_i [i_update_tau_T - I0]))
    print(list("dexp_trunc_log" = dexp_trunc_log( # inf error
      tau_T[i_update_tau_T     ], mu_i [i_update_tau_T - I0],
      low_i[i_update_tau_T - I0], upp_i[i_update_tau_T - I0]
    )))
  }
  ###

  # Contribution of removal
  p_i <- compute_pi(theta, tau_T, iota_dist, t_end)

  i_update_obs     <- setdiff(i_update, which(is.infinite(tau_J)))
  i_update_not_obs <- setdiff(i_update, which(is.finite  (tau_J)))

  contribution_removal <-
    sum(
      log(1 - p_i[i_update_not_obs])
    ) +
    sum(
      log(p_i[i_update_obs]) +
        contribution_observed_removal(
          theta, tau_T, tau_J, i_update_obs, iota_dist, t_end
        )
    )


  # Output
  loglik <- contribution_infection + contribution_removal
  # print(list("contribution_infection" = contribution_infection, "contribution_removal" = contribution_removal))
  return(loglik)

}
