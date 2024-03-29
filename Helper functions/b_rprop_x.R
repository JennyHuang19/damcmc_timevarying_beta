#' Generates PD-SIR process conditionally on the observed data
#'
#' @inheritParams run_DAMCMC
#'
#' @param theta current values of the parameters for the SIR
#' @param x     current configuration of the latent data
#'
#' @return latent data for the SIR
#' @export
#'
b_rprop_x <- function(
  beta_vector, theta, Y,
  gener, b, iota_dist, approx,
  x = NULL, rho = 1, gamma
) {

  # Setup

  T_k    <- Y[["T_k"  ]] # observed infectious at designated times.
  I0     <- Y[["I0"   ]]
  S0     <- Y[["S0"   ]]
  ts     <- Y[["ts"   ]] # observation times.
  t_end  <- Y[["t_end"]]

  theta <- complete_theta(theta, iota_dist, S0)


  beta <- beta_vector
  # beta   <- theta[["beta"  ]]
  # gamma  <- theta[["gamma" ]]
  lambda <- theta[["lambda"]]
  shape  <- theta[["shape" ]]

  K        <- length(T_k)
  T_cumsum <- cumsum(T_k) # total infectious up to a point.
  T_sum    <- T_cumsum[K] + I0 # number infected in total by t_end.
  I_k      <- J_k <- mu_k <- rep(NA, K)
  U_k      <- c(0, T_cumsum) # number removed from S
  S_k      <- (S0 - U_k)[1 : K] # number S
  p_i      <- rep(NA, T_sum)

  if(is.null(x)) { # for initial iteration of Markov chain
    tau_T  <- rep(Inf, T_sum)
    tau_J  <- rep(Inf, T_sum)
    rho    <- 1
  } else {
    tau_T <- x[["tau_T"]] # infection times
    tau_J <- x[["tau_J"]] # recovery times
  }

  # Indices of Updates
  i_update        <- which(as.logical(stats::rbinom(T_sum, 1, rho)))
  tau_J[i_update] <- Inf # tau_J that are updated are set to Inf by default

  # Initially Infectious
  I_k[1]     <- I0
  i_k        <- 1 : I0
  i_k_update <- intersect(i_k, i_update)
  n_k_update <- length(i_k_update)

  if(n_k_update > 0) {

    tau_T[i_k_update] <- 0 # should not be 0 for non-Markovian process

    # Compute p_i
    p_i[i_k_update] <- compute_pi(theta, tau_T[i_k_update], iota_dist, t_end) # uses gamma (the probability of being removed before t_end given an infection at time tau_T). (LATER: take in gamma.)

    # Sample particles recovering before t_end
    is_obs    <- as.logical(stats::rbinom(n_k_update, 1, p_i[i_k_update]))
    tau_J_obs <- i_k_update[is_obs] # indices of particles recovering before t_end

    # Propose tau_J (the tau_J not updated are Inf by default)
    tau_J[tau_J_obs] <- propose_tau_J(theta, tau_T[tau_J_obs], iota_dist, t_end) # propose removal times (LATER: change this to take in gamma.)

  }


  # Update trajectories of particles infected in interval k
  for(k in 1 : K) {

    ### for time-varying beta: keep track of the current week.
    wk = floor(k / (length(T_k) / length(beta) ) ) + 1 # mapping the interval k to the beta+index.

    if(k == K){
      wk = floor( k / (length(T_k) / length(beta) ) ) # adjustment for the last ts value.
    }
    ###

    if(T_k[k] > 0) {

      # Verify that I(t_k) > 0
      if(I_k[k] <= 0)  return(list(compatible = FALSE))

      i_k        <- I0 + (U_k[k] + 1) : U_k[k + 1] # indices of particles infected during interval k
      i_k_update <- intersect(i_k, i_update)
      n_k_update <- length(i_k_update) # number updated in interval k.

      if(n_k_update > 0) {

        # Propose tau_T
        mu_k[k]  <- if(gener)  beta[wk] * I_k[k] * S_k[k]^(-b)  else  beta[wk] * I_k[k] # use the corresponding beta_i.

        tau_T[i_k_update] <- propose_tau_T(n_k_update, mu_k[k], ts[k], ts[k+1], approx)

        # Propose tau_J
        p_i[i_k_update]  <- compute_pi(theta, tau_T[i_k_update], iota_dist, t_end)
        is_obs           <- as.logical(stats::rbinom(n_k_update, 1, p_i[i_k_update])) # particles recovering before t_end
        tau_J_obs        <- i_k_update[is_obs]
        tau_J[tau_J_obs] <- propose_tau_J(theta, tau_T[tau_J_obs], iota_dist, t_end = t_end)

      } # end-if(n_k_update)

    } # end-if(T_k[k])

    # Update I_k
    if(k < K) {
      J_k[k    ] <- sum(dplyr::between(tau_J, ts[k], ts[k + 1]))
      I_k[k + 1] <- I_k[k] + T_k[k] - J_k[k]
    }

  } # end-for

  # Output
  x_new <- list(
    compatible = TRUE, tau_T = tau_T, tau_J = tau_J, i_update = i_update, I_k = I_k, S_k = S_k # updated tau_T and tau_J lists.
  )
  return(x_new)

}
