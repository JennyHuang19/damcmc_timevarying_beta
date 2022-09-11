#' Simulate an SEM using Gillespie's algorithm to propose wait times and types of events. Beta is a piecewise constant.
#'
#'
#' @inheritParams experiment_1_proof_of_concept
#'
#' @param E0 intial exposed population size
#' @param type c("SIR", "SEIR"); type of process to simulate
#'
#' @return a list with all types of useful objects
#' @export
#'
sim_changepoint_sem <- function(
  S0 = 1e3, I0 = 1e1, t_end = 6,
  theta = list(R0 = c(2.5, 1.5, 3.5), gamma = 1, lambda = 1, shape = 1),
  change_day = c(2,4),
  iota_dist = "exponential",
  gener = FALSE, b = 1/2,
  E0 = 0, type = "SIR"
) {

  #
  # Population

  N0 <- S0 + I0


  #
  # Parameters

  theta <- complete_theta(theta, iota_dist, S0)

  beta    <- theta[["beta"   ]] # revised: vector.
  gamma   <- theta[["gamma"  ]]
  lambda  <- theta[["lambda" ]]
  shape   <- theta[["shape"  ]]
  epsilon <- theta[["epsilon"]]


  #
  # Initialization

  tau_T <- tau_F <- tau_J <- rep(Inf, N0)
  S <- E <- I <- t <- W <- X <- I_tau_t_true <- c() # S, E, I, R, time, waiting time, type of event, and number of infectious at infection times

  # Initialize tau's
  if(type == "SIR") {

    tau_T[1 : I0] <- 0 # TODO: relax assumption that individual initially infected at 0; important for non-Markovian process
    iotas <- simulate_iota(I0, iota_dist, gamma, shape, lambda) # time until recovery.
    tau_J[1 : I0] <- tau_T[1 : I0] + iotas

  } else if(type == "SEIR") {

    if(I0 > 0) {
      tau_T[1 : I0] <- -Inf
      tau_F[1 : I0] <- 0 # TODO: relax assumption
      iotas <- simulate_iota(I0, iota_dist, gamma, shape, lambda)
      tau_J[1 : I0] <- tau_F[1 : I0] + iotas
    } # end-if(I0)
    if(E0 > 0) {
      tau_T[I0 + (1 : E0)] <- 0
      epsilons <- rexp(E0, epsilon) # could not find rexp2
      #  epsilons <- stats::rexp(E0, epsilon)
      tau_F[I0 + (1 : E0)] <- tau_T[I0 + (1 : E0)] + epsilons
    } # end-if(E0)

  }# end-if(SIR)

  # Initialize compartments, time, event type and event number
  S[1] <- S0
  E[1] <- E0
  I[1] <- I0
  t[1] <- 0
  X[1] <- "no event"
  n_t  <- n_f <- n_j <- 0


  #
  # Simulation

  j <- 1 # iteration

  ### revised: SET INITIAL BETA TO beta[1].
  beta_curr <- beta[1]
  repeat{

    # Compute next removal time
    not_recovered <- tau_J > t[j]
    tau_J_next    <- min(tau_J[not_recovered])

    # Generate candidate infection time
    tau_T_candidate <- if(S[j] > 0) {

      ### revised: IF t[j] > CHANGE-POINT, beta_curr <- beta[i+1].
      for(i in 1:length(change_day)){
        if(t[j] >= change_day[i]){ # keep adjusting until you reach to change_day that t[j] is larger than.
          beta_curr <- beta[i+1]
        }

      }

      mu_j <- if(gener)  beta_curr * S[j]^(1 - b) * I[j]  else beta_curr * S[j] * I[j]
      t[j] + rexp(1, mu_j)            # candidate infection time (##propose infection time.##)
    } else if(S[j] == 0) {            # susceptible pop depleted
      Inf
    }

    # Event: removal or infection
    if(tau_T_candidate < tau_J_next) { # infection (SIR) / exposition (SEIR)

      I_tau_t_true <- c(I_tau_t_true, I[j]) # sanity check for f_log()

      n_t <- n_t + 1
      X[j + 1] <- "t"
      t[j + 1] <- tau_T_candidate
      S[j + 1] <- S[j] - 1
      if(type == "SEIR") {
        E[j + 1] <- E[j] + 1
        I[j + 1] <- I[j]
        #tau_T[n_t] <- tau_T_candidate # check index
      } else if(type == "SIR") {
        I[j + 1] <- I[j] + 1
        tau_T[I0 + n_t] <- tau_T_candidate
      }

      # simulate removal time of newly infected
      iota <- simulate_iota(1, iota_dist, gamma, shape, lambda)

      tau_J[I0 + n_t] <- tau_T[I0 + n_t] + iota

    } else if(tau_J_next <= tau_T_candidate) { # removal

      n_j <- n_j + 1
      X[j + 1] <- "j"
      t[j + 1] <- tau_J_next
      S[j + 1] <- S[j]
      if(type == "SEIR")  E[j + 1] <- E[j]
      I[j + 1] <- I[j] - 1

      if(I[j + 1] == 0)  break # infectious population depleted
      #TODO: if SEIR, also need to have E[j + 1]] == 0

    } # end-if

    j <- j + 1

  } # end-repeat

  # Events observed before t_end

  t_obs             <- t <= t_end
  tau_T_obs         <- tau_T <= t_end
  tau_J_obs         <- tau_J <= t_end
  tau_F_obs         <- tau_F <= t_end
  tau_T[!tau_T_obs] <- Inf
  tau_J[!tau_J_obs] <- Inf
  tau_F[!tau_F_obs] <- Inf

  x <- list(
    compatible = TRUE, tau_T = tau_T, tau_J = tau_J, tau_F = tau_F
  )

  # SEGMENT each change point interval.

  t_obs_1 <- t >= 0 & t <= 2
  t_obs_2 <- t >= 2 & t <= 4
  t_obs_3 <- t >= 4 & t <= 6

  ###

  # COMPUTE MLE FOR EACH SEGMENT.
  n_t_obs_1    <- sum(tau_T >= 0 & tau_T <= 2)
  n_t_obs_2   <- sum(tau_T >= 2 & tau_T <= 4)
  n_t_obs_3    <- sum(4 <= tau_T & tau_T <= 6)

  dt_1     <- diff(c(t[t_obs_1], 2))
  integral_I_obs_1  <- sum(I[t_obs_1] * dt_1)
  integral_SI_obs_1 <- if(gener) {
    sum(I[t_obs_1] * S[t_obs_1]^(1 - b) * dt_1)
  } else{
    sum(I[t_obs_1] * S[t_obs_1]         * dt_1)
  }

  dt_2     <- diff(c(t[t_obs_2], 4))
  integral_I_obs_2  <- sum(I[t_obs_2] * dt_2)
  integral_SI_obs_2 <- if(gener) {
    sum(I[t_obs_2] * S[t_obs_2]^(1 - b) * dt_2)
  } else{
    sum(I[t_obs_2] * S[t_obs_2]         * dt_2)
  }

  dt_3     <- diff(c(t[t_obs_3], t_end))
  integral_I_obs_3  <- sum(I[t_obs_3] * dt_3)
  integral_SI_obs_3 <- if(gener) {
    sum(I[t_obs_3] * S[t_obs_3]^(1 - b) * dt_3)
  } else{
    sum(I[t_obs_3] * S[t_obs_3]         * dt_3)
  }

  beta_1_MLE   <- n_t_obs_1 / integral_SI_obs_1
  beta_2_MLE   <- n_t_obs_2 / integral_SI_obs_2
  print(list("n_t_obs_3" = n_t_obs_3, "integral_SI_obs_3" = integral_SI_obs_3))
  beta_3_MLE   <- n_t_obs_3 / integral_SI_obs_3

  # MLE FOR ENTIRE TRAJECTORY
  n_t_obs    <- sum(0 < tau_T & is.finite(tau_T))
  n_j_obs    <- sum(is.finite(tau_J))

  dt              <- diff(c(t[t_obs], t_end))
  integral_I_obs  <- sum(I[t_obs] * dt)
  integral_SI_obs <- if(gener) {
    sum(I[t_obs] * S[t_obs]^(1 - b) * dt)
  } else{
    sum(I[t_obs] * S[t_obs]         * dt)
  }


  beta_MLE   <- n_t_obs / integral_SI_obs
  gamma_MLE  <- n_j_obs / integral_I_obs
  R0_MLE     <- S0 * beta_MLE / gamma_MLE


  # Output
  out <- list(
    x = x, t = t, X = X, S = S, E = E, I = I, t_end = t_end, I0 = I0, S0 = S0,
    MLE = c("beta" = beta_MLE, "gamma" = gamma_MLE, "R0" = R0_MLE,
            "beta_1_MLE" = beta_1_MLE, "beta_2_MLE" = beta_2_MLE, "beta_3_MLE" = beta_3_MLE),
    SS = c("n_t" = n_t_obs, "n_j" = n_j_obs, "integral_SI" = integral_SI_obs, "integral_I" = integral_I_obs,
           "n_t_3" = n_t_obs_3,  "integral_SI_3" = integral_SI_obs_3, "integral_I_3" = integral_I_obs_3)
  )

  if(iota_dist == "weibull") {

    lambda_MLE <- n_j_obs / sum( (pmin(tau_J, t_end) - tau_T)[is.finite(tau_T)]^shape )
    out[["MLE"]][["lambda"]] <- lambda_MLE

    # shape_MLE <- optimize(
    # nu_posterior,
    # interval = c(0.01, 1e1),
    # theta = theta_true, SS = SS_true,
    # maximum = TRUE
    # )[["maximum"]]
    # TODO: run optim(), with true values as initial values to obtain MLE

  }

  return(out)

}
