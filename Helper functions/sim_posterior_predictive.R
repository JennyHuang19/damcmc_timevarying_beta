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
sim_posterior_predictive <- function(
    S0, I0 = 1e1, t_end,
    theta = list(R0 = R0, gamma = 0.8, lambda = 1, shape = 1),
    change_day, # beta changes with each integer.
    iota_dist = "exponential",
    gener = FALSE, b = 1/2,
    E0 = 0, type = "SIR"
) {

  #
  # Population

  N0 <- S0 + I0

  #
  # Retrieve Parameters

  theta <- complete_theta(theta, iota_dist, S0)

  beta    <- theta[["beta"   ]] # revised: vector.
  gamma   <- theta[["gamma"  ]]
  lambda  <- theta[["lambda" ]]
  shape   <- theta[["shape"  ]]
  epsilon <- theta[["epsilon"]]


  #
  # Initialize

  tau_T <- tau_F <- tau_J <- rep(Inf, N0)
  S <- E <- I <- t <- W <- X <- I_tau_t_true <- c() # S, E, I, R, time, waiting time, type of event, and number of infectious at infection times

  # Initialize tau's
  tau_T[1 : I0] <- 0
  iotas <- simulate_iota(I0, iota_dist, gamma, shape, lambda) # time until recovery.
  tau_J[1 : I0] <- tau_T[1 : I0] + iotas


  # Initialize compartments, time, event type and event number
  S[1] <- S0
  E[1] <- E0
  I[1] <- I0
  t[1] <- 0
  X[1] <- "no event"
  n_t  <- n_f <- n_j <- 0

  # Keep track of mu_j
  mu_j_list = c()

  # keep track of new infections.

  # keeps track of cumulative cases at each time step
  cumu_cases = c()
  cumu_cases_time = c() # timesteps at which I report cumulative counts of cases.

  #
  # Simulation

  j <- 1 # iteration

  ### revised: SET INITIAL BETA TO beta[1].
  beta_curr <- beta[1]
{
  repeat{

    # Compute next removal time
    not_recovered <- tau_J > t[j]
    tau_J_next    <- min(tau_J[not_recovered])

    # Generate candidate infection time
    tau_T_candidate <- if(S[j] > 0) {

      ### Index to retrieve the current beta value; beta_curr.
      # take the floor of the current event time to be the index_for_beta.
      if (t[j] >= 1 & t[j] < t_end) {
        index_for_beta = floor(t[j])
        beta_curr = beta[index_for_beta+1] # 1.5 -> 1 -> 2nd beta.
      } else if (t[j] < 1) {
        beta_curr <- beta[1] # at the beginning of the epidemic, set beta to be beta[1].
      } else
        beta_curr <- beta[t_end] # beyond the observed data portion of the epidemic, set beta to be beta[t_end].

      mu_j <- beta_curr * S[j] * I[j]

      # debugging purposes: mu_j discrepancy
      mu_j_list <- append(mu_j_list, mu_j)

      if(j %% 100  == 0){
        print(mu_j)
        }

      # if(j %% 100 == 0){
      #   print(list("Iteration" = j, "mu_j" = mu_j, "beta_curr" = beta_curr,  "S(t)" = S[j], "I(t)" = I[j]))
      # }



      t[j] + rexp(1, mu_j)            # candidate infection time (##propose infection time.##)

    } else if(S[j] == 0) {            # susceptible population depleted
      Inf
    }

    # Event: removal or infection
    if(tau_T_candidate < tau_J_next) { # infection

      I_tau_t_true <- c(I_tau_t_true, I[j]) # sanity check for f_log()

      n_t <- n_t + 1
      X[j + 1] <- "t"
      t[j + 1] <- tau_T_candidate
      S[j + 1] <- S[j] - 1
      I[j + 1] <- I[j] + 1
      tau_T[I0 + n_t] <- tau_T_candidate # fill in the infection time, tau_T, for the I0+n_Tth infected individual.

      # simulate removal time of newly infected
      iota <- simulate_iota(1, iota_dist, gamma, shape, lambda)

      tau_J[I0 + n_t] <- tau_T[I0 + n_t] + iota

    } else if(tau_J_next <= tau_T_candidate) { # removal

      n_j <- n_j + 1
      X[j + 1] <- "j"
      t[j + 1] <- tau_J_next
      S[j + 1] <- S[j]
      I[j + 1] <- I[j] - 1

      if(I[j + 1] == 0)  break # infectious population depleted

    } # end-if

    # cum_cases <- append(cum_cases, mu_j)
    if(t[j] <= t_end){
      n_t_current    <- sum(0 < tau_T & is.finite(tau_T)) # is tau_T the list infection events (TRUE/FALSE) which happened before t[j]?
      cumu_cases_time <- append(cumu_cases_time, t[j])
      cumu_cases <- append(cumu_cases, n_t_current)
    }

    j <- j + 1

  } # end-repeat
}
  # Events observed before t_end (why do we only observe to t_end for the calculation of the MLE?)

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


  # MLE FOR ENTIRE TRAJECTORY
  n_t_obs    <- sum(0 < tau_T & is.finite(tau_T))
  n_j_obs    <- sum(is.finite(tau_J))

  dt              <- diff(c(t[t_obs], t_end))
  integral_I_obs  <- sum(I[t_obs] * dt)
  integral_SI_obs <- sum(I[t_obs] * S[t_obs] * dt)


  beta_MLE   <- n_t_obs / integral_SI_obs
  gamma_MLE  <- n_j_obs / integral_I_obs
  R0_MLE     <- S0 * beta_MLE / gamma_MLE


  # Output
  out <- list(
    t = t, I = I,
    final_size = n_t_obs, # the cumulative number of infections.
    beta_MLE = beta_MLE, R0_MLE = R0_MLE, cumu_cases=cumu_cases, cumu_cases_time=cumu_cases_time,
    mu_j_list = mu_j_list
  )

  return(out)

}
