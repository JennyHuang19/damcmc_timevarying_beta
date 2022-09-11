#' Run a DA-MCMC to fit the stochastic SIR model to discretely observed incidence counts with the PD-SIR algorithm.
#' This version models beta as a piece-wise constant, inferring both beta and change point positions.
#'
#' @inheritParams experiment_1_proof_of_concept
#'
#' @param Y observed data
#' @param theta_0 initial value for the parameters
#' @param print_i logical; whether to print the iteration
#' @param save_SS logical; whether to save the sufficient statistics generated each iteration
#' @param par_prior parameters of the prior distributions
#'
#' @return list with the draws for the parameters, the log likelihood, acceptance rate and size of initial susceptible population
#' @export
#'

run_DAMCMC_complete_componenth6 <- function(
  theta, Y, N = 20000,
  rho = 0.5, param = "bg", approx = "ldp",
  iota_dist = "exponential", gener = FALSE, b = 1/2,
  thin = 100, plt = 1000,
  par_prior = list(
    a_beta = 0.01, b_beta = 0.01,
    a_gamma = 1, b_gamma = 1,
    a_R0 = 2, b_R0 = 2,
    a_lambda = 1, b_lambda = 1
  ),
  length_delta = 5,
  a_pi = rep(500, length_delta),
  b_pi = rep(500, length_delta),
  num_iter_pi = 20,
  x_b_ratio = 5,
  gamma=0.1
) {

  run_time <- system.time({

    # Setup
    f_save <- b_f_save <- numeric(N / thin)
    accept <- numeric(N)
    b_accept <- numeric(N)
    ss_save <- vector(mode = "list", length = x_b_ratio*N / thin)
    pi_save <- vector(mode = "list", length = N / num_iter_pi)
    dbeta_save <- beta_save <- vector(mode = "list", length = N / thin)


    # init proposal
    # Propose the initial latent path using Y and initial theta.
    SS_current <- list(compatible = FALSE)
    while(! SS_current[["compatible"]]) {
      x_current <- rprop_x(theta, Y, gener, b, iota_dist, approx, NULL, 1) # original D-A-MCMC code.
      SS_current <- suff_stat(x_current, Y, gener, b)
    }

    # Plot initial proposal
    x_current_plt <- suff_stat(x_current, Y, gener, b, return_SI = TRUE)
    plot(x = x_current_plt$t, y = x_current_plt$I, main="Initial Latent Epidemic Trajectory",
         xlab="Time", ylab="I(t)")


    # PI: Propose Pi
    dc_current <- rep(0, length_delta) # initial SS
    pi_current <- gibbs_pi6(a_pi, b_pi, dc_current, 0)

    # 1. Propose dbeta
    dbeta_current <- Rlab::rbern(length_delta, pi_current)
    dbeta_save[1] <- list(dbeta_current)
    u_current <- 1

    # 2. determine the change-point locations.
    segment_start <- cp_locations(dbeta_current, Y)$segment_start
    segment_end <- cp_locations(dbeta_current, Y)$segment_end


    # 3. calculate unsegmented SS
    SS_current <- suff_stat(x_current, Y, gener, b, return_SI = FALSE)
    # 3. calculate sufficient statistics for each segment of x_current using the change-point locations calculated.
    segmented_SS_current <- suff_stat2(segment_start, segment_end, x_current, Y, gener, b, return_SI = FALSE)


    # 4. Propose (beta | dbeta)
    beta_current <- draw_beta(dbeta_current, segmented_SS_current, iota_dist, gamma, par_prior, param, Y)
    beta_save[1] <- list(beta_current)
    unique_beta_current <- unique(beta_current)

    # 5. Calculate current log likelihood
    b_target_current     <- 0
    for(i in 1:length(unique_beta_current)){
      b_target_current <- b_target_current + loglik_infec(unique_beta_current[i], segmented_SS_current[[i]], gener, b, iota_dist = "exponential")
    }
    b_target_current <- b_target_current + loglik_remov_infectious(gamma, SS_current, gener, b, iota_dist = "exponential")



    # MH
    for(s in 2 : N) {

      # Propose (dbeta, beta | Z) ###########

      # Update sufficient statistics for Gibbs update on pi.
      dc_current <- dc_current + dbeta_current # add-on previous iteration

      # Pi Gibbs Update: update pi and reset dc_current every num_iter_pi iterations
      if(s %% num_iter_pi == 0) {
        print(list("dc_current" = dc_current))
        pi_current <- gibbs_pi6(a_pi, b_pi, dc_current, num_iter_pi)
        # save pi every num_iter_pi to save pi directly following update.
        a <- s / num_iter_pi
        pi_save[a] <- list(pi_current)
        # reset
        dc_current <- rep(0, length_delta)
      }

      # 1. Propose dbeta (LATER: change step size)
      dbeta_new <- draw_delta_bern(pi_current, dbeta_current)[[1]]
      u_new <- draw_delta_bern(pi_current, dbeta_current)[[2]]

      # 2. determine the change-point locations.
      segment_start <- cp_locations(dbeta_new, Y)$segment_start
      segment_end <- cp_locations(dbeta_new, Y)$segment_end

      # 3. calculate sufficient statistics for Z using the start and end points.
      segmented_SS_new <- suff_stat2(segment_start, segment_end, x_current, Y, gener, b, return_SI = FALSE) # segmented_SS_new because of change points.

      # 4. Propose (beta | dbeta)
      beta_new <- draw_beta(dbeta_new, segmented_SS_new, iota_dist, gamma, par_prior, param, Y)
      unique_beta_current <- unique(beta_current)
      unique_beta_new <- unique(beta_new)

      # 5. Compute the Ratio.

      # Target density
      b_target_current     <- 0
      for(i in 1:length(unique_beta_current)){
        b_target_current <- b_target_current + loglik_infec(unique_beta_current[i], segmented_SS_current[[i]], gener, b, iota_dist = "exponential")
      }
      b_target_current <- b_target_current + loglik_remov_infectious(gamma, SS_current, gener, b, iota_dist = "exponential")

      b_target_new     <- 0
      for(i in 1:length(unique_beta_new)){
        b_target_new <- b_target_new + loglik_infec(unique_beta_new[i], segmented_SS_new[[i]], gener, b, iota_dist = "exponential")
      }
      b_target_new <- b_target_new + loglik_remov_infectious(gamma, SS_current, gener, b, iota_dist = "exponential")

      # Proposal density
      b_prop_current <- f_dbeta_func(dbeta_current, pi_current, u_current) + f_beta_func(unique_beta_current, par_prior, segmented_SS_current)
      b_prop_new <- f_dbeta_func(dbeta_new, pi_current, u_new) + f_beta_func(unique_beta_new, par_prior, segmented_SS_new)

      # Priors
      b_prior_current <- prior_beta(beta_current, par_prior) + prior_dbeta(pi_current, dbeta_current)
      b_prior_new <- prior_beta(beta_new, par_prior) + prior_dbeta(pi_current, dbeta_new)

      # MH ratio
      b_R_log <- b_prior_new - b_prior_current + b_target_new - b_target_current + b_prop_current - b_prop_new
      b_R         <- min(1, exp(b_R_log))
      b_accept[s] <- stats::runif(1) < b_R

      ### Debugging Inf values
      if(is.na(b_accept[s]) ){
        print(list("b_target_new" = b_target_new, "b_prior_new" = b_prior_new, "b_prop_new" = b_prop_new,
                   "f_pi_new" = f_pi_func(Pi_new, dbeta_new, par_prior), "f_delta_new" = f_dbeta_func(dbeta_new, Pi_new, u_new),
                   "f_beta_new" = f_beta_func(unique_beta_new, par_prior, segmented_SS_new),
                   "prior_pi_new" = prior_pi(Pi_new, par_prior), "prior_delta_new" = prior_dbeta(Pi_new, dbeta_new),"prior_beta_new" = prior_beta(beta_new, par_prior)))
      }
      ###

      if(s %% 1000 == 0){
        print("Acceptance Ratio (Beta given data)")
        print(b_R)
      }


      # Accept/reject new draws
      if(b_accept[s]) {
        # update to values.
        dbeta_current <- dbeta_new
        beta_current  <- beta_new
        b_target_current <- b_target_new
        segmented_SS_current <- segmented_SS_new
        u_current <- u_new
      }
      ##############

      ##############
      # Propose (x | dbeta, beta) ########

      # 1. propose x_new
      x_new    <- b_rprop_x(beta_current, theta, Y, gener, b, iota_dist, approx, x_current, rho)
      i_update <- x_new[["i_update"]]

      # 2. the current change-point locations.
      segment_start <- cp_locations(dbeta_current, Y)$segment_start
      segment_end <- cp_locations(dbeta_current, Y)$segment_end

      # 3. calculate unsegmented sufficient statistics for x_new
      SS_new <- suff_stat(x_new, Y, gener, b, return_SI = FALSE)
      # 3. calculate sufficient statistics for x_new using the start and end points.
      segmented_SS_new_x <- suff_stat2(segment_start, segment_end, x_new, Y, gener, b, return_SI = FALSE)

      # Target density
      unique_beta_current <- unique(beta_current)

      f_target_current <- b_target_current # (beta_current, x_current)


      f_target_new     <- 0 # (beta_current, x_new)
      for(t in 1:length(unique_beta_current)){
        f_target_new <- f_target_new + loglik_infec(unique_beta_current[t], segmented_SS_new_x[[t]], gener, b, iota_dist = "exponential")
      }
      f_target_new <- f_target_new + loglik_remov_infectious(gamma, SS_new, gener, b, iota_dist = "exponential")



      # Proposal density
      f_prop_current <- b_dprop_x(beta_current, theta, Y, x_current, i_update, gener, b, iota_dist, approx) # Infinity
      f_prop_new     <- b_dprop_x(beta_current, theta, Y, x_new    , i_update, gener, b, iota_dist, approx)

      # MH ratio
      R_log     <- f_target_new - f_target_current - f_prop_new + f_prop_current

      ### Debugging Inf values
      if(is.nan(accept[s])){
        print(list("f_target_new" = f_target_new, "f_target_current" = f_target_current, "f_prop_new" = f_prop_new, "f_prop_current" = f_prop_current))
      }
      ###

      R         <- min(1, exp(R_log))

      # print(R) # NAN

      accept[s] <- stats::runif(1) < R # how do i keep track of the new indexing for the accept vector?

      # Accept/reject new draws
      if(accept[s]) {
        x_current        <- x_new
        SS_current       <- SS_new
        f_target_current <- f_target_new
        segmented_SS_current <- segmented_SS_new_x

      } # end-if

      ##############

      # Plot x accepted.
      if(s %% plt == 0){ # plot the first proposal x.
        x_current_plt <- suff_stat(x_current, Y, gener, b, return_SI = TRUE)
        plot(x = x_current_plt$t, y = x_current_plt$I, main="Latent Epidemic Trajectory",
             xlab="Time", ylab="I(t)")

      }

      # Save thinned draws for output.
      if(s %% thin == 0){
        j <- s / thin

        f_save[[j]] <- f_target_current
        b_f_save[[j]] <- b_target_current
        dbeta_save[j] <- list(dbeta_current)
        beta_save[j] <- list(beta_current)
        ss_save[j] <- list("Iteration"=segmented_SS_current) # already a list of lists

      }

    } # end-for

    # Output
    out <- list(
      segmented_SS = ss_save, pi_lst = pi_save, dbeta_lst = dbeta_save, beta_lst = beta_save, beta_loglik = b_f_save, x_loglik = f_save, theta_rate_accept = mean(b_accept), x_rate_accept = mean(accept)
    )

  }) # end-system.time
  run_time <- as.numeric(run_time)[3]
  out[["run_time"]] <-  run_time

  return(out)
}
