#' Compute the sufficient statistics from the latent space for each segment of beta.
#'
#' These sufficient statistics are used to conduct inference on parameters.
#'
#' @inheritParams rprop_x
#'
#' @param return_SI logical; whether to return the trajectories of S and I
#'
#' @return a list of lists of sufficient statistics, in which each element is a list of sufficient statistics for a segment of constant beta.
#' @export
#'

suff_stat2 <- function(segment_start, segment_end, x, Y, gener, b, return_SI = FALSE) {

  # Verify compatibility of x
  if( !x[["compatible"]] )  return(list(list("!x[[compatible]]", compatible = FALSE)))

  # list to store SS for each PW constant beta.
  SS_list <- list()

  # extract data from epidemic.
  tau_T <- x[["tau_T"]] # list of infection event times (already in order)
  tau_J <- x[["tau_J"]] # list of recovery event time (not originally in order)
  t_end <- Y[["t_end"]]
  I0    <- Y[["I0"   ]]
  S0    <- Y[["S0"   ]]

  # list of I(tau) and S(tau) at each change-point.
  I_tau_changepoints <- c(I0)
  S_tau_changepoints <- c(S0)

  for(i in 1:length(segment_start)){

    # Infected
    infected        <- is.finite(tau_T)
    infected_during <- infected & tau_T > 0 # excludes initially infectious

    infected_during_segment <- (tau_T > segment_start[i]) & (tau_T <= segment_end[i])
    infected_before_segment <- tau_T < segment_start[i]
    # Recovered
    recovered       <- is.finite(tau_J)
    recovered_during_segment <- (recovered & tau_J > segment_start[i]) & (recovered & tau_J < segment_end[i])
    recovered_before_segment <- recovered & tau_J < segment_start[i]

    # Those who remain infectious (used in f_log)
    # infectious      <- infected & (! recovered)
    infectious_after_segment <- (infected_during_segment | infected_before_segment) & (!(recovered_during_segment | recovered_before_segment))

    # Number of events
    n_T_segment <- sum(infected_during_segment)
    n_J_segment <- sum(recovered_during_segment)

    # Iota (the removal wait times)
    # iota_removed    <- tau_J[recovered] - tau_T[recovered]
    iota_removed_segment <- tau_J[recovered_during_segment] - tau_T[recovered_during_segment]

    # wait times for those individuals who remain infectious until t_end.
    # iota_infectious <- t_end            - tau_T[infectious]
    iota_infectious_segment <- segment_end[i]     - tau_T[infectious_after_segment] # Here, we count the (T_end - infection_times) for all individuals who did not recovered by T_end.

    # wait times should be positive.
    if(any(iota_removed_segment < 1e-12))  return(list(list("iota_removed_segment", compatible = FALSE)))
    if(any(iota_infectious_segment < 1e-12))  return(list(list("iota_infectious_segment", compatible = FALSE)))


    # Event time
    tau_T_segment <- tau_T[infected_during_segment] # infection times for those infected during [segment start, segment end].
    tau_J_segment <- tau_J[recovered_during_segment]

    # ALL INFECTION AND RECOVERY TIMES TOGETHER.
    tau_segment <- c(tau_T_segment, tau_J_segment)
    order_tau_segment <- order(tau_segment)

    #--- Compute S(tau), I(tau), I(tau_T)

    # type of event (TRUE: infection, FALSE: recovery)
    chi_segment     <- c(rep(TRUE, n_T_segment), rep(FALSE, n_J_segment))[order_tau_segment] # Sanity check: all.equal(X=="t", chi)

    # the waiting times between each event
    delta_I_segment <- ifelse(chi_segment,  1, -1)
    delta_S_segment <- ifelse(chi_segment,  -1, 0)

    # the number of I(tau), S(tau) at each event.
    I_tau_segment   <- c(I_tau_changepoints[i], I_tau_changepoints[i] + cumsum(delta_I_segment)) # cumsum: sums the the first 1:n terms.
    S_tau_segment   <- c(S_tau_changepoints[i], S_tau_changepoints[i] + cumsum(delta_S_segment))

    #--- Append to the list holding I(t) and S(t) at change points.
    next_I0 <- n_T_segment - n_J_segment + I_tau_changepoints[i]
    I_tau_changepoints <- c(I_tau_changepoints, next_I0)
    next_S0 <- S_tau_changepoints[i] - n_T_segment
    S_tau_changepoints <- c(S_tau_changepoints, next_S0)
    #---

    # the number of I(tau), S(tau) at infection time only.
    I_tau_T_segment <- I_tau_segment[c(  chi_segment, FALSE)]
    S_tau_T_segment <- S_tau_segment[c(  chi_segment, FALSE)]

    if(any(I_tau_T_segment == 0))  return(list(list("I_tau_T_segment == 0", compatible = FALSE))) # depletion of infectious particles

    # the ordered waiting times within the segment
    tau_segment         <- tau_segment[order_tau_segment] # order tau_segment. could also use sort.
    dtau_segment        <- diff(c(segment_start[i], tau_segment, segment_end[i]))

    #-- Compute integrals
    integral_SI_segment <- if (gener) {
      sum(dtau_segment * I_tau_segment * S_tau_segment ^(1 - b))
    } else {
      sum(dtau_segment * I_tau_segment * S_tau_segment)
    }

    integral_I_segment  <- sum(dtau_segment * I_tau_segment)

    # Output
    #return(list(S = S_tau, I = I_tau, t = c(0, tau))) # for plotting trajectories

    #return(list(S = S_tau_segment, I = I_tau_segment, t = c(segment_start, tau_segment)))


    SS_segment <- list( # for MCMC
      compatible = TRUE,
      n_T = n_T_segment, n_J = n_J_segment,
      iota_removed = iota_removed_segment, iota_infectious = iota_infectious_segment,
      integral_SI = integral_SI_segment, integral_I = integral_I_segment,
      I_tau_T = I_tau_T_segment, S_tau_T = S_tau_T_segment
    )

    # append to the list of SS.
    SS_list <- c(SS_list, list(SS_segment))
  }
  return(SS_list)
}
