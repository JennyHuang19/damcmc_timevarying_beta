################################################################################
# R-code to simulate from a Markov and Non-Markov Stochastic epidemic models
################################################################################


# This function assumes 1 initial infective and N-1 initial susceptibles
# Per-person infection rate is beta/N
# beta is a time-dependent vector.

simSIR.Markov <- function(N, I0, t_end, beta, gamma) {

  # initial number of infectives and susceptibles;
  I <- I0
  S <- N-1;

  # recording time;
  t <- 0;
  times <- c(t);
  cumu_cases <- c(I0);
  T_k_times <- c(0);
  cumu_cases_T_k <- c(0);

  # a vector which records the type of event (1=infection, 2=removal)
  type <- c(1);

  while (I > 0 & t <= t_end) {

    ### fetch the current beta value.
    # take the floor of the current time to be the index_for_beta.
    if (t >= 1 & t < t_end) {
      index_for_beta = floor(t)
      beta_curr = beta[index_for_beta+1] # 1.5 -> 1 -> 2nd beta.
    } else if (t < 1) {
      beta_curr <- beta[1] # at the beginning of the epidemic, set beta to be beta[1].
    } else
      beta_curr <- beta[t_end] # beyond the observed data portion of the epidemic, set beta to be beta[t_end].

    # time to next event;
    prev_t  = t
    t <- t + rexp(1, (beta_curr/N)*I*S + gamma*I);
    times <- append(times, t);

    if (runif(1) < beta_curr*S/(beta_curr*S + N*gamma)) {
      # infection
      I <- I+1;
      S <- S-1;
      type <- append(type, 1);
    }
    else {
      #removal
      I <- I-1
      type <- append(type, 2);
    }
    # append current cumulative cases.
    cumu_cases <- append(cumu_cases, sum(type==1));

    ### For new_cases. If t crosses an integer, append the cumu case count.
    if(floor(prev_t) < floor(t)){
      T_k_times <- append(T_k_times, t);
      cumu_cases_T_k <- append(cumu_cases_T_k, sum(type==1));
    }
  }

  # record the final size , i.e. the number of initial susceptlbles who contracted the disease sometime during the epidemic.
  #
  #
  #
  final.size <- sum(type==1) - I0
  duration <- max(times)

  # record the times of events (infections/removals) as well as the type
  res <- list("t"=times, "type"=type, "cumu_cases"=cumu_cases, "final.size"=final.size, "duration" = duration,
              "T_k_times"=T_k_times, "cumu_cases_T_k"=cumu_cases_T_k);
}

# For example, one can type:
res1 <- simSIR.Markov(N=5000, I0=1, t_end=12, beta=rep(3,12), gamma = 1)

### Simulation 1
# new cases
res1$cumu_cases_T_k
T_k <- diff(res1$cumu_cases_T_k) # the difference in cumulative cases
Observation_times <- res1$T_k_times # times at which we cross an integer for the first time.
plot(Observation_times, T_k, col="blue")
# times, types, final size
res1$t
res1$type
res1$final.size
# cumulative cases
plot(res1$t, res1$cumu_cases, type="l")
###

