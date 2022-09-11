################################################################################
# R-code to simulate from a Markov and Non-Markov Stochastic epidemic models
################################################################################


# This function assumes 1 initial infective and N-1 initial susceptibles
# Per-person infection rate is beta/N

simSIR.Markov <- function(N, beta, gamma) {

  # initial number of infectives and susceptibles;
  I <- 1
  S <- N-1;

  # recording time;
  t <- 0;
  times <- c(t);

  # a vector which records the type of event (1=infection, 2=removal)
  type <- c(1);

  while (I > 0) {

    # time to next event;
    t <- t + rexp(1, (beta/N)*I*S + gamma*I);
    times <- append(times, t);

    if (runif(1) < beta*S/(beta*S + N*gamma)) {
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
  }

  # record the final size , i.e. the number of initial susceptlbles who contracted the disease sometime during the epidemic.
  #
  #
  #
  final.size <- sum(type==1) - 1
  duration <- max(times)

  # record the times of events (infections/removals) as well as the type
  res <- list("t"=times, "type"=type, "final.size"=final.size, "duration" = duration);
}

# For example, one can type:

res1 <- simSIR.Markov(N=21, 8, 4)
res1$duration



