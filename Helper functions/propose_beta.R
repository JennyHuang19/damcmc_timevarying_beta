# Update (pi | delta)


#' Gibbs for component-wise model of pi
#'
#' @param delta_count vector of delta counts, where component i is the number of 1's in position i aggregated over the past X iterations.
#' @param num_iter_pi number of iterations between updates for pi
#' @param a_prior vector of 'a' parameter for pi.
#' @param b_prior vector of 'b' parameter for pi.
#'
#' @return new draw of delta
#' @export
#'
gibbs_pi6 <- function(a_prior, b_prior, delta_count, num_iter_pi){
  return(rbeta( length(delta_count), a_prior + delta_count, b_prior + (rep(num_iter_pi, length(delta_count)) - delta_count) ))
}

#' Gibbs update for pi conditional on delta
#'
#' @param delta delta beta vector (current)
#' @param par_prior prior on pi.
#'
#' @return new draw of delta
#' @export
#'
gibbs_pi <- function(par_prior, delta) {

  # retrieve statistics n00, n01, n10, n11 from delta
  n00 <- n01 <- n10 <- n11 <- 0

  for(t in 1:(length(delta)-1)){
    if(delta[t] == 0 & delta[t+1] == 0){
      n00 <- n00 + 1
    }
    if(delta[t] == 0 & delta[t+1] == 1){
      n01 <- n01 + 1
    }
    if(delta[t] == 1 & delta[t+1] == 0){
      n10 <- n10 + 1
    }
    if(delta[t] == 1 & delta[t+1] == 1){
      n11 <- n11 + 1
    }
  }

  # Gibbs update.
  pi01 <- rbeta(1, par_prior['a_pi01'][[1]] + n01, par_prior['b_pi01'][[1]] + n00)
  pi11 <- rbeta(1, par_prior['a_pi11'][[1]] + n11, par_prior['b_pi11'][[1]] + n10)

  return(c(pi01, pi11))
}


# draw from (delta | pi)

#' Propose a new delta through the Markov chain model.
#'
#' @param pi01, pi11 transition probabilities
#' @param length_delta length of delta vector
#'
#' @return new draw of delta
#' @export
#'
draw_mm_delta <- function(pi01, pi11, length_delta) {

  # initialize vector of len(length_delta)
  delta <- rep(0,length_delta)

  # draw initial state
  delta[0] <- Rlab::rbern(1, 0.1) # h0 = 0.1

  for(t in 1:(length(delta)-1)){
    if(delta[t] == 0){
      delta[t+1] <- Rlab::rbern(1, pi01)
    }
    if(delta[t] == 1){
      delta[t+1] <- Rlab::rbern(1, pi11)
    }
  }
  return(delta)
}



# small proposal delta | pi

#' Propose a new delta through the Markov chain model only updating one component of the current delta vector.
#'
#' @param pi01, pi11 transition probabilities
#' @param delta current delta vector
#'
#' @return new draw of delta
#' @export
#'
draw_sm_delta <- function(pi01, pi11, delta) {

  # choose randomly 1 component in the delta vector to change (u+1).
  u <- sample(0:(length(delta)-1), 1)

  if(u == 0){
    delta[u+1] <- Rlab::rbern(1, pi01)
  }
  else if(delta[u] == 0){
    delta[u+1] <- Rlab::rbern(1, pi01)
  }
  else{
    delta[u+1] <- Rlab::rbern(1, pi11)
  }

  return(list(delta,u))
}

# delta | pi (for the component-wise model)

#' Propose a new delta through the Markov chain model only updating one component of the current delta vector.
#'
#' @param pi, vector of cohesions.
#' @param delta current delta vector
#'
#' @return new draw of delta
#' @export
#'
draw_delta_bern <- function(pi, delta) {

  # choose randomly 1 component in the delta vector to change (u+1).
  u <- sample(1:(length(delta)-1), 1)

  delta[u] <- Rlab::rbern(1, pi[u])

  return(list(delta,u))
}






# draw from (dbeta_star | dbeta)

#' Propose a new delta beta through flipping 1 component of the current dbeta vector at random.
#'
#' @param dbeta delta beta vector (current)
#' @param is_prob boolean telling whether to use the method of changing dbeta depending on a bernoulli draw with probability h_hat.
#' @param h_hat prior inclusion probabilities (i.e. the P(dbeta_i = 1))
#'
#' @return new draw of dbeta, u - the index or time instance of dbeta which was changed.
#' @export
#'
draw_dbeta <- function(dbeta, is_prob, h_hat) {

  # choose randomly 1 component in the dbeta vector to change.
  u <- sample(1:length(dbeta), 1) # sample.int(6,1)

  if(is_prob){ # h_hat method
    changed_dbeta <- Rlab::rbern(1, h_hat[u]) # propose delta, beta jointly every fraction of iterations. rather than proposing from the prior. also consider proposing latent data 5x the proposal of theta.
    dbeta[u] <- changed_dbeta
  }

  else{ # uniform flip method
    if(dbeta[u] == 1){
      dbeta[u] <- 0
    }
    else{
      dbeta[u] <- 1
    }
  }


  return(list(dbeta, u))
}



# draw from (beta | dbeta_star, Z)

#' Propose a new beta vector using gibbs_theta conditional on the location of change points (dbeta) and data.
#'
#' @inheritParams b_gibbs_theta Gibb's Sampler for beta_vector.
#' @param dbeta_star delta beta vector (locations of change points)
#' @param segmented_SS sufficient statistics of the current configuration of the latent data

#' @return new draw of beta (where a component beta_i indicates a beta value for a time unit, such as a week)
#' @export
#'

draw_beta <- function(dbeta_star, segmented_SS, iota_dist, gamma, par_prior, param, Y){

  # Function to draw a beta vector given the sufficient statistics and locations of the change-points.

  # 1. initialize beta_star
  beta_star <- c()

  # 2. propose beta for the first segment.
  segment_number = 1
  beta_i <- b_gibbs_theta(segmented_SS[[segment_number]],  iota_dist="exponential", par_prior, gamma, param, Y)
  beta_star <- append(beta_star, beta_i)

  # 3. propose beta for each following segment.
  for (i in 1:length(dbeta_star)) {
    if(dbeta_star[i] == 1){ # propose new beta_i each time the dbeta vector indicates a change point.
      segment_number = segment_number + 1
      beta_i <- b_gibbs_theta(segmented_SS[[segment_number]], iota_dist="exponential", par_prior, gamma, param, Y)
      beta_star <- append(beta_star, beta_i)
    }
    else{ # else, append the previous week's beta_i.
      beta_star <- append(beta_star, beta_i)
    }
  }
  return(beta_star)
}
