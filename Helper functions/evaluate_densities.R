# P(dbeta)

#' LOG PRIOR DENSITY OF DELTA-BETA.
#'
#' @param dbeta delta beta vector (current)
#' @param h_hat prior for P(dbeta_i = 1)
#'
#' @return the prior density evaluated at dbeta
#' @export
#'
prior_dbeta <- function(h_hat, dbeta){ # dbern(x, prob, log = FALSE)
  total_prob <- 0
  for(i in 1:length(dbeta)){
    if(dbeta[i] == 1){
      total_prob <- total_prob + log(h_hat[i]) # h_i = p(dbeta_i = 1)
    }
    else{
      total_prob <- total_prob + log(1-h_hat[i])
    }
  }
  return(total_prob)
}



#' LOG PRIOR DENSITY OF BETA (SUM GAMMA DENSITIES)
#'
#' @param beta_star beta vector (current)
#' @param par_prior parameters of the prior distributions
#'
#' @return the prior density evaluated at beta_star
#' @export
#'
prior_beta <- function(beta_star, par_prior){
  # EVALUATE PRIOR DENSITY OF BETA (PRODUCT OF GAMMA DENSITIES)
  # beta_star is the vector of betas returned by draw_beta

  # evaluate f_beta on the first segment.
  f_beta <- stats::dgamma(beta_star[1], par_prior['a_beta'][[1]], par_prior['b_beta'][[1]], log = TRUE)

  # evaluate f_beta on each of the subsequent segments.
  for(i in 2:length(beta_star)){
    if(beta_star[i] != beta_star[i-1]){
      f_beta <- f_beta + stats::dgamma(beta_star[i], par_prior['a_beta'][[1]], par_prior['b_beta'][[1]], log = TRUE)
    }

  }
  return(f_beta)
}


#' LOG PRIOR DENSITY OF PI (SUM BERNOULLI DENSITIES)
#'
#' @param Pi_star pi vector (current)
#' @param par_prior parameters of the prior distributions
#'
#' @return the prior density evaluated at beta_star
#' @export
#'
prior_pi <- function(Pi_star, par_prior){
  return( sum( dbeta(Pi_star, par_prior['a_pi'][[1]], par_prior['b_pi'][[1]], log = TRUE) ) ) # sum of vector components.
}





### B. Evaluating proposal densities.

#' LOG PROPOSAL DENSITY OF PI (GIBBS PROPOSAL, Be(a_star, b_star))
#'
#' @param Pi_star Pr(Delta_i = 1) vector
#' @param delta delta vector (current)
#' @param par_prior parameters of the prior distributions # note that if the prior > length(delta), then the prior is stronger than the data.
#'
#' @return the proposal density evaluated at Pi
#' @export
#'
f_pi_func <- function(Pi_star, delta, par_prior) {
  return( sum( dbeta(Pi_star, par_prior['a_pi'][[1]] + sum(delta), par_prior['b_pi'][[1]] + (length(delta) - sum(delta)), log = TRUE) ) )
}



#' LOG PROPOSAL DENSITY OF DELTA-BETA (DRAW 1 COMPONENT FROM A UNIFORM DISTRIBUTION THEN DRAWN FROM A BERNOULLI(H_HAT))
#'
#' @param dbeta delta beta vector (current)
#' @param Pi the prior inclusion probabilities (i.e. the P(dbeta_i = 1))
#' @param u index of dbeta that was changed in the current proposal of dbeta_star.
#'
#' @return the proposal density evaluated at dbeta
#' @export
#'
f_dbeta_func <- function(dbeta_j, h_hat, u) {
  return(Rlab::dbern(dbeta_j[u], h_hat[u], log = TRUE))
}




#' LOG PROPOSAL DENSITY OF BETA (SUM GAMMA DENSITIES)
#'
#' @param beta_star beta_star is the beta vector returned by draw_beta.
#' @param par_prior parameters of the priors
#' @param segmented_SS sufficient statistics of the current configuration of latent data
#'
#' @return the proposal density evaluated at beta_star.
#' @export
#'
f_beta_func <- function(beta_star, par_prior, segmented_SS){

  # evaluate f_beta on the first segment.
  f_beta <- stats::dgamma(beta_star[1], par_prior['a_beta'][[1]] + segmented_SS[[1]]$n_T, par_prior['b_beta'][[1]] + segmented_SS[[1]]$integral_SI, log = TRUE)

  # case when there is only 1 segment
  if(length(beta_star) == 1) return(f_beta)

  # evaluate f_beta on each of the subsequent segments.
  for(i in 2:length(beta_star)){
    if(beta_star[i] != beta_star[i-1]){ # an occurrence of a change in the beta vector.
      f_beta <- f_beta + stats::dgamma(beta_star[i], par_prior['a_beta'][[1]] + segmented_SS[[i]]$n_T, par_prior['b_beta'][[1]] + segmented_SS[[i]]$integral_SI, log = TRUE)
    }
  }
  return(f_beta)
}




#' LOG PROPOSAL DENSITY OF DELTA (Markov Model with Transition Probabilities Pi)
#'
#' @param delta delta vector (proposed)
#' @param pi01, pi11 Pi: the transition probabilities
#'
#' @return the proposal density evaluated at delta
#' @export
#'
f_delta_mm <- function(delta, pi01, pi11) {

  # initialize log-likelihood.
  ans <- 0

  for(k in 1:(length(delta) - 1)){
    if(delta[k] == 0 & delta[k+1] == 0){
      ans <- ans + Rlab::dbern(0, pi01, log = TRUE)
    }
    if(delta[k] == 0 & delta[k+1] == 1){
      ans <- ans + Rlab::dbern(1, pi01, log = TRUE)
    }
    if(delta[k] == 1 & delta[k+1] == 0){
      ans <- ans + Rlab::dbern(0, pi11, log = TRUE)
    }
    if(delta[k] == 1 & delta[k+1] == 1){
      ans <- ans + Rlab::dbern(1, pi11, log = TRUE)
    }
  }
  return(ans)
}



#' LOG PROPOSAL DENSITY OF DELTA (Markov Model with Small Proposal)
#'
#' @param delta delta vector (proposed)
#' @param pi01, pi11 Pi: the transition probabilities
#'
#' @return the proposal density evaluated at delta
#' @export
#'
f_delta_sm <- function(delta, pi01, pi11, u) {
  if(u == 0){
    return(Rlab::dbern(delta[u+1], pi01, log = TRUE))
  }

  else if(delta[u] == 0){
    return(Rlab::dbern(delta[u+1], pi01, log = TRUE))
  }
  else{
    return(Rlab::dbern(delta[u+1], pi11, log = TRUE))
  }
}
