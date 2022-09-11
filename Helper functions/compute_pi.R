#' Probability of being removed before next observation time
#'
#' Computes p_i, the probability of being removed before t_end given an infection at time tau_T.
#'
#' @inheritParams rprop_x
#'
#' @param tau_T observed event times
#' @param t_end time before which individuals may be removed
#'
#' @return a vector of the probabilities that individual infected at time tau_T are removed before t_end
#' @export
#'
compute_pi <- function(theta, tau_T, iota_dist, t_end) {

  # Setup
  gamma  <- theta[["gamma" ]]
  lambda <- theta[["lambda"]]
  shape  <- theta[["shape" ]]

  # Compute p_i
  p_i <- if(iota_dist == "exponential") {
    stats::pexp(t_end - tau_T, gamma)
  } else if(iota_dist == "weibull") {
    pweibull2(t_end - tau_T, shape, lambda)
  }

  # Output
  return(p_i)

}
