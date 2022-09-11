#' Include beta to the parameter vector theta
#'
#' @param theta parameters (without beta)
#' @param S0 initial population size
#'
#' @return theta including beta
#' @export
#'

add_beta <- function(theta, S0) {

  gamma <- theta[["gamma"]]
  R0    <- theta[["R0"   ]]

  theta[["beta"]] <- R0 * gamma / S0 #if R0 is a vector, beta is a vector.

  return(theta)

}
