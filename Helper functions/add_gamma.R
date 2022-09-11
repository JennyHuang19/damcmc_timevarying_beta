#' Include gamma to the parameter vector theta
#'
#' @param theta parameters (without beta)
#' @param S0 initial population size
#'
#' @return theta including beta
#' @export
#'
add_gamma <- function(theta) {

  lambda <- theta[["lambda"]]
  shape  <- theta[["shape" ]]

  theta[["gamma"]] <- 1 / (lambda ^ (- 1 / shape) * gamma(1 + 1 / shape))

  return(theta)

}
