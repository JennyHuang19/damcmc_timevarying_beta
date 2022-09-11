#' Return a df where each row is a list of change point positions for one MCMC sample.
#'
#' @param mcmc_output_dbeta delta beta vector (current)
#' @param burnin
#'
#' @return a list where each element is a list of change point positions for a sample.
#' @export
#'
changepoint_locations <- function(mcmc_output_dbeta, burnin) {

  dbeta <- mcmc_output_dbeta[burnin:length(mcmc_output_dbeta)]

  N <- length(dbeta)
  locations <- vector(mode = "list", length = N)

  for(j in 1:N){
    locations[j] <- list(which(dbeta[[j]] == 1))
  }

  locations_df <- plyr::ldply(locations, rbind)

  return(locations_df)
}


