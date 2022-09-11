#' Function to analyze the proportion of 1s in each component of delta
#'
#'
#'
#' @param res_mcmc_dbeta output delta list from DA-MCMC.
#'
#' @return vector containing information about the frequency of
#' change points (proportion of MCMC samples that indicate a change point)
#' in each position.
#' @export
#'

analyze_dbeta <- function(res_mcmc_dbeta){

  proportion_ones <- numeric(length(res_mcmc_dbeta[[1]]))

  for(s in 1:length(res_mcmc_dbeta)){
    proportion_ones <- proportion_ones + res_mcmc_dbeta[[s]]
  }

  counts <- proportion_ones
  props <- proportion_ones/length(res_mcmc_dbeta)

  return(list(counts, props))
}
