
#' create tidy dataframe from beta_list_untidy
#'
#' @param beta_list_untidy beta samples
#'
#' @return tidy dataframe of beta samples.
#' @export
#'
create_mcmc_list <- function(beta_list_untidy){
  # make each beta into a plain list
  all_beta <- list()

  for(k in 1:length(beta_list_untidy[[1]])){
    index <- paste("beta", k, sep="_")
    beta_i = unlist(lapply(beta_list_untidy, `[[`, k)) # beta k
    all_beta <- append(all_beta, list(beta_i))
  }

  return(all_beta)
}
