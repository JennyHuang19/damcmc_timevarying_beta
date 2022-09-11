#' Function to format data for creating trace plots and posterior density estimates using the coda library.
#'
#'
#'
#' @param acc_beta_lst output beta list from run_damcmc_beta.
#'
#' @return mcmc.trace
#' @export
#'

create_tp_data <- function(acc_beta_lst){

  # Change output into a df
  library('coda')
  matrix <- do.call(rbind, acc_beta_lst)
  df <- data.frame(matrix)
  mcmc.trace <- mcmc(df)

  return(mcmc.trace)
}
