#' Generates a tidy data frame for plotting a parameter (pi or beta) at a time instant of interest
#'
#' @param pi_mcmc_list output of create_mcmc_list (row x col = timeseries_length x iteration)
#' @param time_instant parameter at the time instant of interest.
#'
#' @return a tidy dataframe of the parameter at the time instance specified.
#' @export
#'
param_tidy <- function(pi_mcmc_list, time_instant){
  parameter_df <- do.call(rbind, Map(data.frame, Iteration = list(c(1:length(pi_mcmc_list[[1]]))), pi = pi_mcmc_list[time_instant]))
  return(parameter_df)
}
