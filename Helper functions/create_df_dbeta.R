#' Create Dataframe from raw dbeta and beta lists for frequency plot of configurations.
#'
#' @param mcmc_output_dbeta delta beta vector (current)
#' @param mcmc_output_beta prior on pi.
#' @param burnin
#' @param strlastX string value of the last dbeta (change point) index.
#'
#' @return tidy dataframe for frequency plot of configurations
#' @export
#'
create_df_dbeta <- function(mcmc_output_dbeta, mcmc_output_beta, mcmc_output_loglik, burnin, strlastX) {
  # strlastX: string value of the last dbeta (change point) index.
  temp_df <- mcmc_output_dbeta[burnin:length(mcmc_output_dbeta)]
  temp_df <- data.frame(t(sapply(temp_df,c)))
  dbeta_lst_df <- temp_df %>%
    mutate(ID = burnin:length(mcmc_output_dbeta)) %>%
    unite("delta_beta", X1:strlastX, remove=FALSE)

  beta_lst_raw <- mcmc_output_beta[burnin:length(mcmc_output_beta)]
  beta_lst_df <- data.frame(t(sapply(beta_lst_raw,c)))
  beta_lst_df <- beta_lst_df %>%
    mutate(ID = burnin:length(mcmc_output_beta))

  loglik_lst_raw <- mcmc_output_loglik[burnin:length(mcmc_output_loglik)]
  loglik_lst_df <- data.frame(sapply(loglik_lst_raw,c))
  loglik_lst_df <-loglik_lst_df %>%
    mutate(ID = burnin:length(mcmc_output_loglik))

  dbeta_beta_tibble <- merge(dbeta_lst_df, beta_lst_df, by="ID")
  complete_tibble <- merge(dbeta_beta_tibble, loglik_lst_df, by='ID')

  return(complete_tibble)
}




