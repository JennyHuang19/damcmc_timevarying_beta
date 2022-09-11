create_mcmc_df <- function(MC, thin, burnin) {

  # create data frame
  theta         <- MC[["theta"      ]] # a list where each element is 1 iteration of values.
  loglik        <- MC[["loglik"     ]]
  thin          <- 1
  # if(do_SS)  SS <- MC[["SS"         ]]
  rate_accept   <- MC[["rate_accept"]]
  S0            <- MC[["S0"         ]]
  run_time      <- MC[["run_time"   ]]

  # Data Wrangling
  theta_tidy <- data.table::rbindlist(theta) %>%
    add_iteration(thin) %>%
    dplyr::mutate(loglik = loglik) %>%
    dplyr::filter(.data$Iteration <= 20000) %>%
    remove_burnin(burnin) %>%
    dplyr::mutate(expected_infection_length = 1 / gamma)

  return(theta_tidy)

}
