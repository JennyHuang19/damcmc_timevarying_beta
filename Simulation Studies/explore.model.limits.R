# 5 x 5 (size of change points (0.75, 0.5, 0.25, 0.1) x frequency of change points.)
create_epidemic <- function(S0 = 10000, I0 = 25, t_end = 12,
                            theta = list(R0 = c(1.75, 1.25, 0.75), gamma = 1, lambda = 1, shape = 1),
                            change_day = c(3,10),
                            iota_dist = "exponential",
                            gener = FALSE, b = 1/2,
                            E0 = 0, type = "SIR"){

  traj <- sim_changepoint_sem(
    S0, I0, t_end,
    theta,
    change_day,
    iota_dist,
    gener, b,
    E0, type = "SIR"
  )

  # Below, we visualize the I(t) compartment over time.
  plot(x = traj$t[1:6000], y = traj$I[1:6000], main="Epidemic Trajectory",
       xlab="Time", ylab="I(t)")

  obs.Y <- observed_data(traj, K=10)
  plot(x = obs.Y$ts, y = c(0, obs.Y$T_k),
       main="Observed Counts", xlab="Time (in Weeks)", ylab="New Infections"
       , cex = 2.5)

  return(list(traj, obs.Y))
}


# the difficult part is generating epidemics that
# 1. are interesting. (when R0 crosses 0)
# 2. do not die out immediately.
limits.0.5 <- create_epidemic(
  S0 = 10000, I0 = 25, t_end = 12,
  theta = list(R0 = c(1.75, 1.25, 0.75), gamma = 1, lambda = 1, shape = 1),
  change_day = c(3,10),
  iota_dist = "exponential",
  gener = FALSE, b = 1/2,
  E0 = 0, type = "SIR")

limits.0.25 <- create_epidemic(
  S0 = 10000, I0 = 25, t_end = 12,
  theta = list(R0 = c(1.25, 1.00, 0.75), gamma = 1, lambda = 1, shape = 1),
  change_day = c(3,10),
  iota_dist = "exponential",
  gener = FALSE, b = 1/2,
  E0 = 0, type = "SIR"
)

limits.0.1 <- create_epidemic(
  S0 = 10000, I0 = 25, t_end = 12,
  theta = list(R0 = c(1.5, 1.1, 1.0), gamma = 1, lambda = 1, shape = 1),
  change_day = c(3,10),
  iota_dist = "exponential",
  gener = FALSE, b = 1/2,
  E0 = 0, type = "SIR"
)

### Frequency of change points (placed on a regular schedule).
limits.f3 <- create_epidemic(
  S0 = 10000, I0 = 25, t_end = 12,
  theta = list(R0 = c(1.5, 1.25, 1.00, 0.75), gamma = 1, lambda = 1, shape = 1),
  change_day = c(2,5,8),
  iota_dist = "exponential",
  gener = FALSE, b = 1/2,
  E0 = 0, type = "SIR"
)

limits.f2 <- create_epidemic(
  S0 = 10000, I0 = 25, t_end = 12,
  theta = list(R0 = c(1.5, 1.25, 1.00, 0.75), gamma = 1, lambda = 1, shape = 1),
  change_day = c(2,4,6),
  iota_dist = "exponential",
  gener = FALSE, b = 1/2,
  E0 = 0, type = "SIR"
)

limits.f1 <- create_epidemic(
  S0 = 10000, I0 = 25, t_end = 12,
  theta = list(R0 = c(1.5, 1.25, 1.00, 0.75), gamma = 1, lambda = 1, shape = 1),
  change_day = c(2,3,4),
  iota_dist = "exponential",
  gener = FALSE, b = 1/2,
  E0 = 0, type = "SIR"
)

# Save the data object to a file
# saveRDS(clearcut_sim1, file = "clearcut_sim1.rds")
# Restore the object
# change_point_data <- readRDS(file = "clearcut_sim1.rds")

# Run MCMC
theta0 = list(R0 = 1, beta=1.927e-4, lambda = 1, shape = 1, gamma = 1)
## now let's run the DA-MCMC
test.limits.0.1 <- run_DAMCMC_complete_smallh6(
  theta0, limits.0.1[[2]], N = 50000, # you should run for twice as long to accomodate the log_lik traceplot.
  rho = 0.05, param = "bg", approx = "ldp",
  iota_dist = "exponential", gener = FALSE, b = 1/2,
  thin = 1, plt = 6000,
  par_prior = list(
    a_beta = 0.1, b_beta = 0.1,
    a_gamma = 1, b_gamma = 1,
    a_R0 = 2, b_R0 = 2,
    a_lambda = 1, b_lambda = 1,
    a_pi11 = .5, b_pi11 = .5,
    a_pi01 = .5, b_pi01 = .5
  ),
  length_delta = 11,
  x_b_ratio = 2,
  gamma=1
)
###

# Post-Processing
# Outputs the posterior draws of (Delta, Beta) and in a tidy dataframe.
sim_df_stationary <- create_df_dbeta(test.limits.0.1$dbeta_lst, test.limits.0.1$beta_lst, test.limits.0.1$pi01_lst, test.limits.0.1$pi11_lst,
                                     test.limits.0.1$beta_loglik, burnin=1, "X11")

sim_df_stationary <- sim_df_stationary %>%
  rename(Beta_1 = X1.y, Delta_1 = X1.x,
         Beta_2 = X2.y, Delta_2 = X2.x,
         Beta_3 = X3.y, Delta_3 = X3.x,
         Beta_4 = X4.y, Delta_4 = X4.x,
         Beta_5 = X5.y, Delta_5 = X5.x,
         Beta_6 = X6.y, Delta_6 = X6.x,
         Beta_7 = X7.y, Delta_7 = X7.x,
         Beta_8 = X8.y, Delta_8 = X8.x,
         Beta_9 = X9.y, Delta_9 = X9.x,
         Beta_10 = X10.y, Delta_10 = X10.x,
         Beta_11 = X11.y, Delta_11 = X11.x,
         Beta_12 = X12,
         pi_01 = sapply.pi01_lst_raw..c.,
         pi_11 = sapply.pi11_lst_raw..c.,
         theta_log_lik = sapply.loglik_lst_raw..c.)


posterior_delta_beta <- sim_df_stationary %>%
  dplyr::select(Delta_1:Beta_12)

## convert the output into coda format
coda_chain_simstudy = mcmc(posterior_delta_beta)
# the summary of each variable of interest.
summary(coda_chain_simstudy)
effectiveSize(coda_chain_simstudy)
## look at the menu options
codamenu()

### Plot the change point locations.
# configurations of delta.
df_posterior_cp = data.frame(Week = seq(1,11,1),
                             cp_counts = analyze_dbeta(test.limits.0.1$dbeta_lst)[[2]])

ggplot(df_posterior_cp, aes(Week, cp_counts)) +
  geom_col(fill="black") +
  scale_y_continuous(limits = c(0, 1), breaks= seq(0, 1, 0.2)) +
  scale_x_continuous(limits = c(0, 12), breaks = seq(1,11,1)) +
  labs(title = "Posterior Change Point Locations",
       y = "Posterior Probability", x = "Time (in Weeks)") +
  theme_bw() +
  theme(text = element_text(size=23),
        axis.text.x = element_text(size = 19),
        axis.text.y = element_text(size = 19))
###

### Plot the Posterior Transmission Rates
sim1_delta_beta <- sim_df_stationary %>% # get rid of delta_beta column. we won't need it from here on.
  mutate(Beta_13 = Beta_12) %>%
  dplyr::select(ID, Beta_1:Beta_13, -theta_log_lik, -pi_01, -pi_11)

sim1_time_series <- sim1_delta_beta  %>%
  pivot_longer(cols = -ID,
               names_to = "time_unit", names_prefix = "V",
               # names_transform = list(Time = as.integer),
               values_to = "Beta") %>%
  filter(ID %% 10 == 0)

#preserve the order of time units.
sim1_time_series$time_unit <-
  factor(sim1_time_series$time_unit, levels=unique(sim1_time_series$time_unit))

# spagetti plot over time
sim1_time_series %>%
  ggplot(aes(x = time_unit, y = Beta*10000, group = ID)) +
  geom_step(alpha = 0.1, color = "red") +
  labs( y = "Transmission Rate", x = "Time (in Weeks)") +
  scale_x_discrete(labels=seq(0,12,1)) +
  expand_limits(y = 0) +
  theme_bw() +
  theme(text = element_text(size=25),
        axis.text.x = element_text(size = 21),
        axis.text.y = element_text(size = 21))
###

