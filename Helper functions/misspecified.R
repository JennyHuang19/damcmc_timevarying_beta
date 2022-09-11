# Run sim_study_1 with cp location fixed but fixed incorrectly. 3 --> 2, 10 --> 9.
test_6.2 <- run_DAMCMC_complete_smallh6.2(
  theta0, cp_Y, N = 50000,
  rho = 0.2, param = "bg", approx = "ldp",
  iota_dist = "exponential", gener = FALSE, b = 1/2,
  thin = 10, plt = 6000,
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
  gamma=1,
  fixed_delta = c(1,1,0,1,0,0,0,0,1,0,0)
)
# create marginal plot
print(analyze_dbeta(test_6.2$dbeta_lst))

barplot(analyze_dbeta(test_6.2$dbeta_lst)[[2]], main="Posterior Change Point Locations", xlab="Change Point Location",
        ylab="Proportion of MCMC Samples",  border="blue", names.arg=c("1","2","3","4","5","6","7","8","9","10","11"),ylim=c(0,1))

# create spaghetti plot
# This chunk takes in the posterior draws of (Delta, Beta) and outputs it in a tidy dataframe.
sim_df_station <- create_df_dbeta(test_6.2$dbeta_lst, test_6.2$beta_lst, test_6.2$beta_loglik, burnin=1, "X11")
sim_df_station <- sim_df_station %>%
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
         theta_log_lik = sapply.loglik_lst_raw..c.)

sim1_delta_beta <- sim_df_station %>% # get rid of delta_beta column. we won't need it from here on.
  mutate(Beta_13 = Beta_12) %>%
  select(ID, Beta_1:Beta_13, -theta_log_lik)

sim1_time_series <- sim1_delta_beta  %>%
  pivot_longer(cols = -ID,
               names_to = "time_unit", names_prefix = "V",
               # names_transform = list(Time = as.integer),
               values_to = "Beta") %>%
  filter(ID %% 10 == 0)

#preserve the order of time units.
sim1_time_series$time_unit <- factor(sim1_time_series$time_unit, levels=unique(sim1_time_series$time_unit))

#mutate(r0)
sim1_time_series <- sim1_time_series %>%
  mutate(R0 = Beta * cp_Y$S0 / 1)
# Plot r0
# spaghetti plot of beta over time.
sim1_time_series %>%
  ggplot(aes(x = time_unit, y = Beta, group = ID)) +
  geom_step(alpha = 0.1, color = "red") +
  labs(title = "Posterior Estimate of Transmission Rate", y = "B(t)", x = "Time") +
  scale_x_discrete(labels=seq(0,12,1)) +
  scale_y_continuous(limits=c(0,2.5e-4)) +
  # expand_limits(y = c(0, 2.5e-4)) +
  theme_bw()
###

### Posterior Summary Stats.
posterior_delta_beta <- sim_df_station %>%
  select(Delta_1:theta_log_lik)
## convert the output into coda format
coda_chain_simstudy = mcmc(posterior_delta_beta)

coda_chain_simstudy[1,] # look at the first row, all columns.

summary(coda_chain_simstudy) # the summary of each variable of interest.

# Change point on time unit 3 is missed and change point on time unit 9 occurs a unit earlier.
# misspecified 1.
# 1.213e-04 [1.143e-04, 1.285e-04]

# misspecified 2.
# before beta_9: 1.314e-04 [1.232e-04 1.402e-04]
# after beta_9: 9.212e-05 [8.784e-05 1.046e-04]

