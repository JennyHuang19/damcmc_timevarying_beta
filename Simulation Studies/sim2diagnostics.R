
change_point_data2 <- sim_changepoint_sem(
  S0 = 1000, I0 = 10, t_end = 10,
  theta = list(R0 = c(1.75, 1.5, 1.25, 1.00, 0.75), gamma = 1, lambda = 1, shape = 1),
  change_day = c(2,4,6,8),
  iota_dist = "exponential",
  gener = FALSE, b = 1/2,
  E0 = 0, type = "SIR"
)

cp_data2 <- change_point_data2$x

print(change_point_data2$MLE)

plot(x = change_point_data2$t[1:950], y = change_point_data2$I[1:950], main="Epidemic Trajectory",
     xlab="Time", ylab="I(t)")

cp_Y2 <- observed_data(change_point_data2, K=10)
plot(x = cp_Y2$ts, y = c(0, cp_Y2$T_k), main="Observed Counts", xlab="Time", ylab="New Infections")


# Save an object to a file
# saveRDS(change_point_data2, file = "hard_sim2.rds")

# Restore the object
hard_sim2 <- readRDS(file = "hard_sim2.rds")

plot(x = hard_sim2$t[1:950], y = hard_sim2$I[1:950], main="Epidemic Trajectory",
     xlab="Time", ylab="I(t)")
###
pdf(file = "sim2_obs_data.pdf",   # The directory you want to save the file in
    width = 7, # The width of the plot in inches
    height = 5) # The height of the plot in inches

cp_Y2 <- observed_data(hard_sim2, K=12)
plot(x = cp_Y2$ts, y = c(0, cp_Y2$T_k), main="Observed Counts", xlab="Time", ylab="New Infections")

### new plot.
df_cp_Y2 = data.frame(cp_Y2$ts, c(0, cp_Y2$T_k))

ggplot(df_cp_Y2, aes(x=cp_Y2.ts, y=c.0..cp_Y2.T_k.)) + 
  geom_point(size=6) +
  scale_y_continuous(limits = c(0, 90), breaks= seq(0, 100, 10)) +
  scale_x_continuous(limits = c(0, 10), breaks = seq(0,10,1)) +
  labs(title = "Observed Counts", y = "New Infections", x = "Time (in Weeks)") +
  theme_bw() +
  theme(text = element_text(size=25), 
        axis.text.x = element_text(size = 21),
        axis.text.y = element_text(size = 21))
###
dev.off()

## now let's run the DA-MCMC
theta0 = list(R0 = 1, beta=.00025, lambda = 1, shape = 1, gamma = 1)

test_2 <- run_DAMCMC_complete_smallh6(
  theta0, cp_Y2, N = 100000, # you should run for twice as long to accomodate the log_lik traceplot.
  rho = 0.15, param = "bg", approx = "ldp",
  iota_dist = "exponential", gener = FALSE, b = 1/2,
  thin = 10, plt = 6000,
  par_prior = list(
    a_beta = 0.1, b_beta = 0.1,
    a_gamma = 1, b_gamma = 1,
    a_R0 = 2, b_R0 = 2,
    a_lambda = 1, b_lambda = 1,
    a_pi11 = .5, b_pi11 = .5,
    a_pi01 = 50, b_pi01 = 50
  ),
  length_delta = 9,
  x_b_ratio = 2,
  gamma=1
)

# This chunk takes in the posterior draws of (Delta, Beta) and outputs it in a tidy dataframe.
sim_df_station <- create_df_dbeta(test_2$dbeta_lst, test_2$beta_lst, test_2$pi01_lst, test_2$pi11_lst,
                                  test_2$beta_loglik, burnin=1, "X9")
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
         Beta_10 = X10, 
         pi_01 = sapply.pi01_lst_raw..c.,
         pi_11 = sapply.pi11_lst_raw..c.,
         theta_log_lik = sapply.loglik_lst_raw..c.)


posterior_delta_beta <- sim_df_station %>%
  dplyr::select(Delta_1:theta_log_lik)

### PLOT OF CHANGE POINT LOCATIONS
# configurations of delta.
barplot(analyze_dbeta(test_2$dbeta_lst)[[2]], main="Posterior Change Point Locations", xlab="Change Point Location",
        ylab="Proportion of MCMC Samples",  border="blue", names.arg=c("1","2","3","4","5","6","7","8","9"))
abline(h=1/2, col="blue", lty="dashed", lwd=4)

df_posterior_cp = data.frame(Week = seq(1,9,1),
                             cp_counts = analyze_dbeta(test_2$dbeta_lst)[[2]])

ggplot(df_posterior_cp, aes(Week, cp_counts)) +
  geom_col(fill="black") +
  scale_y_continuous(limits = c(0, 0.5), breaks= seq(0, 0.5, 0.05)) +
  scale_x_continuous(limits = c(0, 10), breaks = seq(1,9,1)) +
  labs(title = "Posterior Change Point Locations", 
       y = "Posterior Probability", x = "Time (in Weeks)") +
  theme_bw() +
  theme(text = element_text(size=23), 
        axis.text.x = element_text(size = 19),
        axis.text.y = element_text(size = 19))
###

## convert the output into coda format
coda_chain_simstudy = mcmc(posterior_delta_beta)

coda_chain_simstudy[1,] # look at the first row, all columns.

summary(coda_chain_simstudy) # the summary of each variable of interest.

## we can also compute Monte Carlo error for the quantiles of each variable in the coda_gibbs_chain mat # ola.
mcse.q.mat(coda_chain_simstudy, 0.025)
mcse.q.mat(coda_chain_simstudy, 0.975)

effectiveSize(coda_chain_simstudy)

## plot traceplots and posterior densities
plot(coda_chain_simstudy)

## look at what coda can do
help(package=coda)

## look at the menu options
codamenu()

### Spaghetti Plot of Posterior R0
sim2_delta_beta <- sim_df_station %>% # get rid of delta_beta column. we won't need it from here on.
  mutate(Beta_11 = Beta_10) %>%
  dplyr::select(ID, Beta_1:Beta_11,
                -pi_01, -pi_11,
                -theta_log_lik)

sim2_time_series <- sim2_delta_beta  %>%
  pivot_longer(cols = -ID,
               names_to = "time_unit", names_prefix = "V",
               # names_transform = list(Time = as.integer),
               values_to = "Beta") %>%
  filter(ID %% 10 == 0)

#preserve the order of time units.
sim2_time_series$time_unit <- factor(sim2_time_series$time_unit, levels=unique(sim2_time_series$time_unit))

#mutate(r0)
sim2_time_series <- sim2_time_series %>%
  mutate(R0 = Beta * cp_Y2$S0 / 1)

# Plot r0
# spagetti plot of R0 over time
sim2_time_series %>%
  ggplot(aes(x = time_unit, y = Beta*1000, group = ID)) +
  geom_step(alpha = 0.1, color = "red") +
  labs(y = "Transmission Rate", x = "Time (in Weeks)") +
  scale_x_discrete(labels=seq(0,12,1)) +
  scale_y_continuous(limits = c(0, 3.5)) +
  expand_limits(y = 0) +
  theme_bw() +
  theme(text = element_text(size=25), 
        axis.text.x = element_text(size = 21),
        axis.text.y = element_text(size = 21))
###



