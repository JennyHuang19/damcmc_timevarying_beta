## This script runs mcmc diagnostics for simulation study 1.

## first we need to load coda and mcmcse packages (you need to install them first)
library(coda)
library(mcmcse)

### Restore the epidemic object
clearcut_sim1 <- readRDS(file = "clearcut_sim1.rds")

plot(x = clearcut_sim1$t[1:5000], y = clearcut_sim1$I[1:5000], main="Epidemic Trajectory",
     xlab="Time", ylab="I(t)")
cp_Y <- observed_data(clearcut_sim1, K=11)
plot(x = cp_Y$ts, y = c(0, cp_Y$T_k), 
     main="Observed Counts", xlab="Time (in Weeks)", ylab="New Infections"
     , cex = 2.5)

df_cp_Y = data.frame(cp_Y$ts, c(0, cp_Y$T_k))

ggplot(df_cp_Y, aes(x=cp_Y.ts, y=c.0..cp_Y.T_k.)) + 
  geom_point(size=6) +
  scale_y_continuous(limits = c(0, 350), breaks= seq(0, 350, 50)) +
  scale_x_continuous(limits = c(0, 12), breaks = seq(0,12,1)) +
  labs(title = "Observed Counts", y = "New Infections", x = "Time (in Weeks)") +
  theme_bw() +
  theme(text = element_text(size=25), 
        axis.text.x = element_text(size = 21),
        axis.text.y = element_text(size = 21))

###

### Now let's decide on the hyperparameters.
# pi's

#define range
bvals = seq(0,1, length=100)
#plot Gamma distributions
plot(p, dgamma(bvals, 1, 1), ylab='Density', xlab='Value for Beta', type ='l', col='black')
#lines(p, dbeta(p, 0.5, 0.5), col='red') 
#lines(p, dbeta(p, 5, 2), col='blue')

#define range
p = seq(0,1, length=100)
#plot Beta distributions
plot(p, dbeta(p, 0.5, 0.5), ylab='Density', xlab='Value for Pi', type ='l', col='black')
#lines(p, dbeta(p, 0.5, 0.5), col='red') 
#lines(p, dbeta(p, 5, 2), col='blue')

#add legend
legend(.7, 4, c('Beta(2, 10)','Beta(2, 2)','Beta(1,1)'),
       lty=c(1,1,1),col=c('purple', 'red', 'blue'))

# beta (we provide a weakly informative prior, 
# meaning that we know little about the transmission rate beforehand.)
# 


theta0 = list(R0 = 1, beta=1.927e-4, lambda = 1, shape = 1, gamma = 1)

## now let's run the DA-MCMC
test_6 <- run_DAMCMC_complete_smallh6(
  theta0, cp_Y, N = 100000, # you should run for twice as long to accomodate the log_lik traceplot.
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


clearcut_sim1 <- readRDS(file = "clearcut_sim1.rds")

plot(x = clearcut_sim1$t[1:5000], y = clearcut_sim1$I[1:5000], main="Epidemic Trajectory",
     xlab="Time", ylab="I(t)")
cp_Y <- observed_data(clearcut_sim1, K=12)
plot(x = cp_Y$ts, y = c(0, cp_Y$T_k), main="Observed Counts", xlab="Time", ylab="New Infections")

###

# This chunk takes in the posterior draws of (Delta, Beta) and outputs it in a tidy dataframe.
sim_df_stationary <- create_df_dbeta(test_6$dbeta_lst, test_6$beta_lst, test_6$pi01_lst, test_6$pi11_lst,
                                     test_6$beta_loglik, burnin=1, "X11")

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

## all commands are availabe outside of the menu
codamenu()

### Plot of change point location.
# configurations of delta.
barplot(analyze_dbeta(test_6$dbeta_lst)[[2]], main="Posterior Change Point Locations", 
        xlab="Time (in Weeks)", ylab="Estimated Posterior Probability",  
        border="black", col="black", names.arg=c("1","2","3","4","5","6","7","8","9","10","11"))

df_posterior_cp = data.frame(Week = seq(1,11,1),
                             cp_counts = analyze_dbeta(test_6$dbeta_lst)[[2]])

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

### Spaghetti Plot of Posterior R0
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
sim1_time_series$time_unit <- factor(sim1_time_series$time_unit, levels=unique(sim1_time_series$time_unit))

#mutate(r0)
sim1_time_series <- sim1_time_series %>%
  mutate(R0 = Beta * cp_Y$S0 / 1)


# Plot r0
# spagetti plot of R0 over time
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


## plot autocorrelations plots
autocorr.plot(coda_chain_simstudy)
### autocorrelations
autocorr.diag(coda_chain_simstudy)
## calculate effective sample size
effectiveSize(coda_chain_simstudy)


## Run 3 chains with overdispersed starting values
coda_sim1_chains = list()

# list of overdispersed starting values
R0_init = c(10, 1, 0.1)
beta_init = c(1.9e-4*10, 1.9e-4*1, 1.9e-4*0.1)

for (i in 1:3){
  print(list("current chain.." = i))
  ## now let's load our M-H and Gibbs sampler examples
  theta0 = list(R0 = R0_init[i], beta=beta_init[i], lambda = 1, shape = 1, gamma = 1)

  chain_i <- run_DAMCMC_complete_smallh6(
    theta0, cp_Y, N = 100000,
    rho = 0.2, param = "bg", approx = "ldp",
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
  ) # choosing a different population size each time?

  sim_df_i <- create_df_dbeta(chain_i$dbeta_lst, chain_i$beta_lst, 
                              chain_i$pi01_lst, chain_i$pi11_lst,
                              chain_i$beta_loglik, burnin = 1, "X11") # add log likelihood
  sim_df_i <- sim_df_i %>%
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
           theta_loglik = sapply.loglik_lst_raw..c.)

  post_samples_i <- sim_df_i %>%
    dplyr::select(Delta_1:theta_loglik)

  ## convert the output into coda then append it to mcmc list.
  coda_sim1_chains[[i]] = mcmc(post_samples_i)
  print(list("current chain end" = i))

}

coda_sim1_list = mcmc.list(coda_sim1_chains)

plot(coda_sim1_list) # codamenu()

## compute and plot Gelman-Rubin potential reduction factor
gelman.diag(coda_sim1_list)
gelman.plot(coda_sim1_list)



