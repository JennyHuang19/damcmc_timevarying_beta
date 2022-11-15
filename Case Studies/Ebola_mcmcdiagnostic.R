## This script runs mcmc diagnostics for simulation study 1.

## first we need to load coda and mcmcse packages (you need to install them first)
library(coda)
library(mcmcse)

theta0 = list(R0 = 1, beta=3e-7, lambda = 1, shape = 1, gamma = 1)

## now let's load our M-H and Gibbs sampler examples
ebola_6 <- run_DAMCMC_complete_smallh6(
  theta0, Y_ebola, N = 50000,
  rho = 0.10, param = "bg", approx = "ldp",
  iota_dist = "exponential", gener = FALSE, b = 1/2,
  thin = 10, plt = 6000,
  par_prior = list(
    a_beta = 0.1, b_beta = 0.1,
    a_gamma = 1, b_gamma = 1,
    a_R0 = 2, b_R0 = 2,
    a_lambda = 1, b_lambda = 1,
    a_pi11 = .5, b_pi11 = .5,
    a_pi01 = 5, b_pi01 = 5
  ),
  length_delta = 11,
  x_b_ratio = 2,
  gamma=0.8
)

ebola_6$theta_rate_accept # theta = (dbeta, beta).
ebola_6$x_rate_accept # latent data.

print(analyze_dbeta(ebola_6$dbeta_lst))
barplot(analyze_dbeta(ebola_6$dbeta_lst)[[2]], main="Delta Configurations", xlab="Change Point Location",
        ylab="Proportion of MCMC Samples",names.arg = seq(1,11),
        border="blue")
# abline(h=1/2, col="blue", lty="dashed", lwd=4)

# This chunk takes in the posterior draws of (Delta, Beta) and outputs it in a tidy dataframe.
ebola_df <- create_df_dbeta(ebola_6$dbeta_lst, ebola_6$beta_lst, ebola_6$beta_loglik, burnin = 1, "X11")
ebola_df <- ebola_df %>%
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
         theta_loglik = sapply.loglik_lst_raw..c.)


### the top configurations in MCMC samples.
ebola_df %>%
  count(delta_beta) %>%
  mutate(freq = n / sum(n)) %>%
  arrange(desc(n)) %>%
  head(5) %>%
  ggplot(data = ., aes(x = delta_beta, y = freq)) +
  geom_bar(stat = "identity") +
  labs(title = "Frequency of Change Point Configurations", y = "Frequency", x = "Configuration")+
  theme(text = element_text(size=10), axis.text.x = element_text(angle = 90))
###

###
ebola_delta_beta <- ebola_df %>% # get rid of delta_beta column. we won't need it from here on.
  mutate(Beta_13 = Beta_12) %>%
  select(ID, Beta_1:Beta_13, -theta_loglik)

time_series_df <- ebola_delta_beta  %>%
  pivot_longer(cols = -ID,
               names_to = "time_unit", names_prefix = "V",
               # names_transform = list(Time = as.integer),
               values_to = "Beta") %>%
  filter(ID %% 10 == 0)

#preserve the order of time units.
time_series_df$time_unit <- factor(time_series_df$time_unit, levels=unique(time_series_df$time_unit))

#mutate(r0)
time_series_df <- time_series_df %>%
  mutate(R0 = Beta * Y_ebola$S0 / 0.8)


# Plot r0
# spagetti plot of R0 over time
time_series_df %>%
  ggplot(aes(x = time_unit, y = R0, group = ID)) +
  geom_step(alpha = 0.1, color = "red") +
  labs(title = "Reproduction number of Ebola", y = "R0", x = "Month (Jan-Dec. 2014)", color='Change Point Configuration') +
  scale_x_discrete(labels=seq(1,13,1)) +
  expand_limits(y = 0) +
  theme_bw()
###




## convert the output into coda format
coda_chain_ebola = mcmc(ebola_delta_beta)

coda_chain_ebola[1,] # look at the first row, all columns.

summary(coda_chain_ebola) # the summary of each variable of interest.

## we can also compute Monte Carlo error for the quantiles of each variable in the coda_gibbs_chain mat # ola.
mcse.q.mat(coda_chain_ebola, 0.025)
mcse.q.mat(coda_chain_ebola, 0.975)

## plot traceplots and posterior densities
plot(coda_chain_ebola)

## look at what coda can do
help(package=coda)

## look at the menu options
codamenu()

## all commands are availabe outside of the menu

## plot autocorrelations plots
autocorr.plot(coda_chain_ebola)
### autocorrelations
autocorr.diag(coda_chain_ebola)
## calculate effective sample size
effectiveSize(coda_chain_ebola)


