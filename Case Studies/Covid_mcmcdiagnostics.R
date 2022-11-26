
plot(x = duke_Y2$ts, y = c(0, duke_Y2$T_k), main="Duke COVID-19 Spr. 2021", sub="Jan. 11-17th to Apr. 19-25th", xlab="Time (each point represents a week)", ylab="New Infection Counts")

### new plot of observed data.
df_cp_duke = data.frame(duke_Y2$ts, c(0, duke_Y2$T_k))
ggplot(df_cp_duke, aes(x=duke_Y2.ts - 0.5, y=c.0..duke_Y2.T_k.)) + 
  geom_point(size=6) +
  scale_y_continuous(limits = c(0, 250), breaks= seq(0, 250, 50)) +
  scale_x_continuous(limits = c(0, 15), breaks = seq(0,15,1)) +
  labs(title = "COVID-19 Incidence Observed", y = "New Infections", x = "Time (in Weeks)") +
  theme_bw() +
  theme(text = element_text(size=25), 
        axis.text.x = element_text(size = 21),
        axis.text.y = element_text(size = 21))+
  geom_vline(xintercept=c(6, 9, 12, 14),lwd=1,color="red") # shift to the left 0.5.
###

### Plot of Draws from the posterior.
covidbeta_over_time6 <- covid_df %>% 
  mutate(X16 = X15) %>%  # add the stopping week
  dplyr::select(ID, X1.y:X16) %>% 
  dplyr::select(-sapply.pi01_lst_raw..c., -sapply.pi11_lst_raw..c., -sapply.loglik_lst_raw..c.)

covidbeta_over_time6 <- covidbeta_over_time6  %>% 
  pivot_longer(cols = -ID,
               names_to = "time_unit", names_prefix = "V",
               # names_transform = list(Time = as.integer),
               values_to = "Beta") %>%
  filter(ID %% 10 == 0)

#preserve the order of time units.
covidbeta_over_time6$time_unit <- factor(covidbeta_over_time6$time_unit, levels=unique(covidbeta_over_time6$time_unit))

#mutate(r0)
covidbeta_over_time6 <- covidbeta_over_time6 %>%
  mutate(R0 = Beta / gamma * duke_Y2$S0)

# Save as PDF
pdf(file = "Covid_beta_posterior.pdf",   # The directory you want to save the file in
    width = 5, # The width of the plot in inches
    height = 6) # The height of the plot in inches

covidbeta_over_time6 %>% 
  ggplot(aes(x = time_unit, y = R0, group = ID)) +
  geom_step(alpha = 0.1, color = "red") +
  labs(y = "Transmission Rate", x = "Weeks (Jan. 11th - April 19th)", color='Change Point Configuration') + 
  scale_x_discrete(labels=seq(0,16,1)) + # add week 6 in.
  expand_limits(y = 0)+
  theme_bw()+
  theme(text = element_text(size=25), 
        axis.text.x = element_text(size = 21),
        axis.text.y = element_text(size = 21))
###

# This chunk takes in the posterior draws of (Delta, Beta) and outputs it in a tidy dataframe.
covid_df <- create_df_dbeta(dk_observations6$dbeta_lst, dk_observations6$beta_lst, dk_observations6$pi01_lst, dk_observations6$pi11_lst, dk_observations6$beta_loglik, burnin = 1, "X14")
cov_df <- covid_df %>%
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
         Beta_12 = X12.y, Delta_12 = X12.x,
         Beta_13 = X13.y, Delta_13 = X13.x,
         Beta_14 = X14.y, Delta_14 = X14.x,
         Beta_15 = X15,
         pi_01 = sapply.pi01_lst_raw..c.,
         pi_11 = sapply.pi11_lst_raw..c.,
         theta_loglik = sapply.loglik_lst_raw..c.)


### the top configurations in MCMC samples.
cov_df %>%
  count(delta_beta) %>%
  mutate(freq = n / sum(n)) %>%
  arrange(desc(n)) %>%
  head(5) %>%
  ggplot(data = ., aes(x = delta_beta, y = freq)) +
  geom_bar(stat = "identity") +
  labs(title = "Frequency of Change Point Configurations", y = "Frequency", x = "Configuration")+
  theme(text = element_text(size=10), axis.text.x = element_text(angle = 90))
###


### Change point posteriors

df_posterior_covid = data.frame(Week = seq(1,14,1),
                                cp_counts = analyze_dbeta(dk_observations6$dbeta_lst)[[2]])

ggplot(df_posterior_covid, aes(Week, cp_counts)) +
  geom_col(fill="black") +
  scale_y_continuous(limits = c(0, 1), breaks= seq(0, 1, 0.2)) +
  scale_x_continuous(limits = c(0, 15), breaks = seq(1,14,1)) +
  labs(title = "Posterior Change Point Locations", 
       y = "Posterior Probability", x = "Time (in Weeks)") +
  theme_bw() +
  theme(text = element_text(size=23), 
        axis.text.x = element_text(size = 19),
        axis.text.y = element_text(size = 19))

###

###
covid_delta_beta <- cov_df %>% # get rid of delta_beta column. we won't need it from here on.
  mutate(Beta_15 = Beta_14) %>%
  select(ID, Delta_1:Beta_15, theta_loglik)


## convert the output into coda format
coda_chain_covid= mcmc(covid_delta_beta)

coda_chain_covid[1,] # look at the first row, all columns.

summary(coda_chain_covid) # the summary of each variable of interest.

## we can also compute Monte Carlo error for the quantiles of each variable in the coda_gibbs_chain mat # ola.
mcse.q.mat(coda_chain_covid, 0.025)
mcse.q.mat(coda_chain_covid, 0.975)

## plot traceplots and posterior densities
plot(coda_chain_covid)

## look at what coda can do
help(package=coda)

## look at the menu options
codamenu()

## all commands are availabe outside of the menu

## plot autocorrelations plots
autocorr.plot(coda_chain_covid)
### autocorrelations
autocorr.diag(coda_chain_covid)
## calculate effective sample size
effectiveSize(coda_chain_covid)

