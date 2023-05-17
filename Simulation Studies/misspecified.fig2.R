# Run sim_study_1 with cp location fixed but fixed incorrectly. 3 --> 2, 10 --> 9.
# (may 16th) run existing methods and compare on metrics of final size and peak infection numbers.

### Run MCMC.
test_6.fig2A <- run_DAMCMC_complete_smallh6.2(
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
  fixed_delta = c(0,0,0,0,0,0,0,0,0,0,0)
)

test_6.fig2B <- run_DAMCMC_complete_smallh6.2(
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
  fixed_delta = c(0,0,0,0,0,0,0,0,1,0,0)
)

test_6.fig2C <- run_DAMCMC_complete_smallh6.2(
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

test_6.fig2D <- run_DAMCMC_complete_smallh6(
  theta0, cp_Y, N = 50000,
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

# create marginal plot
print(analyze_dbeta(test_6.fig2C$dbeta_lst))

barplot(analyze_dbeta(test_6.fig2C$dbeta_lst)[[2]], main="Posterior Change Point Locations", xlab="Change Point Location",
        ylab="Proportion of MCMC Samples",  border="blue", names.arg=c("1","2","3","4","5","6","7","8","9","10","11"),ylim=c(0,1))

### 2. Compare on final size and peak infection count.
# Simulate 1 draw from the predictive distribution

### Posterior Predictive Check on Covid-19 Case Study.
T_k_list <- list()
final_sizes <- list()
for(i in 100:500){
  # draw from posterior
  R0 <- test_6.fig2C$beta_lst[[i*10]]*cp_Y$S0
  
  
  # simulate
  trajectory.i <- simSIR.Markov(N=cp_Y$S0, I0=cp_Y$I0, t_end=cp_Y$t_end, beta=R0, gamma = 1) # gamma has big influence.
  # print(trajectory.i)
  # plot
  T_k <- c(0, diff(trajectory.i$cumu_cases_T_k)) # the difference in cumulative cases (if the epidemic dies out, fill the remaining vector with 0's.)
  Observation_times <- trajectory.i$T_k_times # times at which it crosses an integer for the first time.
  # print(T_k)
  T_k_list <- lappend(T_k_list, T_k)
  final_sizes <- append(final_sizes, trajectory.i$final.size)
  # lines(Observation_times, T_k, col="blue")
  print(i*10)
}
###


# Histogram of final size. Bayesian p-value (fraction of simulations greater than the observed data)
hist(unlist(final_sizes), breaks = 20, col = "blue",
     xlab="Final Size",main="Final Size of Simulated Outbreaks")
abline(v = sum(cp_Y$T_k), col="black", lwd=4, lty=2)
#
obs_final_size = sum(cp_Y$T_k)
bayesian_pval = length(final_sizes[final_sizes<obs_final_size])/length(final_sizes)
bayesian_pval # 0.304


# Post processing to get cases over time.
maxlength <- max(sapply(T_k_list, length)) # ndatapoints

T_k_list <- sapply(T_k_list, FUN = function(x, ml) { # for formatting, add 0's to the epidemics that die out early.
  difference <- ml - length(x)
  c(x, rep(0, difference))
}, ml = maxlength, simplify = FALSE)

# List of lists into a matrix for fbplots.
mat_Tk <- do.call("cbind",T_k_list)
dim(mat_Tk) # ndatapoints x nsimulations = 16 x 400

#functional boxplot
fbplot(mat_Tk,x = seq(0,12,1), method='MBD',
       color = "light blue",
       ylim=c(0,400), xlim = c(0,12), outliercol="white", 
       cex.axis = 1.5, cex.lab=1.25, cex.main = 1.9,
       main = "Fixed Change Point Model",
       xlab = "Time (in Weeks)", ylab = "New Cases") # unable to deal with starting at 0. unable to take 95% CI.
#add observed case data.
points(cp_Y$ts, c(0,cp_Y$T_k), col = "black", lwd=4, lty=2)
