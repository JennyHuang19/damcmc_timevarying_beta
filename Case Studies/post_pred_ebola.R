#
# predictive trajectory (one draw.)

profvis::profvis({
  ebola_pp <- sim_posterior_predictive(S0 = 3000000, I0 = 10, t_end = 12,
                                  theta = list(R0 = ebola_R0, gamma = 0.8, lambda = 1, shape = 1),
                                  change_day = seq(1,11,1),
                                  iota_dist = "exponential",
                                  gener = FALSE, b = 1/2,
                                  E0 = 0, type = "SIR") # 291823
})

ebola_pp$beta_MLE
ebola_pp$final_size # 10,000 vs. 50,000
mean(ebola_pp$mu_j_list) # 50,000 = 191.3949 # 300000  = 170.7944
length(ebola_pp$cumu_cases)


#
# Initialize an empty list to store the posterior predictive trajectories and the final_size of each trajectory.
pred_trajectories <- vector("list", 100)
final_sizes <- vector(length = 100)
mu_j_avgs <- vector("list", 100)

# Draw samples from P(y_hat | theta).
j = 1 # counter to index the posterior samples.
for(i in 1:100){

  print(j)
  # R0_j = ebola_6$beta_lst[[j+1000]]*S0/gamma # sample a draw from the posterior distribution.
  pred_trajectories[[i]] <- sim_posterior_predictive(S0 = 300000, I0 = 5, t_end = 12,
                                                     theta = list(R0 = ebola_6$beta_lst[[j+2000]]*S0/gamma, gamma = 0.8, lambda = 1, shape = 1),
                                                     change_day = seq(1,11,1),
                                                     iota_dist = "exponential",
                                                     gener = FALSE, b = 1/2,
                                                     E0 = 0, type = "SIR")

  final_sizes[i] <- pred_trajectories[[i]]$final_size
  j <- j + 25

}

# Histogram of final size. Bayesian p-value (fraction of simulations greater than the observed data)
hist(final_sizes, breaks = 20, col = "blue",
     xlab="Final Size",main="Final Size of Simulated Ebola Epidemics")
abline(v = sum(Y_ebola$T_k), col="black", lwd=4, lty=2)
#
obs_final_size = sum(Y_ebola$T_k)
bayesian_pval = length(final_sizes[final_sizes>obs_final_size])/length(final_sizes)
bayesian_pval # 0.42
# Kernel Density Plot
d <- density(final_sizes) # returns the density data
plot(d) # plots the results
abline(v = sum(Y_ebola$T_k), col="black", lwd=4, lty=2)

#
### Plot Cumulative cases over time
plot(x = Y_ebola$ts, y = c(0, cumsum(Y_ebola$T_k)), main="Cumulative Cases Over Time", col="red",lwd=2,
     xlab = "Time (months)", ylab = "Cumulative Cases", ylim=c(0,1000))

for(i in 2:100){
  # posterior draws
  lines(x = pred_trajectories[[i]]$cumu_cases_time, y = pred_trajectories[[i]]$cumu_cases, type="l",
        col = "blue", lwd=0.5, main="Cumulative Cases Over Time",
        xlab = "Time (months)", ylab = "Cumulative Cases")
} # much variability.



### Plot incidence cases over time w observed data.
plot(x = Y_ebola$ts, y = c(0, Y_ebola$T_k), main="New Cases Over Time", col="red",lwd=3,
     xlab = "Time (months)", ylab = "Observed New Cases")
T_k_list <- list()
for(i in 100:1900){
  # draw from posterior
  ebola_R0 <- ebola_6$beta_lst[[i*10]]*291823 # N=5000 for the posterior samples. num iter.= 5000
  
  
  # simulate
  trajectory.i <- simSIR.Markov(N=291823, I0=5, t_end=12, beta=ebola_R0*gamma, gamma = 0.8)
  # plot
  T_k <- c(0, diff(trajectory.i$cumu_cases_T_k)) # the difference in cumulative cases (if the epidemic dies out, fill the remaining vector with 0's.)
  Observation_times <- trajectory.i$T_k_times # times at which it crosses an integer for the first time.
  # print(T_k)
  T_k_list <- lappend(T_k_list, T_k)
  # lines(Observation_times, T_k, col="blue")
  print(i*10)
}
###



### Function to append any vector to a list. Keep track of T_k at each iter.
lappend <- function (lst, ...){
  lst <- c(lst, list(...))
  return(lst)
}

# Post processing to get list of [2.5%, 97.5%] summary stats. (https://stackoverflow.com/questions/31399954/function-for-same-length-of-vectors-in-r)
maxlength <- max(sapply(T_k_list, length))

T_k_list <- sapply(T_k_list, FUN = function(x, ml) {
  difference <- ml - length(x)
  c(x, rep(0, difference))
}, ml = maxlength, simplify = FALSE)


# List of lists into a matrix for fbplots.
mat_Tk <- do.call("cbind",T_k_list)
dim(mat_Tk) # n x p = 10 x 151

#functional boxplot
fbplot(mat_Tk, x = seq(0,12,1), method='MBD',
       color = "light blue",
       ylim=c(0,60), xlim = c(0,12), outliercol="white", 
       cex.axis = 1.5, cex.lab=1.25, cex.main=1.9,
       main = "Ebola",
       xlab = "Time (in Months)", ylab = "New Cases") # unable to deal with starting at 0. unable to take 95% CI.
#observed case data.
points(Y_ebola$ts, c(0,Y_ebola$T_k), col = "black", lwd=4, lty=1)


# Take the element-wise summary statistics of the list of lists (simulated trajectories).

stats_sim_trajectory=apply(simplify2array(T_k_list), 1, quantile, probs = c(.025, .975)) # take column means.
mean_sim_trajectory=apply(simplify2array(T_k_list), 1, mean)
sd_sim_trajectory=apply(simplify2array(T_k_list), 1, sd)

sim_trajectory_2.5=stats_sim_trajectory[1,]
sim_trajectory_97.5=stats_sim_trajectory[2,]


### Plot incidence cases over time and observed data.
plot(x = Y_ebola$ts, y = c(0, Y_ebola$T_k), col="red",lwd=3,
     xlab = "Time (months)", ylab = "Observed New Cases", ylim=c(0,100))
lines(seq(0,12,1), sim_trajectory_2.5, col="blue") # note: the observation times for the simulated trajectories are chosen as the first time after each integer.   .
lines(seq(0,12,1), mean_sim_trajectory, col="black")
lines(seq(0,12,1), sim_trajectory_97.5, col="blue")




