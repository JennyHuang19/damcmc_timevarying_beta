###


### Plot incidence cases over time w observed data.
plot(x = duke_Y2$ts, y = c(0, duke_Y2$T_k), main="New Cases Over Time", col="red",lwd=3,
     xlab = "Time (weeks)", ylab = "Observed New Cases")
T_k_list <- list()
final_sizes <- list()
for(i in 1:2500){
  # draw from posterior
  covid_R0 <- dk_observations6$beta_lst[[i+1000]]*duke_Y2$S0 # num iter.= 5000


  # simulate
  trajectory.i <- simSIR.Markov(N=duke_Y2$S0, I0=duke_Y2$I0, t_end=duke_Y2$t_end, beta=covid_R0*gamma, gamma = 1)
  # plot
  T_k <- c(0, diff(trajectory.i$cumu_cases_T_k)) # the difference in cumulative cases (if the epidemic dies out, fill the remaining vector with 0's.)
  Observation_times <- trajectory.i$T_k_times # times at which it crosses an integer for the first time.
  # print(T_k)
  T_k_list <- lappend(T_k_list, T_k)
  final_sizes <- append(final_sizes, trajectory.i$final.size)
  # lines(Observation_times, T_k, col="blue")
  print(i)
}
###


# Histogram of final size. Bayesian p-value (fraction of simulations greater than the observed data)
hist(unlist(final_sizes), breaks = 20, col = "blue",
     xlab="Final Size",main="Final Size of Simulated Covid Outbreaks")
abline(v = sum(duke_Y2$T_k), col="black", lwd=4, lty=2)
#
obs_final_size = sum(duke_Y2$T_k)
bayesian_pval = length(final_sizes[final_sizes<obs_final_size])/length(final_sizes)
bayesian_pval



### Function to append any vector to a list. Keep track of T_k at each iter.
lappend <- function (lst, ...){
  lst <- c(lst, list(...))
  return(lst)
}
a <- list()
a <- lappend(a, c(1,2,3))
a <- lappend(a, c(7,8,9,10))
a

# Post processing to get list of [2.5%, 97.5%] summary stats. (https://stackoverflow.com/questions/31399954/function-for-same-length-of-vectors-in-r)
maxlength <- max(sapply(T_k_list, length))

T_k_list <- sapply(T_k_list, FUN = function(x, ml) {
  difference <- ml - length(x)
  c(x, rep(0, difference))
}, ml = maxlength, simplify = FALSE)

# Take the element-wise summary statistics of the list of lists (simulated trajectories).

stats_sim_trajectory=apply(simplify2array(T_k_list), 1, quantile, probs = c(.025, .975)) # take column means.
mean_sim_trajectory=apply(simplify2array(T_k_list), 1, mean)
sd_sim_trajectory=apply(simplify2array(T_k_list), 1, sd)

sim_trajectory_2.5=stats_sim_trajectory[1,]
sim_trajectory_97.5=stats_sim_trajectory[2,]


### Plot incidence cases over time and observed data.
plot(x = duke_Y2$ts, y = c(0, duke_Y2$T_k), col="red",lwd=3,
     xlab = "Time (months)", ylab = "Observed New Cases", ylim=c(0,500))
lines(seq(0,15,1), sim_trajectory_2.5, col="blue") # note: the observation times for the simulated trajectories are chosen as the first time after each integer.   .
lines(seq(0,15,1), mean_sim_trajectory, col="black")
lines(seq(0,15,1), sim_trajectory_97.5, col="blue")

