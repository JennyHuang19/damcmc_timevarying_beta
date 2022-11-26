### Posterior Predictive Check on Covid-19 Case Study.

### Helper Function to append any vector to a list. Keep track of T_k at each iter.
lappend <- function (lst, ...){
  lst <- c(lst, list(...))
  return(lst)
}
a <- list()
a <- lappend(a, c(1,2,3))
a <- lappend(a, c(7,8,9,10))
a


### Plot incidence cases over time w observed data.
plot(x = duke_Y2$ts, y = c(0, duke_Y2$T_k), main="New Cases Over Time", col="red",lwd=3,
     xlab = "Time (weeks)", ylab = "Observed New Cases")
T_k_list <- list()
# final_sizes <- list()
for(i in 100:500){
  # draw from posterior
  covid_R0 <- dk_observations6$beta_lst[[i*10]]*duke_Y2$S0 # num iter.= 5000


  # simulate
  trajectory.i <- simSIR.Markov(N=duke_Y2$S0, I0=duke_Y2$I0, t_end=duke_Y2$t_end, beta=covid_R0*0.95, gamma = 1) # gamma has big influence.
  # print(trajectory.i)
  # plot
  T_k <- c(0, diff(trajectory.i$cumu_cases_T_k)) # the difference in cumulative cases (if the epidemic dies out, fill the remaining vector with 0's.)
  Observation_times <- trajectory.i$T_k_times # times at which it crosses an integer for the first time.
  # print(T_k)
  T_k_list <- lappend(T_k_list, T_k)
  # final_sizes <- append(final_sizes, trajectory.i$final.size)
  # lines(Observation_times, T_k, col="blue")
  print(i*10)
}
###


# Histogram of final size. Bayesian p-value (fraction of simulations greater than the observed data)
hist(unlist(final_sizes), breaks = 20, col = "blue",
     xlab="Final Size",main="Final Size of Simulated Covid Outbreaks")
abline(v = sum(duke_Y2$T_k), col="black", lwd=4, lty=2)
#
obs_final_size = sum(duke_Y2$T_k)
bayesian_pval = length(final_sizes[final_sizes<obs_final_size])/length(final_sizes)
bayesian_pval # 0.304


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
fbplot(mat_Tk,x = seq(0,15,1), method='MBD',
       color = "light blue",
       ylim=c(0,1400), xlim = c(0,15), outliercol="white", 
       cex.axis = 1.5, cex.lab=1.25, cex.main = 1.9,
       main = "Covid-19",
       xlab = "Time (in Weeks)", ylab = "New Cases") # unable to deal with starting at 0. unable to take 95% CI.
#observed case data.
points(duke_Y2$ts, c(0,duke_Y2$T_k), col = "black", lwd=4, lty=2)

### OLD CODE BELOW.

# Take the element-wise summary statistics of the list of lists (simulated trajectories).
stats_sim_trajectory=apply(simplify2array(T_k_list), 1, quantile, probs = c(.025, .975)) # take column means.
mean_sim_trajectory=apply(simplify2array(T_k_list), 1, mean)
sd_sim_trajectory=apply(simplify2array(T_k_list), 1, sd)


sim_trajectory_2.5=stats_sim_trajectory[1,]
sim_trajectory_97.5=stats_sim_trajectory[2,]

### Plot incidence cases over time and observed data.
plot(x = duke_Y2$ts, y = c(0, duke_Y2$T_k), col="red",lwd=3,
     xlab = "Weeks (Jan. 11th - April 19th)", ylab = "Observed New Cases", ylim=c(0,500))

plot(x = duke_Y2$ts, y = c(0, duke_Y2$T_k), col="red",lwd=3,
     xlab = "Weeks (Jan. 11th - April 19th)", ylab = "Observed New Cases", ylim=c(0,500))
lines(seq(0,15,1), sim_trajectory_2.5, col="blue") # note: the observation times for the simulated trajectories are chosen as the first time after each integer.   .
lines(seq(0,15,1), mean_sim_trajectory, col="black")
lines(seq(0,15,1), sim_trajectory_97.5, col="blue")

# 2.5% trajectory
temp_lo <- loess(sim_trajectory_2.5~seq(0,15,1))
smooth_Dist2.5 <- predict(temp_lo, seq(0,15,1))
lines(seq(0,15,1),smooth_Dist2.5, col="brown")
# Mean trajectory
temp_lo <- loess(mean_sim_trajectory~seq(0,15,1), span=0.50)
smooth_mean <- predict(temp_lo, seq(0,15,1))
lines(seq(0,15,1),smooth_mean, col="brown")



# Mean trajectory
mean_trajectory = smooth.spline(seq(151,1), mean_sim_trajectory, spar=0)
plot(seq(151,1),mean_sim_trajectory)
lines(mean_trajectory)
# 2.5% trajectory
two.five = smooth.spline(seq(151,1), sim_trajectory_2.5, spar=0)
plot(seq(151,1),sim_trajectory_2.5)
lines(two.five)
# 97.5% trajectory
ninetyseven.five = smooth.spline(seq(151,1), sim_trajectory_97.5, spar=0)
plot(seq(151,1),sim_trajectory_97.5)
lines(ninetyseven.five)




smooth_Dist <- matrix(NA, length(mat_Tk[,1]), length(mat_Tk[1,])) # t by IDs

for(i in 1:length(mat_Tk[1,])){
  
  temp <-  mat_Tk[,i]
  
  temp_lo <- loess(mat_Tk[,i]~seq(0,15,1))
  
  smooth_Dist[,i] <- predict(temp_lo, seq(0,15,1))
  
}

#the following adds zeros to all the endpoints of trajectories, where the interpolant becomes negative or gives NA
for(i in 1:length(IDs)){
  smallestRow <-  which.min(smooth_Dist[,i])
  if(smallestRow<length(t)){ 
    if(smooth_Dist[smallestRow,i] > 0.00001){smooth_Dist[(smallestRow+1),i] <- 0} else {smooth_Dist[smallestRow,i] <- 0}} 
}

smooth_Dist[smooth_Dist < 0] <- 0 #replaces negative smoothed values just with zero
smooth_Dist[is.na(smooth_Dist)] <- 0 #replace NA with zero: not sure if this is correct





