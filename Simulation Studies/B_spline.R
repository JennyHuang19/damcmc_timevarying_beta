library(splines)
library(data.table)
library(ggplot2)
library(broom)


genSpline <- function(x, knots, degree, theta) {
  basis <- bs(x = x, knots = knots, degree = degree,
              Boundary.knots = c(0,12), intercept = TRUE)
  y.spline <- basis %*% theta
  dt <- data.table(x, y.spline = 2*as.vector(y.spline)) # scale by a factor of 2.
  return(list(dt = dt, basis = basis, knots = knots))
}

plot.basis <- function(basisdata) {

  dtbasis <- as.data.table(basisdata$basis)
  dtbasis[, x := seq(0, 12, length.out = .N)]
  dtmelt <- melt(data = dtbasis, id = "x",
                 variable.name = "basis", variable.factor = TRUE)
  ggplot(data=dtmelt, aes(x=x, y=value, group = basis)) +
    geom_line(aes(color=basis), size = 1) +
    theme(legend.position = "none") +
    scale_x_continuous(limits = c(0, 12),
                       breaks = c(0, basisdata$knots, 1)) +
    theme(panel.grid.minor = element_blank())
}

plot.spline <- function(basisdata, points = FALSE) {
  p <- ggplot(data = basisdata$dt)

  if (points) p <- p + geom_point(aes(x=x, y = y), color = "grey75")

  p <- p +
    geom_line(aes(x = x, y = y.spline), color = "red", size = 1) +
    scale_y_continuous(limits = c(0, 3.2)) +
    scale_x_continuous(limits = c(0, 12), breaks = seq(0,12,1)) +
    labs(y = "Transmission Rate", x = "Time (in Weeks)") +
    theme(panel.grid.minor = element_blank())+
    theme_bw() +
    theme(text = element_text(size=25), 
          axis.text.x = element_text(size = 21),
          axis.text.y = element_text(size = 21))

  return(p)

}


# Piecewise function

fx <- ifelse(x >= 0 & x < 3, 1.92,
             ifelse(x > 3 & x < 10,  1.21,
                    ifelse(x > 10 & x <= 12, 0.75, 0.75)))

ggplot(data = spline_points) +
  geom_point(data = spline_points, aes(x=time, y = fx), color = "red", size = 0.8) + 
  scale_y_continuous(limits = c(0, 3.2)) +
  scale_x_continuous(limits = c(0, 12), breaks = seq(0,12,1)) +
  labs(y = "Transmission Rate", x = "Time (in Weeks)") +
  theme(panel.grid.minor = element_blank())+
  theme_bw() +
  theme(text = element_text(size=25), 
        axis.text.x = element_text(size = 21),
        axis.text.y = element_text(size = 21))



s# quadratic spline
x <- seq(0, 12, length.out = 1200)
knots <- c(3,10)
theta = c(0.99, 0.9, 0.2, 0.9, 0.1) # red, yellow, green, blue, pink.
sdata <- genSpline(x, knots, 2, theta)
plot.basis(sdata)
# View(sdata$dt)
plot.spline(sdata)

# cubic spline
x <- seq(0, 12, length.out = 1200)
knots <- c(2.0, 2.5, 3, 3.5, 4, 9, 9.5,10, 10.5, 11)
# theta = c(0.99, 0.99, 0.99, 0.6, 0.6, 0.6, 0.6, 0.4, 0.4, 0.4) # red, yellow, green, sky, blue, pink.
theta = c(0.95, 0.96, 0.96, 0.96, 0.7, 0.65, 0.6, 0.7, 0.6, 0.36, 0.36, 0.36, 0.36, 0.34) # red, yellow, green, sky, blue, pink.
sdata <- genSpline(x, knots, 3, theta)
plot.spline(sdata)
spline_points <- sdata$dt %>%
  mutate(time = round(x, 2)) %>%
  select(time, y.spline)
# head(spline_points)
# plot.basis(sdata)
#view(sdata$dt)




# Simulate an epidemic trajectory with spline_points B(t).
change_point_data6 <- sim_smooth_cp_sem(
  S0 = 10000, I0 = 10, t_end = 12,
  spline_points = spline_points,
  theta = list(R0 = c(1.75, 1.5, 1.25, 1.00, 0.75), gamma = 1, lambda = 1, shape = 1),
  iota_dist = "exponential",
  gener = FALSE, b = 1/2,
  E0 = 0, type = "SIR"
)
# Save an object to a file
# saveRDS(change_point_data6, file = "smooth_sim.rds")
plot(x = change_point_data6$t, y = change_point_data6$I, main="Epidemic Trajectory",
     xlab="Time", ylab="I(t)")
cp_Y3 <- observed_data(change_point_data6, K=11)
plot(x = cp_Y3$ts, y = c(0, cp_Y3$T_k), main="Observed Counts", xlab="Time", ylab="New Infections")

### Design a smooth function for Beta.
# https://mycurvefit.com/

# t=seq(0, 10, length.out = 50)
# y = 0.75 + 2.44*(t+1.5) - 1.146*(t+1.5)^2 + 0.2159*(t+1.5)^3 - 0.01785*(t+1.5)^4 + 0.0005316*(t+1.5)^5
# plot(t,y/10000, ylab = "B(t)", main = "Transmission Rate Over Time")
#
# smooth_function <- function(t){
#   y = 0.75 + 2.44*(t+1.5) - 1.146*(t+1.5)^2 + 0.2159*(t+1.5)^3 - 0.01785*(t+1.5)^4 + 0.0005316*(t+1.5)^5
#   return(y)
# }
