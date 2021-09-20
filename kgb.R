set.seed(555)

possible.ns <- seq(from=100, to=2000, by=50)     # The sample sizes we'll be considering
powers_1 <- rep(NA, length(possible.ns))           # Empty object to collect simulation estimates
alpha <- 0.05                                    # Standard significance level
sims <- 500                                      # Number of simulations to conduct for each N


# ==================================================
# Effect size 1
# ==================================================

#### Outer loop to vary the number of subjects ####
for (j in 1:length(possible.ns)){
  N <- possible.ns[j]                              # Pick the jth value for N
  
  significant.experiments <- rep(NA, sims)         # Empty object to count significant experiments
  
  #### Inner loop to conduct experiments "sims" times over for each N ####
  for (i in 1:sims){
    Y0 <-  rnorm(n = N, mean = -0.17, sd = 0.49)              # control potential outcome
    tau <- 0.02-(-0.17)                             # Hypothesize treatment effect
    Y1 <- Y0 + tau                                 # treatment potential outcome
    Z.sim <- rbinom(n=N, size=1, prob=.5)          # Do a random assignment
    Y.sim <- Y1*Z.sim + Y0*(1-Z.sim)               # Reveal outcomes according to assignment
    fit.sim <- lm(Y.sim ~ Z.sim)                   # Do analysis (Simple regression)
    p.value <- summary(fit.sim)$coefficients[2,4]  # Extract p-values
    significant.experiments[i] <- (p.value <= alpha) # Determine significance according to p <= 0.05
  }
  
  powers_1[j] <- mean(significant.experiments)       # store average success rate (power) for each N
}

# ==================================================
# Effect size 2
# ==================================================
set.seed(555)
powers_2 <- rep(NA, length(possible.ns))           # Empty object to collect simulation estimates


#### Outer loop to vary the number of subjects ####
for (j in 1:length(possible.ns)){
  N <- possible.ns[j]                              # Pick the jth value for N
  
  significant.experiments <- rep(NA, sims)         # Empty object to count significant experiments
  
  #### Inner loop to conduct experiments "sims" times over for each N ####
  for (i in 1:sims){
    Y0 <-  rnorm(n = N, mean = -0.17, sd = 0.49)              # control potential outcome
    tau <- (0.02-(-0.17))/2                             # Hypothesize treatment effect
    Y1 <- Y0 + tau                                 # treatment potential outcome
    Z.sim <- rbinom(n=N, size=1, prob=.5)          # Do a random assignment
    Y.sim <- Y1*Z.sim + Y0*(1-Z.sim)               # Reveal outcomes according to assignment
    fit.sim <- lm(Y.sim ~ Z.sim)                   # Do analysis (Simple regression)
    p.value <- summary(fit.sim)$coefficients[2,4]  # Extract p-values
    significant.experiments[i] <- (p.value <= alpha) # Determine significance according to p <= 0.05
  }
  
  powers_2[j] <- mean(significant.experiments)       # store average success rate (power) for each N
}

# ==================================================
# Gather in one dataset
# ==================================================
power_df <- data.frame(sample_size = rep(possible.ns, 2),
                       power = c(powers_1, powers_2),
                       group = factor(c(rep("Fuld effektstørrelse fra Hainmueller et al. (2017)", 39), 
                                        rep("Halv effektstørrelse fra Hainmueller et al. (2017)", 39))))

power_df <- data.frame(sample_size = possible.ns,
                       power = powers_1)


scaleFUN <- function(x) sprintf("%.1f", x)


# ==================================================
# Plot it
# ==================================================
library(ggplot2)

ggplot(power_df, aes(x = sample_size, y = power)) +
  geom_segment(aes(x = 0, xend = 2000, y = 0.9, yend = 0.9), linetype = 2, size = 0.5) +
  geom_line(data = subset(power_df, group == "Fuld effektstørrelse fra Hainmueller et al. (2017)")) +
  geom_line(data = subset(power_df, group == "Halv effektstørrelse fra Hainmueller et al. (2017)")) +
  geom_point(size = 3, fill = "white", color = "black", aes(shape = group),
             alpha = 1) +
  scale_shape_manual(values = c(21, 24)) +
  labs(x = "Stikprøvestørrelse",
       y = "Power") +
  scale_x_continuous(breaks = seq(0,2000, 250),
                     labels = c("0", "250", "500", "750", "1.000", "1.250", "1.500", "1.750", "2.000"),
                     limits = c(0, 2000)) +
  scale_y_continuous(breaks = seq(0.1,1, .1),
                     labels = scaleFUN,
                     limits = c(0.1, 1)) +
  theme_bw() +
  theme(legend.position = "bottom",
        legend.title = element_blank(),
        legend.direction = "vertical")


ggplot(power_df, aes(x = sample_size, y = power)) +
  geom_segment(aes(x = 0, xend = 2000, y = 0.9, yend = 0.9), linetype = 2, size = 0.5) +
  geom_line(data = power_df) +
  geom_point(size = 3, fill = "white", color = "black",
             shape = 21,
             alpha = 1) +
  labs(x = "Stikprøvestørrelse",
       y = "Power") +
  scale_x_continuous(breaks = seq(0,2000, 250),
                     labels = c("0", "250", "500", "750", "1.000", "1.250", "1.500", "1.750", "2.000"),
                     limits = c(0, 2000)) +
  scale_y_continuous(breaks = seq(0.1,1, .1),
                     labels = scaleFUN,
                     limits = c(0.1, 1)) +
  theme_bw() +
  theme(legend.position = "bottom",
        legend.title = element_blank(),
        legend.direction = "vertical")
