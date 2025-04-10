library(ltm)
library(difR)

set.seed(123)

# Group sizes
n_total <- 1000
ratio <- c(1, 1)
group_sizes <- function(n_total, ratio) {
  total_parts <- sum(ratio)
  n2 <- round(n_total * ratio[1] / total_parts)  # focal group
  n1 <- n_total - n2                             # reference group
  return(list(n1 = n1, n2 = n2))
}
sizes <- group_sizes(n_total, ratio)
n1 <- sizes$n1
n2 <- sizes$n2

# Abilities
theta1 <- rnorm(n1, mean = 0, sd = 1)    # Reference
theta2 <- rnorm(n2, mean = -1, sd = 1)   # Focal

# Item parameters (5 items)
generate_param <- function(I) {
  b <- rnorm(I, mean = 0, sd = 1)
  z <- rnorm(I, mean = 0, sd = sqrt(0.1225))
  a <- exp(z)
  c <- rep(0.2, I)
  param <- cbind(b, a, c)
  rownames(param) <- paste0("Item", 1:I)
  colnames(param) <- c("b", "a", "c")
  return(param)
}
I <- 5
params <- generate_param(I)

# Data generation
data1 <- rmvlogis(n1, params, IRT = TRUE, link = "logit", z.vals = matrix(theta1))
data2 <- rmvlogis(n2, params, IRT = TRUE, link = "logit", z.vals = matrix(theta2))
data <- rbind(data1, data2)
colnames(data) <- paste0("Item", 1:I)  # make sure names match expectations
group <- c(rep(0, n1), rep(1, n2))

# Run SIBTEST (by default, first item is tested)
res_sib <- difSIBTEST(data, group, focal.name = 1, purify = FALSE)[3]

difMH(data, group, focal.name = 1, correct = TRUE, MHstat = "MHChisq", exact = FALSE)
difSIBTEST(data, group, focal.name = 1,, purify = FALSE)$stat
difSIBTEST(data, group, focal.name = 1,, purify = FALSE)[[3]][1]
difMH(data, group, focal.name = 1, correct = TRUE, MHstat = "MHChisq", exact = FALSE)[[1]][1]
res_bd <- difBD(data, group, focal.name = 1, BDstat = "trend")[[1]][1, 3]
