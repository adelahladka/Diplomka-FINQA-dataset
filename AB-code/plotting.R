#Libraries-----------------------------------------------------------------------
rm(list = ls())
library(ggplot2)
library(dplyr)
library(tidyr)
library(gridExtra)
################################################################################
################################################################################
#A. Theoretical graphs


#Introduction-------------------------------------------------------------------
# ICC function
logit <- function(x) 1 / (1 + exp(-x))
icc <- function(theta, a, b, c) c + (1 - c) * logit(a * (theta - b))

# Set constants
theta <- seq(-4, 4, length.out = 500)
a0 <- 1
b0 <- 0
c <- 0.2

# ----- UNIFORM DIF (varying d) -----
d_vals <- c(0.5, 1, 2, 4)
uniform_data <- lapply(d_vals, function(d) {
  data.frame(
    theta = theta,
    prob = icc(theta, a0, b0 + d, c),
    group = paste0("d = ", d)
  )
}) %>% bind_rows() %>%
  mutate(type = "Focal")

ref_uniform <- data.frame(
  theta = theta,
  prob = icc(theta, a0, b0, c),
  group = "Reference",
  type = "Reference"
)

plot_data_uniform <- bind_rows(uniform_data, ref_uniform)

# ----- NON-UNIFORM DIF (varying Δ) -----
delta_vals <- c(0.4, 0.6, 0.8, 1.0)
compute_a1 <- function(a0, c, delta) {
  numerator <- 2 * a0 * log(2)
  denominator <- 2 * log(2) + (delta * a0) / (1 - c)
  numerator / denominator
}

nonuniform_data <- lapply(delta_vals, function(delta) {
  a1 <- compute_a1(a0, c, delta)
  data.frame(
    theta = theta,
    prob = icc(theta, a1, b0, c),
    group = paste0("Δ = ", delta)
  )
}) %>% bind_rows() %>%
  mutate(type = "Focal")

ref_nonuniform <- data.frame(
  theta = theta,
  prob = icc(theta, a0, b0, c),
  group = "Reference",
  type = "Reference"
)

plot_data_nonuniform <- bind_rows(nonuniform_data, ref_nonuniform)

# ----- PLOTTING -----

# Uniform DIF plot
p1 <- ggplot(plot_data_uniform, aes(x = theta, y = prob, color = group, linetype = type)) +
  geom_line(size = 1.2) +
  labs(title = "Uniform DIF: Shift in Difficulty",
       x = expression(theta), y = "P(θ)") +
  theme_minimal() +
  scale_color_brewer(palette = "Set1")

# Non-uniform DIF plot
p2 <- ggplot(plot_data_nonuniform, aes(x = theta, y = prob, color = group, linetype = type)) +
  geom_line(size = 1.2) +
  labs(title = "Non-Uniform DIF: Change in Discrimination",
       x = expression(theta), y = "P(θ)") +
  theme_minimal() +
  scale_color_brewer(palette = "Set2")

# Display both
p1
p2
################################################################################
###############################################################################¨
#B.Simulations plots

#1. DIF when parameter n is varied---------------------------------------------
#1.1 UNIFORM

#Data loading
load("results/variable-n-res-U.RData")
load("results/variable-n-U.RData")


#Extracting information from tab
tab_df <- as.data.frame(matrix(unlist(tab), nrow = 8,ncol=45, byrow = TRUE))
tab_df<- t(tab_df)
notimes <- nrow(tab_df) / 3
methods <- rep(c("Mantel", "Log", "SIB"), notimes)
colnames(tab_df) <- c('N', 'n', 'I', 'type', 'gamma', 'dif-eff', 'power', 'rejection')
tab_df <- data.frame(tab_df)
tab_df$method <- methods

small_tab <- tab_df[, c("n", "power", "rejection", 'method')]


#Plotting
nUniformPower <- ggplot(small_tab, aes(x = as.numeric(n), y = as.numeric(power), color = method)) +
  geom_line() +
  geom_point() +
  labs(x = "n", y = "Power", color = "Method") +
  theme_minimal()

nUniformRejection <- ggplot(tab_df, aes(x = as.numeric(n), y = as.numeric(rejection), color = method)) +
  geom_line() +
  geom_point() +
  labs(x = "n", y = "Power", color = "Method") +
  ylim(0, 0.1) + #To see better what is happening
  theme_minimal()

#Saving-----------------------------------------------------------------------
# — save individually —
ggsave(
  filename = "plots/nUniformPower.png",
  plot     = nUniformPower,
  width    = 6, 
  height   = 4, 
  dpi      = 300
)

ggsave(
  filename = "plots/nUniformRejection.png",
  plot     = nUniformRejection,
  width    = 6, 
  height   = 4, 
  dpi      = 300
)


# turn them into a single grob
nUniform_combined <- grid.arrange(
  nUniformPower,
  nUniformRejection,
  ncol = 2  # two columns side by side
)

# save the combined plot
ggsave(
  filename = "plots/nUniform-combined.png",
  plot     = nUniform_combined,
  width    = 12,   # adjust as you like
  height   = 6,    # adjust as you like
  dpi      = 300
)
nUniform_combined


# Itemwise plotting ------------------------------------------------------

# 1. get unique n’s (as characters, to index res_total)
unique_ns <- unique(tab_df$n)

# 2. for each n, run the function and grab only method + vote_ratio
results_extracted <- do.call(rbind, lapply(unique_ns, function(n_val) {
  # compute
  tmp <- calculate_itemwise_detection(res_total[[ as.character(n_val) ]])
  # select and add an 'n' column
  data.frame(
    n          = n_val,
    method     = tmp$method,
    vote_ratio = tmp$vote_ratio,
    stringsAsFactors = FALSE
  )
}))

# 3. inspect
head(results_extracted)

#same for one item
ggplot(results_extracted, aes(x = as.numeric(n), y = as.numeric(vote_ratio), color = method)) +
  geom_line() +
  geom_point() +
  labs(x = "n", y = "Vote Ratio", color = "Method") +
  theme_minimal()


#-------------------------------------------------------------------------------
#1.2 NON-UNIFORM
rm(list = ls())
#Data loading
load("results/variable-n-res-NU.RData")
load("results/variable-n-NU-U.RData")


#Extracting information from tab
tab_df <- as.data.frame(matrix(unlist(tab), nrow = 8,ncol=90, byrow = TRUE))
tab_df<- t(tab_df)
notimes <- nrow(tab_df) / 6
methods <- rep(c("Mantel","MantelNormal","MantelLow", "MantelHigh", "Log", "CSIB"), notimes)
colnames(tab_df) <- c('N', 'n', 'I', 'type', 'gamma', 'dif-eff', 'power', 'rejection')
tab_df <- data.frame(tab_df)
tab_df$method <- methods

small_tab <- tab_df[, c("n", "power", "rejection", 'method')]


#Plotting
nNUniformPower <- ggplot(small_tab, aes(x = as.numeric(n), y = as.numeric(power), color = method)) +
  geom_line() +
  geom_point() +
  labs(x = "n", y = "Power", color = "Method") +
  theme_minimal()

nNUniformRejection <- ggplot(tab_df, aes(x = as.numeric(n), y = as.numeric(rejection), color = method)) +
  geom_line() +
  geom_point() +
  labs(x = "n", y = "Power", color = "Method") +
  ylim(0, 0.1) + #To see better what is happening
  theme_minimal()

#Saving-----------------------------------------------------------------------
# — save individually —
ggsave(
  filename = "plots/nNUniformPower.png",
  plot     = nNUniformPower,
  width    = 6, 
  height   = 4, 
  dpi      = 300
)

ggsave(
  filename = "plots/nNUniformRejection.png",
  plot     = nNUniformRejection,
  width    = 6, 
  height   = 4, 
  dpi      = 300
)


# turn them into a single grob
nNUniform_combined <- grid.arrange(
  nNUniformPower,
  nNUniformRejection,
  ncol = 2  # two columns side by side
)

# save the combined plot
ggsave(
  filename = "plots/nNUniform-combined.png",
  plot     = nNUniform_combined,
  width    = 12,   # adjust as you like
  height   = 6,    # adjust as you like
  dpi      = 300
)

################################################################################
#2. DIF when parameter I is varied---------------------------------------------

rm(list = ls())
#-------------------------------------------------------------------------------
#1.1 UNIFORM
#Data loading
load("results/variable-n-I-U.RData")
load("results/variable-I-res-U.RData")

#Extracting information from tab
tab_df <- as.data.frame(matrix(unlist(tab), nrow = 8,ncol=111, byrow = TRUE))
tab_df<- t(tab_df)
colnames(tab_df) <- c('N', 'n', 'I', 'type', 'gamma', 'dif-eff', 'power', 'rejection')
tab_df <- data.frame(tab_df)
tab_df <- tab_df[91:111,]
notimes <- nrow(tab_df) / 3
methods <- rep(c("Mantel","Log", "SIB"), notimes)
tab_df$method <- methods

small_tab <- tab_df[, c("I", "power", "rejection", 'method')]


#Plotting
IUniformPower <- ggplot(small_tab, aes(x = as.numeric(I), y = as.numeric(power), color = method)) +
  geom_line() +
  geom_point() +
  labs(x = "I", y = "Power", color = "Method") +
  theme_minimal()

IUniformRejection <- ggplot(tab_df, aes(x = as.numeric(I), y = as.numeric(rejection), color = method)) +
  geom_line() +
  geom_point() +
  labs(x = "I", y = "Power", color = "Method") +
  #ylim(0, 0.1) + #To see better what is happening
  theme_minimal()

#Saving-----------------------------------------------------------------------
# — save individually —
ggsave(
  filename = "plots/IUniformPower.png",
  plot     = IUniformPower,
  width    = 6, 
  height   = 4, 
  dpi      = 300
)

ggsave(
  filename = "plots/IUniformRejection.png",
  plot     = IUniformRejection,
  width    = 6, 
  height   = 4, 
  dpi      = 300
)


# turn them into a single grob
IUniform_combined <- grid.arrange(
  IUniformPower,
  IUniformRejection,
  ncol = 2  # two columns side by side
)

# save the combined plot
ggsave(
  filename = "plots/IUniform-combined.png",
  plot     = IUniform_combined,
  width    = 12,   # adjust as you like
  height   = 6,    # adjust as you like
  dpi      = 300
)

#-------------------------------------------------------------------------------
#Item wise inspection
#TODO - forgot to save
#Look how many items were detected etc


#################################################################################
rm(list = ls())
#2.2 NON-UNIFORM 
#Data loading
load("results/variable-n-I-NU-U.RData")
load("results/variable-I-res-NU.RData")

#Extracting information from tab
tab_df <- as.data.frame(matrix(unlist(tab), nrow = 8,ncol=153, byrow = TRUE))
tab_df<- t(tab_df)
colnames(tab_df) <- c('N', 'n', 'I', 'type', 'gamma', 'dif-eff', 'power', 'rejection')
tab_df <- data.frame(tab_df)
tab_df <- tab_df[112:153,]
notimes <- nrow(tab_df) / 6
methods <- rep(c("Mantel","MantelNormal","MantelLow", "MantelHigh", "Log", "CSIB"), notimes)
tab_df$method <- methods

small_tab <- tab_df[, c("I", "power", "rejection", 'method')]


#Plotting
INUniformPower <- ggplot(small_tab, aes(x = as.numeric(I), y = as.numeric(power), color = method)) +
  geom_line() +
  geom_point() +
  labs(x = "I", y = "Power", color = "Method") +
  theme_minimal()

INUniformRejection <- ggplot(tab_df, aes(x = as.numeric(I), y = as.numeric(rejection), color = method)) +
  geom_line() +
  geom_point() +
  labs(x = "I", y = "Power", color = "Method") +
  #ylim(0, 0.1) + #To see better what is happening
  theme_minimal()

#Saving-----------------------------------------------------------------------
# — save individually —
ggsave(
  filename = "plots/INUniformPower.png",
  plot     = INUniformPower,
  width    = 6, 
  height   = 4, 
  dpi      = 300
)

ggsave(
  filename = "plots/INUniformRejection.png",
  plot     = INUniformRejection,
  width    = 6, 
  height   = 4, 
  dpi      = 300
)


# turn them into a single grob
INUniform_combined <- grid.arrange(
  INUniformPower,
  INUniformRejection,
  ncol = 2  # two columns side by side
)

# save the combined plot
ggsave(
  filename = "plots/INUniform-combined.png",
  plot     = INUniform_combined,
  width    = 12,   # adjust as you like
  height   = 6,    # adjust as you like
  dpi      = 300
)

#Itemwise----------------------------------------------------------------------
# 1. get unique n’s (as characters, to index res_total)
unique_Is <- unique(tab_df$I) #TODO-generation again

# 2. for each n, run the function and grab only method + vote_ratio
results_extracted <- do.call(rbind, lapply(unique_Is, function(I_val) {
  # compute
  tmp <- calculate_itemwise_detection(res_total[[ as.character(I_val) ]])
  # select and add an 'n' column
  data.frame(
    item       = tmp$item,
    I          = I_val,
    method     = tmp$method,
    vote_ratio = tmp$vote_ratio,
    stringsAsFactors = FALSE
  )
}))

# 3. inspect
head(results_extracted)

#same for one item
ggplot(results_extracted, aes(x = as.numeric(n), y = as.numeric(vote_ratio), color = method)) +
  geom_line() +
  geom_point() +
  labs(x = "n", y = "Vote Ratio", color = "Method") +
  theme_minimal()


calculate_itemwise_detection(res_total$`70`)


################################################################################¨
#3. PENFIELD RATIOS
rm(list = ls())
#3.1 UNIFORM-------------------------------------------------------------------
load("results/penfieldRatio.RData")
load("results/penfieldRatio-resU.RData")

str(tab2)

#Extracting information from tab
tab_df <- as.data.frame(matrix(unlist(tab2), nrow = 8,ncol=12, byrow = TRUE))
tab_df<- t(tab_df)
colnames(tab_df) <- c('N', 'n', 'I', 'type', 'gamma', 'dif-eff', 'power', 'rejection')
tab_df <- data.frame(tab_df)
#tab_df <- tab_df[112:153,]
notimes <- nrow(tab_df) / 3
methods <- rep(c("Mantel", "Log", "SIB"), notimes)
tab_df$method <- methods

small_tab <- tab_df[, c("N", "power", "rejection", 'method')]


#Plotting
PUniformPower <- ggplot(small_tab, aes(x = as.numeric(N), y = as.numeric(power), color = method)) +
  geom_line() +
  geom_point() +
  labs(x = "n (ratio)", y = "Power", color = "Method") +
  theme_minimal()

PUniformRejection <- ggplot(tab_df, aes(x = as.numeric(N), y = as.numeric(rejection), color = method)) +
  geom_line() +
  geom_point() +
  labs(x = "n (ratio)", y = "Power", color = "Method") +
  #ylim(0, 0.1) + #To see better what is happening
  theme_minimal()

#Saving-----------------------------------------------------------------------
# — save individually —
ggsave(
  filename = "plots/PUniformPower.png",
  plot     = PUniformPower,
  width    = 6, 
  height   = 4, 
  dpi      = 300
)

ggsave(
  filename = "plots/PUniformRejection.png",
  plot     = PUniformRejection,
  width    = 6, 
  height   = 4, 
  dpi      = 300
)


# turn them into a single grob
PUniform_combined <- grid.arrange(
  PUniformPower,
  PUniformRejection,
  ncol = 2  # two columns side by side
)

# save the combined plot
ggsave(
  filename = "plots/PUniform-combined.png",
  plot     = PUniform_combined,
  width    = 12,   # adjust as you like
  height   = 6,    # adjust as you like
  dpi      = 300
)


#3.2 NON-UNIFORM-------------------------------------------------------------------
rm(list = ls())
load("results/penfieldRatioNU-tab.RData")
load("results/penfieldRatio-resNU.RData")

str(tab2)
#Extracting information from tab
tab_df <- as.data.frame(matrix(unlist(tab2), nrow = 8,ncol=24, byrow = TRUE))
tab_df<- t(tab_df)
colnames(tab_df) <- c('N', 'n', 'I', 'type', 'gamma', 'dif-eff', 'power', 'rejection')
tab_df <- data.frame(tab_df)
#tab_df <- tab_df[112:153,]
notimes <- nrow(tab_df) / 6
methods <- rep(c("Mantel", "MantelNormal","MantelLow","MantelHigh", "Log", "CSIB"), notimes)
tab_df$method <- methods

small_tab <- tab_df[, c("N", "power", "rejection", 'method')]


#Plotting
PNUniformPower <- ggplot(small_tab, aes(x = as.numeric(N), y = as.numeric(power), color = method)) +
  geom_line() +
  geom_point() +
  labs(x = "n (ratio)", y = "Power", color = "Method") +
  theme_minimal()

PNUniformRejection <- ggplot(tab_df, aes(x = as.numeric(N), y = as.numeric(rejection), color = method)) +
  geom_line() +
  geom_point() +
  labs(x = "n (ratio)", y = "Power", color = "Method") +
  #ylim(0, 0.1) + #To see better what is happening
  theme_minimal()

#Saving-----------------------------------------------------------------------
# — save individually —
ggsave(
  filename = "plots/PNUniformPower.png",
  plot     = PNUniformPower,
  width    = 6, 
  height   = 4, 
  dpi      = 300
)

ggsave(
  filename = "plots/PNUniformRejection.png",
  plot     = PNUniformRejection,
  width    = 6, 
  height   = 4, 
  dpi      = 300
)


# turn them into a single grob
PNUniform_combined <- grid.arrange(
  PNUniformPower,
  PNUniformRejection,
  ncol = 2  # two columns side by side
)

# save the combined plot
ggsave(
  filename = "plots/PNUniform-combined.png",
  plot     = PNUniform_combined,
  width    = 12,   # adjust as you like
  height   = 6,    # adjust as you like
  dpi      = 300
)
