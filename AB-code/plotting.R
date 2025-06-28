#Libraries-----------------------------------------------------------------------
rm(list = ls())
library(ggplot2)
library(dplyr)
library(tidyr)
library(gridExtra)
library(scales)  # for colour values
################################################################################
################################################################################
#IRT models showcase
okabe_ito <- c("#0072B2", "#E69F00", "#009E73", "#D55E00")
cb_palette <- c("Reference" = "#56B4E9", "Focal" = "#E69F00")

# Define theta range
theta <- seq(-4, 4, length.out = 300)

# Define logistic function and ICC
logit <- function(x) 1 / (1 + exp(-x))
icc <- function(theta, a, b, c, d) c + (d - c) * logit(a * (theta - b))


#######################################################
#4PL model
# Define 3 sets of item parameters: a, b, c, d
params <- data.frame(
  item = factor(c("Item 1", "Item 2", "Item 3")),
  a = c(0.8, 1.2, 1.5),
  b = c(-1, 0, 1),
  c = c(0.1, 0.2, 0.25),
  d = c(0.9, 0.95, 1)
)

# Generate ICC data
icc_data <- params %>%
  rowwise() %>%
  do({
    a <- .$a
    b <- .$b
    c <- .$c
    d <- .$d
    item <- .$item
    tibble(
      theta = theta,
      P = icc(theta, a, b, c, d),
      item = item
    )
  }) %>%
  bind_rows()



# Plot with colorblind-friendly palette
theoryIRT4PL <- ggplot(icc_data, aes(x = theta, y = P, color = item)) +
  geom_line(size = 1.2) +
  scale_color_manual(values = okabe_ito[1:3]) +
  labs(
    x = expression(paste("Latent Ability ", theta)),
    y = "Probability of Correct Response",
    color = "Item"
  ) +
  theme_minimal() +
  theme(text = element_text(size = 14)) +
  legend_inside

ggsave(
  filename = "plots/theoryIRT4PL.png",
  plot     = theoryIRT4PL,
  width    = 6, 
  height   = 4, 
  dpi      = 300
)

###############################################################
###############################################################
# Common parameters
a0 <- 1
b0 <- 0
c <- 0.2

# --------------------
# UNIFORM DIF (shift in difficulty b)
# --------------------
b1 <- b0 + 1.1  # shift difficulty
uniform_data <- data.frame(
  theta = rep(theta, 2),
  prob = c(icc(theta, a0, b0, c,1), icc(theta, a0, b1, c, 1)),
  group = rep(c("Reference", "Focal"), each = length(theta)),
  type = "Uniform DIF"
)

# --------------------
# NON-UNIFORM DIF (change in discrimination a)
# --------------------
delta1 <- 1
a1 <- 2 * a0/ (2 + delta1 * a0 / ((1 - c) * log(2)))  # change discrimination
nonuniform_data <- data.frame(
  theta = rep(theta, 2),
  prob = c(icc(theta, a0, b0, c, 1), icc(theta, a1, b0, c, 1)),
  group = rep(c("Reference", "Focal"), each = length(theta)),
  type = "Non-Uniform DIF"
)

# --------------------
# PLOTS
# --------------------
# Colorblind-friendly palette

legend_inside <- theme(
  legend.position = c(0.99,0.01),
  legend.position.inside = "bottom-right",
  legend.justification = c("right", "bottom"),
  legend.background = element_rect(fill = alpha("white", 0.8), color = NA),
  legend.title = element_text(size = 10),
  legend.text = element_text(size = 9)  # ← smaller legend labels
)


intro_unif <- ggplot(uniform_data, aes(x = theta, y = prob, color = group)) +
  geom_line(size = 1.2) +
  scale_color_manual(values = cb_palette) +
  labs(
    x = expression(paste("Latent Ability ", theta)),
    y = "Probability of Correct Response"
  ) +
  theme_minimal(base_size = 14) +
  legend_inside

intro_nonunif <- ggplot(nonuniform_data, aes(x = theta, y = prob, color = group)) +
  geom_line(size = 1.2) +
  scale_color_manual(values = cb_palette) +
  labs(
    x = expression(paste("Latent Ability ", theta)),
    y = "Probability of Correct Response"
  ) +
  theme_minimal(base_size = 14) +
  legend_inside


# Show plots
theoryIntro <- grid.arrange(
  intro_unif,
  intro_nonunif,
  ncol = 2  # two columns side by side
)


ggsave(
  filename = "plots/theoryIntro.png",
  plot     = theoryIntro,
  width    = 6, 
  height   = 4, 
  dpi      = 300
)

############################################################
a0 <- 1
b0 <- 0
c <- 0.2


# --- UNIFORM DIF with different difficulty shifts ---
b_shifts <- c(0.5, 1.0, 1.5)
uniform_data_multi <- lapply(b_shifts, function(shift) {
  data.frame(
    theta = theta,
    prob = icc(theta, a0, b0 + shift, c, 1),
    group = paste0("Focal (δ =" ,shift,")"),
    type = "Focal"
  )
}) %>% bind_rows()

ref_uniform <- data.frame(
  theta = theta,
  prob = icc(theta, a0, b0, c, 1),
  group = "Reference",
  type = "Reference"
)

plot_data_uniform <- bind_rows(ref_uniform, uniform_data_multi)

# --- NON-UNIFORM DIF with different discrimination changes ---
deltas <- c(0.4, 0.8, 1.2)
nonuniform_data_multi <- lapply(deltas, function(delta) {
  a1 <- 2 * a0 / (2 + delta * a0 / ((1 - c) * log(2)))
  data.frame(
    theta = theta,
    prob = icc(theta, a1, b0, c, 1),
    group = paste0("Focal (Δ = ", delta, ")"),
    type = "Focal"
  )
}) %>% bind_rows()

ref_nonuniform <- data.frame(
  theta = theta,
  prob = icc(theta, a0, b0, c, 1),
  group = "Reference",
  type = "Reference"
)

plot_data_nonuniform <- bind_rows(ref_nonuniform, nonuniform_data_multi)

plot_data_uniform <- bind_rows(ref_uniform, uniform_data_multi) %>%
  mutate(linetype = ifelse(group == "Reference", "Reference", "Focal"))
plot_data_nonuniform <- bind_rows(ref_uniform, nonuniform_data_multi) %>%
  mutate(linetype = ifelse(group == "Reference", "Reference", "Focal"))

# --- Plot: Uniform DIF Magnitude ---
theoryIRTUniformUnifb<- ggplot(plot_data_uniform, aes(x = theta, y = prob, color = group, linetype = linetype)) +
  geom_line(size = 1.2) +
  scale_color_manual(values = okabe_ito ) +
  scale_linetype_manual(
    values = c(
      "Reference" = "dashed",
      "Focal"     = "solid"
    ),guide = "none") +
  labs(
    x = expression(paste("Latent Ability ", theta)),
    y = "Probability of Correct Response"
  ) +
  theme_minimal(base_size = 14)+
  legend_inside

ggsave(
  filename = "plots/theoryIRTUniformUnifb.png",
  plot     = theoryIRTUniformUnifb,
  width    = 6, 
  height   = 4, 
  dpi      = 300
)


# --- Plot: Non-uniform DIF Magnitude ---
theoryIRTNONUniformUnifa <- ggplot(plot_data_nonuniform, aes(x = theta, y = prob, color = group,linetype = linetype)) +
  geom_line(size = 1.2) +
  scale_color_manual(values = okabe_ito) +
  scale_linetype_manual(
    values = c(
      "Reference" = "dashed",
      "Focal"     = "solid"
    ),guide = "none") +
  labs(
    x = expression(paste("Latent Ability ", theta)),
    y = "Probability of Correct Response"
  ) +
  theme_minimal(base_size = 14) +
  legend_inside

ggsave(
  filename = "plots/theoryIRTNONUniformUnifa.png",
  plot     = theoryIRTNONUniformUnifa,
  width    = 6, 
  height   = 4, 
  dpi      = 300
)



################################################################################
###############################################################################¨
#B.Simulations plots

#1. DIF when parameter n is varied---------------------------------------------
#1.1 UNIFORM

#Data loading
load("results/variable-n-res-U.RData")
load("results/variable-n-U.RData")
res_total$`10`$metadata
tab
#Extracting information from tab
tab_df <- as.data.frame(matrix(unlist(tab), nrow = 8,ncol=45, byrow = TRUE))
tab_df<- t(tab_df)
notimes <- nrow(tab_df) / 3
methods <- rep(c("Mantel", "Log", "SIB"), notimes)
colnames(tab_df) <- c('N', 'n', 'I', 'type', 'gamma', 'dif-eff', 'power', 'rejection')
tab_df <- data.frame(tab_df)
tab_df$method <- methods

small_tab <- tab_df[, c("n", "power", "rejection", 'method')]
tab_df[, c("n", "power", 'method')]

#Plotting
nUniformPower <- ggplot(small_tab, aes(x = as.numeric(n), y = as.numeric(power), color = method)) +
  geom_line() +
  geom_point() +
  labs(x = "n", y = "Power Rate", color = "Method") +
  theme_minimal()

nUniformRejection <- ggplot(tab_df, aes(x = as.numeric(n), y = as.numeric(rejection), color = method)) +
  geom_line() +
  geom_point() +
  labs(x = "n", y = "Rejection Rate", color = "Method") +
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
load("results/variable-n-res-NONU.RData")
load("results/variable-n-NONU.RData")


#Extracting information from tab
tab_df <- as.data.frame(matrix(unlist(tab), nrow = 8,ncol=90, byrow = TRUE))
tab_df<- t(tab_df)
notimes <- nrow(tab_df) / 6
methods <- rep(c("Mantel","MantelNormal","MantelLow", "MantelHigh", "Log", "CSIB"), notimes)
colnames(tab_df) <- c('N', 'n', 'I', 'type', 'gamma', 'dif-eff', 'power', 'rejection')
tab_df <- data.frame(tab_df)
tab_df$method <- methods

small_tab <- tab_df[, c("n", "power", "rejection", 'method')]
tab_df[, c("n", "power", 'method')]

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
