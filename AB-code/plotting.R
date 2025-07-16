#Libraries-----------------------------------------------------------------------
rm(list = ls())
library(ggplot2)
library(dplyr)
library(tidyr)
library(gridExtra)
library(scales)  # for colour values
################################################################################
################################################################################
################################################################################
#PART 1: Theory-----------------------------------------------------------------

# Colors
okabe_ito <- c("#0072B2", "#E69F00", "#009E73", "#D55E00")
cb_palette <- c("Reference" = "#0072B2", "Focal" = "#E69F00")

# Define theta range
theta <- seq(-6, 6, length.out = 300)

# Define logistic function and ICC
logit <- function(x) 1 / (1 + exp(-x))
icc <- function(theta, a, b, c, d) c + (d - c) * logit(a * (theta - b))

legend_inside <- theme(
  legend.position = c(0.99,0.01),
  legend.position.inside = "bottom-right",
  legend.justification = c("right", "bottom"),
  legend.background = element_rect(fill = alpha("white", 0.8), color = NA),
  legend.title = element_text(size = 10),
  legend.text = element_text(size = 9)  # ← smaller legend labels
)

################################################################################
#4PL model (Figure 2)
# Define 3 sets of item parameters: a, b, c, d
params <- data.frame(
  item = factor(c("4PL", "3PL", "2PL", "1PL")),
  a = c(0.8, 0.8, 0.8, 1),
  b = c(-2, -1, 0,1),
  c = c(0.2, 0.2, 0,0),
  d = c(0.95, 1, 1, 1)
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



#4PL model (Figure 2)
theoryIRT4PL <- ggplot(icc_data, aes(x = theta, y = P, color = item)) +
  geom_line(size = 1.2) +
  scale_color_manual(values = okabe_ito[1:4], name= NULL) +
  labs(
    x = expression(paste("Latent Ability ", theta)),
    y = "Probability of Correct Response",
    color = "Item"
  ) +
  theme_minimal() +
  ylim(0,1)+
  theme(text = element_text(size = 14)) +
  legend_inside

ggsave(
  filename = "plots/theoryIRT4PL.png",
  plot     = theoryIRT4PL,
  width    = 6, 
  height   = 4, 
  dpi      = 300
)

################################################################################
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
intro_unif <- ggplot(uniform_data, aes(x = theta, y = prob, color = group)) +
  geom_line(size = 1.2) +
  scale_color_manual(values = cb_palette,name = NULL) +
  labs(
    x = expression(paste("Latent Ability ", theta)),
    y = "Probability of Correct Response"
  ) +
  ylim(0,1) + 
  theme_minimal(base_size = 14) +
  legend_inside

intro_nonunif <- ggplot(nonuniform_data, aes(x = theta, y = prob, color = group)) +
  geom_line(size = 1.2) +
  scale_color_manual(values = cb_palette, name = NULL) +
  labs(
    x = expression(paste("Latent Ability ", theta)),
    y = "Probability of Correct Response"
  ) +
  ylim(0,1)+
  theme_minimal(base_size = 14) +
  legend_inside


# Figure 1: unif and non unif DIF
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

################################################################################
# How is DIF generated using IRT

#_-----------------------------------------------------
# --- UNIFORM DIF with different difficulty shifts
#-------------------------------------------------------

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
plot_data_uniform <- bind_rows(ref_uniform, uniform_data_multi) %>%
  mutate(linetype = ifelse(group == "Reference", "Reference", "Focal"))
#-------------------------------------------------------------
# --- NON-UNIFORM DIF with different discrimination changes ---
#--------------------------------------------------------------
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
plot_data_nonuniform <- bind_rows(ref_uniform, nonuniform_data_multi) %>%
  mutate(linetype = ifelse(group == "Reference", "Reference", "Focal"))

#-------------------------------------------------------------------------------
# --- Plot: Uniform DIF Magnitude --- (Figure 3)
theoryIRTUniformUnifb <- ggplot(plot_data_uniform, aes(x = theta, y = prob, color = group, linetype = type)) +
  geom_line(size = 1.2) +
  scale_color_manual(values = okabe_ito, name = NULL) +
  scale_linetype_manual(
    values = c(
      "Reference" = "dashed",
      "Focal" = "solid"
    ), name = NULL) +
  guides(
    color = guide_legend(
      override.aes = list(
        linetype = c("solid", "solid", "solid", "dashed")
      )
    ),
    linetype = "none"
  ) +
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

#---------------------------------------------------------------------------------
# --- Plot: Non-uniform DIF Magnitude --- (Figure 4)
theoryIRTNONUniformUnifa <- ggplot(plot_data_nonuniform, aes(x = theta, y = prob, color = group, linetype = linetype)) +
  geom_line(size = 1.2) +
  scale_color_manual(values = okabe_ito, name = NULL) +
  scale_linetype_manual(
    values = c(
      "Reference" = "dashed",
      "Focal"     = "solid"
    ), name = NULL) +
  guides(
    color = guide_legend(
      override.aes = list(
        linetype = c("solid", "solid", "solid", "dashed")
      )
    ),
    linetype = "none"
  ) +
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



theoryParams <- grid.arrange(
  theoryIRTUniformUnifb,
  theoryIRTNONUniformUnifa,
  ncol = 2  # two columns side by side
)


ggsave(
  filename = "plots/theoryParams.png",
  plot     = theoryParams,
  width    = 6, 
  height   = 4, 
  dpi      = 300
)

###############################################################################
################################################################################
###SIMULATION -----------------------------------------------------------

legend_inside <- theme(
  legend.position = c(0.99,0.01),
  legend.position.inside = "bottom-right",
  legend.justification = c("right", "bottom"),
  legend.background = element_rect(fill = alpha("white", 0.8), color = NA),
  legend.title = element_text(size = 10),
  legend.text = element_text(size = 9)  # ← smaller legend labels
)

#Scenario 1: n var
#Scenario 2: I var (gamma stays) (Item wise analysis)
#Scenario 3: I var, gamma var
#Scenario 4: ratio rho var
#Scenario 5: DIF size var 

################################################################################
#SCENARIO 1: DIF when parameter N is varied-------------------------------------
#1.1 UNIFORM------------------------------------------------------------------
rm(list = ls())
#Data loading
load("results/variable-n-res-U.RData")
load("results/variable-n-U.RData")
okabe_ito <- c("#0072B2", "#E69F00", "#009E73", "#D55E00", "#CC79A7", "#F0E442")
#Extracting information from tab
tab_df <- as.data.frame(matrix(unlist(tab), nrow = 8,ncol=45, byrow = TRUE))
tab_df<- t(tab_df)
notimes <- nrow(tab_df) / 3
methods <- rep(c("MH", "LR", "SIB"), notimes)
colnames(tab_df) <- c('M', 'N', 'I', 'type', 'gamma', 'difSize', 'power', 'rejection')
tab_df <- data.frame(tab_df)
tab_df$method <- methods

small_tab <- tab_df[, c("N", "power", "rejection", 'method')]
tab_df[, c("N", "power", 'method')]





#Plotting#PlottingMantelStat
legend_inside <- theme(
  axis.text = element_text(size = 12),    # bigger x/y tick labels
  axis.title = element_text(size = 16), 
  axis.text.x = element_text(angle = 45, hjust = 1),
  # bigger axis titles
  #legend.title = element_text(size = 14),
  #legend.text = element_text(size = 13),
  legend.position = c(0.99,0.01),
  legend.position.inside = "bottom-right",
  legend.justification = c("right", "bottom"),
  legend.background = element_rect(fill = alpha("white", 0.8), color = NA),
  legend.title = element_text(size = 10),
  legend.text = element_text(size = 9)  # ← smaller legend labels
)


nUniformPower <- ggplot(small_tab, aes(x = as.numeric(N), y = as.numeric(power), color = method)) +
  geom_line() +
  geom_point() +
  labs(x = "Number of respondents N", y = "Power Rate w", color = "Method") +
  theme_minimal(base_size = 14) +
  scale_x_continuous(
    breaks = c(0, 100, 200, 400,500,1000, 1200, 2000),
    labels = c("0","100","200", "400","500", "1000", "1200", "2000")
  )+
  ylim(0, 1) +
  scale_color_manual(values = okabe_ito[1:3], name= NULL) +
  legend_inside

nUniformRejection <- ggplot(small_tab, aes(x = as.numeric(N), y = as.numeric(rejection), color = method)) +
  geom_line() +
  geom_point() +
  labs(x = "Number of respondents N", y = "Rejection Rate r", color = "Method") +
  ylim(0, 0.1) + #To see better what is happening
  #scale_x_continuous(breaks = unique(as.numeric(small_tab$N)))+
  theme_minimal(base_size = 14) +
  scale_x_continuous(
    breaks = c(0, 100, 200, 400,500,1000, 1200, 2000),
    labels = c("0","100","200", "400","500", "1000", "1200", "2000")
  )+
  scale_color_manual(values = okabe_ito[1:3], name= NULL) +
  legend_inside


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
  nUniformRejection,
  nUniformPower,
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

##############################################################################xx
nUniformRejectionSmall <- ggplot(small_tab, aes(x = as.numeric(N), y = as.numeric(rejection), color = method)) +
  geom_line() +
  geom_point() +
  labs(x = "Number of respondents N", y = "Rejection Rate r", color = "Method") +
  xlim(10, 100) + #To see better what is happening
  theme_minimal(base_size = 14) +
  ylim(0, 1) +
  scale_x_continuous(
    breaks = c(0, 100, 200, 400,500,1000, 1200, 2000),
    labels = c("0","100","200", "400","500", "1000", "1200", "2000")
  )+
  scale_color_manual(values = okabe_ito[1:3], name= NULL) +
  legend_inside

nUniformPowerSmall <- ggplot(small_tab, aes(x = as.numeric(N), y = as.numeric(power), color = method)) +
  geom_line() +
  geom_point() +
  labs(x = "Number of respondents N", y = "Power Rate w", color = "Method") +
  xlim(10,100) +
  scale_x_continuous(
    breaks = c(0, 100, 200, 400,500,1000, 1200, 2000),
    labels = c("0","100","200", "400","500", "1000", "1200", "2000")
  )+
  theme_minimal(base_size = 14) +
  scale_color_manual(values = okabe_ito[1:3], name= NULL) +
  legend_inside

nUniformSmall_combined <- grid.arrange(
  nUniformPowerSmall,
  nUniformRejectionSmall,
  ncol = 2  # two columns side by side
)
ggsave(
  filename = "plots/nUniform-combinedSmall.png",
  plot     = nUniformSmall_combined,
  width    = 12,   # adjust as you like
  height   = 6,    # adjust as you like
  dpi      = 300
)


#1.2 NON-UNIFORM----------------------------------------------------------------
rm(list = ls())
#Data loading
load("results/variable-n-res-NONU.RData")
load("results/variable-n-NONU.RData")



#Extracting information from tab
tab_df <- as.data.frame(matrix(unlist(tab), nrow = 8,ncol=90, byrow = TRUE))
tab_df<- t(tab_df)
notimes <- nrow(tab_df) / 6
methods <- rep(c("MH","MH-NU","MH-L", "MH-H", "LR", "CSIB"), notimes)
colnames(tab_df) <- c('M', 'N', 'I', 'type', 'gamma', 'difSize', 'power', 'rejection')
tab_df <- data.frame(tab_df)
tab_df$method <- methods

small_tab <- tab_df[, c("N", "power", "rejection", 'method')]
small_tab <- small_tab[!small_tab$method == "MH", ]

#Plottin
legend_inside <- theme(
  axis.text = element_text(size = 12),    # bigger x/y tick labels
  axis.title = element_text(size = 16), 
  axis.text.x = element_text(angle = 45, hjust = 1),
  # bigger axis titles
  #legend.title = element_text(size = 14),
  #legend.text = element_text(size = 13),
  legend.position = c(0.99,0.01),
  legend.position.inside = "bottom-right",
  legend.justification = c("right", "bottom"),
  legend.background = element_rect(fill = alpha("white", 0.8), color = NA),
  legend.title = element_text(size = 10),
  legend.text = element_text(size = 9)  # ← smaller legend labels
)
 okabe_ito <- c("#0072B2", "#E69F00", "#009E73", "#D55E00", "#CC79A7", "#F0E442")


nNUniformPower <- ggplot(small_tab, aes(x = as.numeric(N), y = as.numeric(power), color = method)) +
  geom_line() +
  geom_point() +
  labs(x = "Number of respondents N", y = "Power Rate w", color = "Method") +
  theme_minimal(base_size = 14) +
  ylim(0, 1) +
  scale_x_continuous(
    breaks = c(0, 100, 200, 400,500,1000, 1200, 2000),
    labels = c("0","100","200", "400","500", "1000", "1200", "2000")
  )+
  scale_color_manual(values = okabe_ito[1:5], name= NULL) +
  legend_inside

nNUniformRejection <- ggplot(small_tab, aes(x = as.numeric(N), y = as.numeric(rejection), color = method)) +
  geom_line() +
  geom_point() +
  labs(x = "Number of respondents N", y = "Rejection Rate r", color = "Method") +
  ylim(0, 0.1) + #To see better what is happening
  theme_minimal(base_size = 14) +
  scale_x_continuous(
    breaks = c(0, 100, 200, 400,500,1000, 1200, 2000),
    labels = c("0","100","200", "400","500", "1000", "1200", "2000")
  )+
  scale_color_manual(values = okabe_ito[1:5], name= NULL) +
  legend_inside


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
  nNUniformRejection,
  nNUniformPower,
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



# turn them into a single grob
nRejection_combined <- grid.arrange(
  nUniformRejection,
  nNUniformRejection,
  ncol = 2  # two columns side by side
)

# turn them into a single grob
nPower_combined <- grid.arrange(
  nUniformPower,
  nNUniformPower,
  ncol = 2  # two columns side by side
)


# save the combined plot
ggsave(
  filename = "plots/nRejection-combined.png",
  plot     = nRejection_combined,
  width    = 12,   # adjust as you like
  height   = 6,    # adjust as you like
  dpi      = 300
)



# save the combined plot
ggsave(
  filename = "plots/nPower-combined.png",
  plot     = nPower_combined,
  width    = 12,   # adjust as you like
  height   = 6,    # adjust as you like
  dpi      = 300
)


###############################################################################
legend_inside <- theme(
  legend.position = c(0.99,0.60),
  legend.position.inside = "bottom-right",
  legend.justification = c("right", "bottom"),
  legend.background = element_rect(fill = alpha("white", 0.8), color = NA),
  legend.title = element_text(size = 10),
  legend.text = element_text(size = 9)  # ← smaller legend labels
)

nUniformRejectionSmall <- ggplot(small_tab, aes(x = as.numeric(N), y = as.numeric(rejection), color = method)) +
  geom_line() +
  geom_point() +
  labs(x = "Number of respondents N", y = "Rejection Rate r", color = "Method") +
  xlim(10, 100) + #To see better what is happening
  theme_minimal(base_size = 14) +
  scale_color_manual(values = okabe_ito[1:6], name= NULL) +
  legend_inside

nUniformPowerSmall <- ggplot(small_tab, aes(x = as.numeric(N), y = as.numeric(power), color = method)) +
  geom_line() +
  geom_point() +
  labs(x = "Number of respondents N", y = "Power Rate w", color = "Method") +
  xlim(10,100) +
  theme_minimal(base_size = 14) +
  scale_color_manual(values = okabe_ito[1:6], name= NULL) +
  legend_inside

nUniformSmall_combined <- grid.arrange(
  nUniformPowerSmall,
  nUniformRejectionSmall,
  ncol = 2  # two columns side by side
)
ggsave(
  filename = "plots/nUniform-combinedSmall.png",
  plot     = nUniformSmall_combined,
  width    = 12,   # adjust as you like
  height   = 6,    # adjust as you like
  dpi      = 300
)



################################################################################
# SCENARIO 2: DIF when parameter I is varied-------------------------------------

rm(list = ls())
# 2.1 UNIFORM-------------------------------------------------------------------
#Data loading
load("results/variable-I-U.RData")
load("results/variable-I-res-U.RData")

#Extracting information from tab
tab_df <- as.data.frame(matrix(unlist(tab), nrow = 8,ncol=24, byrow = TRUE))
tab_df<- t(tab_df)
colnames(tab_df) <- c('M', 'N', 'I', 'type', 'gamma', 'difSize', 'power', 'rejection')
tab_df <- data.frame(tab_df)
#tab_df <- tab_df[91:192,]
notimes <- nrow(tab_df) / 3
methods <- rep(c("MH", "LR", "SIB"), notimes)
tab_df$method <- methods
tab_df <- tab_df %>%
  mutate(DIFItems = as.numeric(I) * as.numeric(gamma))
small_tab <- tab_df[, c("I","DIFItems", "power", "rejection", 'method')]
small_tab <- small_tab %>%
  mutate(xlabel = paste0(I, " (", DIFItems, ")")) %>%
  arrange(as.numeric(I))  # make sure it's ordered

tab_df[, c("I", "power", 'method')]
#Plotting
legend_inside <- theme(
  axis.text = element_text(size = 12),    # bigger x/y tick labels
  axis.title = element_text(size = 16),
  # bigger axis titles
  #legend.title = element_text(size = 14),
  #legend.text = element_text(size = 13),
  legend.position = c(0.99,0.01),
  legend.position.inside = "bottom-right",
  legend.justification = c("right", "bottom"),
  legend.background = element_rect(fill = alpha("white", 0.8), color = NA),
  legend.title = element_text(size = 10),
  legend.text = element_text(size = 9))  # ← smaller legend labels
#Plotting
 okabe_ito <- c("#0072B2", "#E69F00", "#009E73", "#D55E00", "#CC79A7", "#F0E442")




# Create unique x-axis labels in the correct order
xlabel_levels <- unique(paste0(small_tab$I, " (", small_tab$DIFItems, ")"))

IUniformPower <- ggplot(
  small_tab,
  aes(
    x = factor(xlabel, levels = xlabel_levels),
    y = as.numeric(power),
    color = method
  )
) +
  geom_line(aes(group = method)) +
  geom_point() +
  labs(x = "Number of Items I (Number of DIF items γ*I)", y = "Power Rate w", color = "Method") +
  theme_minimal(base_size = 14) +
  scale_color_manual(values = okabe_ito[1:3], name= NULL) +
  legend_inside


IUniformRejection <- ggplot(
  small_tab,
  aes(
    x = factor(xlabel, levels = xlabel_levels),
    y = as.numeric(rejection),
    color = method
  )
) +
  geom_line(aes(group = method)) +
  geom_point() +
  labs(x = "Number of Items I (Number of DIF items γ*I)", y = "Rejection Rate r", color = "Method") +
  theme_minimal(base_size = 14) +
  scale_color_manual(values = okabe_ito[1:3], name= NULL) +
  legend_inside

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
# 1. get unique n’s (as characters, to index res_total)
unique_ns <- unique(tab_df$I)

# 2. for each n, run the function and grab only method + vote_ratio
results_extracted <- do.call(rbind, lapply(unique_ns, function(n_val) {
  # compute
  tmp <- calculate_itemwise_detection(res_total[[ as.character(n_val) ]])
  # select and add an 'n' column
  data.frame(
    I          = n_val,
    method     = tmp$method,
    vote_ratio = tmp$vote_ratio,
    stringsAsFactors = FALSE
  )
}))

# 3. inspect
head(results_extracted)

#same for one item
ggplot(results_extracted, aes(x = as.numeric(I), y = as.numeric(vote_ratio), color = method)) +
  geom_line() +
  geom_point() +
  labs(x = "Number of respondents N", y = "Vote Ratio", color = "Method") +
  theme_minimal(base_size = 14) +
  legend_inside

#TODO - forgot to save
#Look how many items were detected etc


#2.2 NON-UNIFORM---------------------------------------------------------------
#Data loading
rm(list = ls())
load("results/variable-I-NONU.RData")
load("results/variable-I-res-NONU.RData")

#Extracting information from tab
tab_df <- as.data.frame(matrix(unlist(tab), nrow = 8,ncol=48, byrow = TRUE))
tab_df<- t(tab_df)
colnames(tab_df) <- c('M', 'N', 'I', 'type', 'gamma', 'difSize', 'power', 'rejection')
tab_df <- data.frame(tab_df)
#tab_df <- tab_df[112:153,]
notimes <- nrow(tab_df) / 6
methods <- rep(c("MH","MH-NU","MH-L", "MH-H", "LR", "CSIB"), notimes)
tab_df$method <- methods
tab_df <- tab_df %>%
  mutate(DIFItems = as.numeric(I) * as.numeric(gamma))
small_tab <- tab_df[, c("I","DIFItems", "power", "rejection", 'method')]
small_tab <- small_tab %>%
  mutate(xlabel = paste0(I, " (", DIFItems, ")")) %>%
  arrange(as.numeric(I))
xlabel_levels <- unique(paste0(small_tab$I, " (", small_tab$DIFItems, ")"))

tab_df[, c("I", "power", 'method')]

#Plotting
legend_inside <- theme(
  legend.position = c(0.99,0.01),
  legend.position.inside = "bottom-right",
  legend.justification = c("right", "bottom"),
  legend.background = element_rect(fill = alpha("white", 0.8), color = NA),
  legend.title = element_text(size = 10),
  legend.text = element_text(size = 9)  # ← smaller legend labels
)

 okabe_ito <- c("#0072B2", "#E69F00", "#009E73", "#D55E00", "#CC79A7", "#F0E442")

INUniformPower <- ggplot(
  small_tab,
  aes(
    x = factor(xlabel, levels = xlabel_levels),
    y = as.numeric(power),
    color = method
  )
) +
  geom_line(aes(group = method)) +
  geom_point() +
  labs(x = "Number of Items I (Number of DIF items γ*I)", y = "Power Rate w", color = "Method") +
  theme_minimal(base_size = 14) +
  scale_color_manual(values = okabe_ito[1:6], name= NULL) +
  legend_inside


INUniformRejection <- ggplot(
  small_tab,
  aes(
    x = factor(xlabel, levels = xlabel_levels),
    y = as.numeric(rejection),
    color = method
  )
) +
  geom_line(aes(group = method)) +
  geom_point() +
  labs(x = "Number of Items I (Number of DIF items γ*I)", y = "Rejection Rate r", color = "Method") +
  theme_minimal(base_size = 14) +
  scale_color_manual(values = okabe_ito[1:6], name= NULL) +
  legend_inside


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
ggplot(results_extracted, aes(x = as.numeric(N), y = as.numeric(vote_ratio), color = method)) +
  geom_line() +
  geom_point() +
  labs(x = "Number of respondents N", y = "Vote Ratio", color = "Method") +
  theme_minimal(base_size = 14) +
  legend_inside


calculate_itemwise_detection(res_total$`70`)
################################################################################
#SCENARIO 3: I and gamma are varied
#3.1 UNIFORM--------------------------------------------------------------------
#Data loading
rm(list = ls())
# Load the first dataset
load("results/variable-I-rat-U.RData")
load("results/variable-I-rat-res-U.RData")
tab1 <- as.data.frame(matrix(unlist(tab2), ncol = 24, byrow = TRUE))
tab1 <- t(tab1)
colnames(tab1) <- c('M', 'ratio', 'N', 'I', 'type', 'gamma', 'difSize', 'power', 'rejection')
tab1 <- data.frame(tab1)

# Load the second dataset
load("results/variable-I-rat-U-5-10.RData")
load("results/variable-I-rat-res-U-5-10.RData")
tab2 <- as.data.frame(matrix(unlist(tab2), ncol = 6, byrow = TRUE))
tab2 <- t(tab2)
colnames(tab2) <- c('M', 'ratio', 'N', 'I', 'type', 'gamma', 'difSize', 'power', 'rejection')
tab2 <- data.frame(tab2)

# Combine the two tables row-wise
tab_df <- rbind(tab1, tab2)

# Add the method column
notimes <- nrow(tab_df) / 3
methods <- rep(c("MH", "LR", "SIB"), notimes)
tab_df$method <- methods
#Extracting information from tab

small_tab <- tab_df[, c("I", "power", "rejection", 'method')]
tab_df[, c("I", "power", 'method')]

#Plotting
legend_inside <- theme(
  axis.text = element_text(size = 12),    # bigger x/y tick labels
  axis.title = element_text(size = 16),
  # bigger axis titles
  #legend.title = element_text(size = 14),
  #legend.text = element_text(size = 13),
  legend.position = c(0.99,0.01),
  legend.position.inside = "bottom-right",
  legend.justification = c("right", "bottom"),
  legend.background = element_rect(fill = alpha("white", 0.8), color = NA),
  legend.title = element_text(size = 10),
  legend.text = element_text(size = 9))  # ← smaller legend labels
 okabe_ito <- c("#0072B2", "#E69F00", "#009E73", "#D55E00", "#CC79A7", "#F0E442")


IRATUniformPower <- ggplot(small_tab, aes(x = as.numeric(I), y = as.numeric(power), color = method)) +
  geom_line() +
  geom_point() +
  labs(x = "Number of Items I", y = "Power Rate w", color = "Method") +
  theme_minimal(base_size = 14) +
  ylim(0, 1) +
  scale_x_continuous(breaks = unique(as.numeric(small_tab$I)))+
  scale_color_manual(values = okabe_ito[1:3], name= NULL) +
  legend_inside

IRATUniformRejection <- ggplot(tab_df, aes(x = as.numeric(I), y = as.numeric(rejection), color = method)) +
  geom_line() +
  geom_point() +
  labs(x = "Number of Items I", y = "Rejection Rate r", color = "Method") +
  ylim(0, 0.1) + #To see better what is happening
  scale_x_continuous(breaks = unique(as.numeric(small_tab$I)))+
  theme_minimal(base_size = 14) +
  scale_color_manual(values = okabe_ito[1:3], name= NULL) +
  legend_inside
tab2

# — save individually —
ggsave(
  filename = "plots/IRATUniformPower.png",
  plot     = IRATUniformPower,
  width    = 6, 
  height   = 4, 
  dpi      = 300
)

ggsave(
  filename = "plots/IRATUniformRejection.png",
  plot     = IRATUniformRejection,
  width    = 6, 
  height   = 4, 
  dpi      = 300
)


# turn them into a single grob
IRATUniform_combined <- grid.arrange(
  IRATUniformRejection,
  IRATUniformPower,
  ncol = 2  # two columns side by side
)

# save the combined plot
ggsave(
  filename = "plots/IRATUniform-combined.png",
  plot     = IRATUniform_combined,
  width    = 12,   # adjust as you like
  height   = 6,    # adjust as you like
  dpi      = 300
)


#3.2 NON-UNIFORM----------------------------------------------------------------
# Load the first dataset
rm(list = ls())
load("results/variable-I-rat-NONU.RData")
load("results/variable-I-rat-res-NONU.RData")
tab1 <- as.data.frame(matrix(unlist(tab2), ncol = 40, byrow = TRUE))
tab1 <- t(tab1)
colnames(tab1) <- c('M', 'ratio', 'N', 'I', 'type', 'gamma', 'difSize', 'power', 'rejection')
tab1 <- data.frame(tab1)

res_total2[[10]]$metadata
# Load the second dataset
load("results/variable-I-rat-NONU-5-10.RData")
load("results/variable-I-rat-res-NONU-5-10.RData")
tab2 <- as.data.frame(matrix(unlist(tab2), ncol = 10, byrow = TRUE))
tab2 <- t(tab2)
colnames(tab2) <- c('M', 'ratio', 'N', 'I', 'type', 'gamma', 'difSize', 'power', 'rejection')
tab2 <- data.frame(tab2)

# Combine the two tables row-wise
tab_df <- rbind(tab1, tab2)

# Add method labels (5 methods this time)
notimes <- nrow(tab_df) / 5
methods <- rep(c("MH-NU", "MH-L", "MH-H", "LR", "CSIB"), notimes)
tab_df$method <- methods
small_tab <- tab_df[, c("I", "power", "rejection", 'method')]
#small_tab <- small_tab[!small_tab$method == "MH", ]
tab_df[, c("I", "power", 'method')]

legend_inside <- theme(
  axis.text = element_text(size = 12),    # bigger x/y tick labels
  axis.title = element_text(size = 16),
  # bigger axis titles
  #legend.title = element_text(size = 14),
  #legend.text = element_text(size = 13),
  legend.position = c(0.99,0.01),
  legend.position.inside = "bottom-right",
  legend.justification = c("right", "bottom"),
  legend.background = element_rect(fill = alpha("white", 0.8), color = NA),
  legend.title = element_text(size = 10),
  legend.text = element_text(size = 9))

 okabe_ito <- c("#0072B2", "#E69F00", "#009E73", "#D55E00", "#CC79A7", "#F0E442")


IRATNUniformPower <- ggplot(small_tab, aes(x = as.numeric(I), y = as.numeric(power), color = method)) +
  geom_line() +
  geom_point() +
  labs(x = "Number of items I", y = "Power Rate w", color = "Method") +
  theme_minimal(base_size = 14) +
  ylim(0, 1) +
  scale_x_continuous(breaks = unique(as.numeric(small_tab$I)))+
  scale_color_manual(values = okabe_ito[1:5], name= NULL) +
  legend_inside

IRATNUniformRejection <- ggplot(small_tab, aes(x = as.numeric(I), y = as.numeric(rejection), color = method)) +
  geom_line() +
  geom_point() +
  labs(x = "Number of Items I", y = "Rejection Rate r", color = "Method") +
  ylim(0, 0.1) + #To see better what is happening
  scale_x_continuous(breaks = unique(as.numeric(small_tab$I)))+
  theme_minimal(base_size = 14) +
  scale_color_manual(values = okabe_ito[1:5], name= NULL) +
  legend_inside


# — save individually —
ggsave(
  filename = "plots/IRATNUniformPower.png",
  plot     = IRATNUniformPower,
  width    = 6, 
  height   = 4, 
  dpi      = 300
)

ggsave(
  filename = "plots/IRATNUniformRejection.png",
  plot     = IRATNUniformRejection,
  width    = 6, 
  height   = 4, 
  dpi      = 300
)


# turn them into a single grob
IRATNUniform_combined <- grid.arrange(
  IRATNUniformRejection,
  IRATNUniformPower,
  ncol = 2  # two columns side by side
)

# save the combined plot
ggsave(
  filename = "plots/IRATNUniform-combined.png",
  plot     = IRATNUniform_combined,
  width    = 12,   # adjust as you like
  height   = 6,    # adjust as you like
  dpi      = 300
)

################################
# turn them into a single grob
IRATRejection_combined <- grid.arrange(
  IRATUniformRejection,
  IRATNUniformRejection,
  ncol = 2  # two columns side by side
)

# turn them into a single grob
IRATPower_combined <- grid.arrange(
  IRATUniformPower,
  IRATNUniformPower,
  ncol = 2  # two columns side by side
)

# save the combined plot
ggsave(
  filename = "plots/IRATRejection-combined.png",
  plot     = IRATRejection_combined,
  width    = 12,   # adjust as you like
  height   = 6,    # adjust as you like
  dpi      = 300
)

# save the combined plot
ggsave(
  filename = "plots/IRATPower-combined.png",
  plot     = IRATPower_combined,
  width    = 12,   # adjust as you like
  height   = 6,    # adjust as you like
  dpi      = 300
)


##########################################################################################
#SCENARIO 4: ratio rho is varied

#4.1 UNIFORM----------------------------------------------------------------------
rm(list = ls())
load("results/variable-ratio-res-U.RData")
load("results/variable-ratio-U.RData")



#Extracting information from tab
tab_df <- as.data.frame(matrix(unlist(tab), nrow = 9,ncol=21, byrow = TRUE))
tab_df<- t(tab_df)
notimes <- nrow(tab_df) / 3
methods <- rep(c("MH", "LR", "SIB"), notimes)
colnames(tab_df) <- c('M', 'N','rho','I', 'type','gamma', 'difSize', 'power', 'rejection')
tab_df <- data.frame(tab_df)
tab_df$method <- methods

small_tab <- tab_df[, c("rho", "power", "rejection", 'method')]
tab_df[, c("rho", "rejection", 'method')]

#Plotting
legend_inside <- theme(
  axis.text = element_text(size = 12),    # bigger x/y tick labels
  axis.title = element_text(size = 16),
  # bigger axis titles
  #legend.title = element_text(size = 14),
  #legend.text = element_text(size = 13),
  legend.position = c(0.99,0.01),
  legend.position.inside = "bottom-right",
  legend.justification = c("right", "bottom"),
  legend.background = element_rect(fill = alpha("white", 0.8), color = NA),
  legend.title = element_text(size = 10),
  legend.text = element_text(size = 9))

 okabe_ito <- c("#0072B2", "#E69F00", "#009E73", "#D55E00", "#CC79A7", "#F0E442")


 # reorder rho levels for plotting
 small_tab$rho <- factor(small_tab$rho, 
                         levels = c("1:4", "1:3", "1:2", "1:1", "2:1", "3:1", "4:1"))
 
 
rhoUniformPower <- ggplot(small_tab, aes(x = rho, y = as.numeric(power), color = method)) +
  geom_line(aes(group = method)) +  # explicitly connect lines by method
  geom_point() +
  labs(x = "Ratio ρ", y = "Power Rate w", color = "Method") +
  theme_minimal(base_size = 14) +
  ylim(0, 1) +
  scale_color_manual(values = okabe_ito[1:3], name= NULL) +
  legend_inside

rhoUniformRejection <- ggplot(small_tab, aes(x =rho, y = as.numeric(rejection), color = method)) +
  geom_line(aes(group = method)) +
  geom_point() +
  labs(x = "Ratio ρ", y = "Rejection Rate r", color = "Method") +
  theme_minimal(base_size = 14) +
  ylim(0, 0.1) +
  scale_color_manual(values = okabe_ito[1:3], name= NULL) +
  legend_inside

#Saving-----------------------------------------------------------------------
# — save individually —
ggsave(
  filename = "plots/rhoUniformPower.png",
  plot     = rhoUniformPower,
  width    = 6, 
  height   = 4, 
  dpi      = 300
)

ggsave(
  filename = "plots/rhoUniformRejection.png",
  plot     = rhoUniformRejection,
  width    = 6, 
  height   = 4, 
  dpi      = 300
)


# turn them into a single grob
rhoUniform_combined <- grid.arrange(
  rhoUniformRejection,
  rhoUniformPower,
  ncol = 2  # two columns side by side
)

# save the combined plot
ggsave(
  filename = "plots/rhoUniform-combined.png",
  plot     = rhoUniform_combined,
  width    = 12,   # adjust as you like
  height   = 6,    # adjust as you like
  dpi      = 300
)

#4.2 NON-UNIFORM----------------------------------------------------------------
rm(list = ls())
#Data loading
load("results/variable-ratio-res-NONU.RData")
load("results/variable-ratio-NONU.RData")

#Extracting information from tab
tab_df <- as.data.frame(matrix(unlist(tab), nrow = 9,ncol=42, byrow = TRUE))
tab_df<- t(tab_df)
notimes <- nrow(tab_df) / 6
methods <- rep(c("MH","MH-NU","MH-L", "MH-H", "LR", "CSIB"), notimes)
colnames(tab_df) <- c('M', 'N',"rho", 'I', 'type', 'gamma', 'difSize', 'power', 'rejection')
tab_df <- data.frame(tab_df)
tab_df$method <- methods

small_tab <- tab_df[, c("rho", "power", "rejection", 'method')]
small_tab <- small_tab[!small_tab$method == "MH", ]
tab_df[, c("rho", "rejection", 'method')]
#Plotting
legend_inside <- theme(
  axis.text = element_text(size = 12),    # bigger x/y tick labels
  axis.title = element_text(size = 16),
  # bigger axis titles
  #legend.title = element_text(size = 14),
  #legend.text = element_text(size = 13),
  legend.position = c(0.99,0.01),
  legend.position.inside = "bottom-right",
  legend.justification = c("right", "bottom"),
  legend.background = element_rect(fill = alpha("white", 0.8), color = NA),
  legend.title = element_text(size = 10),
  legend.text = element_text(size = 9))

 okabe_ito <- c("#0072B2", "#E69F00", "#009E73", "#D55E00", "#CC79A7", "#F0E442")
 small_tab$rho <- factor(small_tab$rho, 
                         levels = c("1:4", "1:3", "1:2", "1:1", "2:1", "3:1", "4:1"))
 
 


rhoNUniformPower <- ggplot(small_tab, aes(x = rho, y = as.numeric(power), color = method)) +
  geom_line(aes(group = method)) +  # explicitly connect lines by method
  geom_point() +
  labs(x = "Ratio ρ", y = "Power Rate w", color = "Method") +
  theme_minimal(base_size = 14) +
  ylim(0, 1) +
  scale_color_manual(values = okabe_ito[1:5], name= NULL) +
  legend_inside

rhoNUniformRejection <- ggplot(small_tab, aes(x =rho, y = as.numeric(rejection), color = method)) +
  geom_line(aes(group = method)) +
  geom_point() +
  labs(x = "Ratio ρ", y = "Rejection Rate r", color = "Method") +
  theme_minimal(base_size = 14) +
  ylim(0, 0.1) +
  scale_color_manual(values = okabe_ito[1:5], name= NULL) +
  legend_inside

#Saving-----------------------------------------------------------------------
# — save individually —
ggsave(
  filename = "plots/rhoNUniformPower.png",
  plot     = rhoNUniformPower,
  width    = 6, 
  height   = 4, 
  dpi      = 300
)

ggsave(
  filename = "plots/rhoNUniformRejection.png",
  plot     = rhoNUniformRejection,
  width    = 6, 
  height   = 4, 
  dpi      = 300
)


# turn them into a single grob
rhoNUniform_combined <- grid.arrange(
  rhoNUniformRejection,
  rhoNUniformPower,
  ncol = 2  # two columns side by side
)

# save the combined plot
ggsave(
  filename = "plots/rhoNUniform-combined.png",
  plot     = rhoNUniform_combined,
  width    = 12,   # adjust as you like
  height   = 6,    # adjust as you like
  dpi      = 300
)

# turn them into a single grob
rhoRejection_combined <- grid.arrange(
  rhoUniformRejection,
  rhoNUniformRejection,
  ncol = 2  # two columns side by side
)

# turn them into a single grob
rhoPower_combined <- grid.arrange(
  rhoUniformPower,
  rhoNUniformPower,
  ncol = 2  # two columns side by side
)

# save the combined plot
ggsave(
  filename = "plots/rhoRejection-combined.png",
  plot     = rhoRejection_combined,
  width    = 12,   # adjust as you like
  height   = 6,    # adjust as you like
  dpi      = 300
)

# save the combined plot
ggsave(
  filename = "plots/rhoPower-combined.png",
  plot     = rhoPower_combined,
  width    = 12,   # adjust as you like
  height   = 6,    # adjust as you like
  dpi      = 300
)





################################################################################
#SCENARIO 5: DIFsize

#5.1 UNIFORM--------------------------------------------------------------
#Data loading
rm(list = ls())
load("results/variable-difSize-res-U.RData")
load("results/variable-difSize-U.RData")

#Extracting information from tab
tab_df <- as.data.frame(matrix(unlist(tab), nrow = 8,ncol=24, byrow = TRUE))
tab_df<- t(tab_df)
notimes <- nrow(tab_df) / 3
methods <- rep(c("MH", "LR", "SIB"), notimes)
colnames(tab_df) <- c('M', 'N', 'I', 'type', 'gamma', 'difSize', 'power', 'rejection')
tab_df <- data.frame(tab_df)
tab_df$method <- methods

small_tab <- tab_df[, c("difSize", "power", "rejection", 'method')]
small_tab <-small_tab[7:24,]
tab_df[, c("difSize", "rejection", 'method')]

#Plotting
legend_inside <- theme(
  axis.text = element_text(size = 12),    # bigger x/y tick labels
  axis.title = element_text(size = 16),
  # bigger axis titles
  #legend.title = element_text(size = 14),
  #legend.text = element_text(size = 13),
  legend.position = c(0.99,0.01),
  legend.position.inside = "bottom-right",
  legend.justification = c("right", "bottom"),
  legend.background = element_rect(fill = alpha("white", 0.8), color = NA),
  legend.title = element_text(size = 10),
  legend.text = element_text(size = 9))

 okabe_ito <- c("#0072B2", "#E69F00", "#009E73", "#D55E00", "#CC79A7", "#F0E442")


difSizeUniformPower <- ggplot(small_tab, aes(x = as.numeric(difSize), y = as.numeric(power), color = method)) +
  geom_line() +
  geom_point() +
  labs(x = "Size of DIF δ", y = "Power Rate w", color = "Method") +
  theme_minimal(base_size = 14) +
  ylim(0, 1) +
  scale_x_continuous(breaks = unique(as.numeric(small_tab$difSize)))+
  scale_color_manual(values = okabe_ito[1:3], name= NULL) +
  legend_inside

difSizeUniformRejection <- ggplot(small_tab, aes(x = as.numeric(difSize), y = as.numeric(rejection), color = method)) +
  geom_line() +
  geom_point() +
  labs(x = "Size of DIF δ", y = "Rejection Rate r", color = "Method") +
  ylim(0, 0.1) + #To see better what is happening
  scale_x_continuous(breaks = unique(as.numeric(small_tab$difSize)))+
  theme_minimal(base_size = 14) +
  scale_color_manual(values = okabe_ito[1:3], name= NULL) +
  legend_inside


# — save individually —
ggsave(
  filename = "plots/difSizeUniformPower.png",
  plot     = difSizeUniformPower,
  width    = 6, 
  height   = 4, 
  dpi      = 300
)

ggsave(
  filename = "plots/difSizeUniformRejection.png",
  plot     = difSizeUniformRejection,
  width    = 6, 
  height   = 4, 
  dpi      = 300
)


# turn them into a single grob
difSizeUniform_combined <- grid.arrange(
  difSizeUniformRejection,
  difSizeUniformPower,
  ncol = 2  # two columns side by side
)

# save the combined plot
ggsave(
  filename = "plots/difSizeUniform-combined.png",
  plot     = difSizeUniform_combined,
  width    = 12,   # adjust as you like
  height   = 6,    # adjust as you like
  dpi      = 300
)

#5.2 NON-UNIFORM----------------------------------------------------------------
rm(list = ls())
#Data loading
load("results/variable-difSize-res-NONU.RData")
load("results/variable-difSize-NONU.RData")

tab
#Extracting information from tab
tab_df <- as.data.frame(matrix(unlist(tab), nrow = 8,ncol=36,byrow = TRUE))
tab_df<- t(tab_df)
notimes <- nrow(tab_df) / 6
methods <- rep(c("MH", "MH-NU","MH-L", "MH-H", "LR", "CSIB"), notimes)
colnames(tab_df) <- c('M', 'N', 'I', 'type', 'gamma', 'difSize', 'power', 'rejection')
tab_df <- data.frame(tab_df)
tab_df$method <- methods

small_tab <- tab_df[, c("difSize", "power", "rejection", 'method')]
small_tab <- small_tab[!small_tab$method == "MH", ]
tab_df[, c("difSize", "rejection", 'method')]
#Plotting
legend_inside <- theme(
  axis.text = element_text(size = 12),    # bigger x/y tick labels
  axis.title = element_text(size = 16),
  # bigger axis titles
  #legend.title = element_text(size = 14),
  #legend.text = element_text(size = 13),
  legend.position = c(0.99,0.01),
  legend.position.inside = "bottom-right",
  legend.justification = c("right", "bottom"),
  legend.background = element_rect(fill = alpha("white", 0.8), color = NA),
  legend.title = element_text(size = 10),
  legend.text = element_text(size = 9))

 okabe_ito <- c("#0072B2", "#E69F00", "#009E73", "#D55E00", "#CC79A7", "#F0E442")



difSizeNUniformPower <- ggplot(small_tab, aes(x = as.numeric(difSize), y = as.numeric(power), color = method)) +
  geom_line() +
  geom_point() +
  labs(x = "Size of DIF Δ", y = "Power Rate w", color = "Method") +
  theme_minimal(base_size = 14) +
  ylim(0, 1) + 
  scale_x_continuous(breaks = unique(as.numeric(small_tab$difSize)))+
  scale_color_manual(values = okabe_ito[1:5], name= NULL) +
  legend_inside

difSizeNUniformRejection <- ggplot(small_tab, aes(x = as.numeric(difSize), y = as.numeric(rejection), color = method)) +
  geom_line() +
  geom_point() +
  labs(x = "Size of DIF Δ", y = "Rejection Rate r", color = "Method") +
  ylim(0, 0.1) + #To see better what is happening
  theme_minimal(base_size = 14) +
  scale_x_continuous(breaks = unique(as.numeric(small_tab$difSize)))+
  scale_color_manual(values = okabe_ito[1:5], name= NULL) +
  legend_inside


# — save individually —
ggsave(
  filename = "plots/difSizeNUniformPower.png",
  plot     = difSizeNUniformPower,
  width    = 6, 
  height   = 4, 
  dpi      = 300
)

ggsave(
  filename = "plots/difSizeNUniformRejection.png",
  plot     = difSizeNUniformRejection,
  width    = 6, 
  height   = 4, 
  dpi      = 300
)


# turn them into a single grob
difSizeNUniform_combined <- grid.arrange(
  difSizeNUniformRejection,
  difSizeNUniformPower,
  ncol = 2  # two columns side by side
)

# save the combined plot
ggsave(
  filename = "plots/difSizeNUniform-combined.png",
  plot     = difSizeNUniform_combined,
  width    = 12,   # adjust as you like
  height   = 6,    # adjust as you like
  dpi      = 300
)



# turn them into a single grob
difSizeRejection_combined <- grid.arrange(
  difSizeUniformRejection,
  difSizeNUniformRejection,
  ncol = 2  # two columns side by side
)

# turn them into a single grob
difSizePower_combined <- grid.arrange(
  difSizeUniformPower,
  difSizeNUniformPower,
  ncol = 2  # two columns side by side
)

# save the combined plot
ggsave(
  filename = "plots/difSizeRejection-combined.png",
  plot     = difSizeRejection_combined,
  width    = 12,   # adjust as you like
  height   = 6,    # adjust as you like
  dpi      = 300
)



# save the combined plot
ggsave(
  filename = "plots/ifSizePower-combined.png",
  plot     = difSizePower_combined,
  width    = 12,   # adjust as you like
  height   = 6,    # adjust as you like
  dpi      = 300
)



##################################################################################
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
colnames(tab_df) <- c('M', 'N', 'I', 'type', 'gamma', 'difSize', 'power', 'rejection')
tab_df <- data.frame(tab_df)
#tab_df <- tab_df[112:153,]
notimes <- nrow(tab_df) / 3
methods <- rep(c("MH", "LR", "SIB"), notimes)
tab_df$method <- methods

small_tab <- tab_df[, c("N", "power", "rejection", 'method')]


#Plotting
PUniformPower <- ggplot(small_tab, aes(x = as.numeric(N), y = as.numeric(power), color = method)) +
  geom_line() +
  geom_point() +
  labs(x = "Number of respondents N (ratio)", y = "Power Rate w", color = "Method") +
  theme_minimal(base_size = 14) +
  legend_inside

PUniformRejection <- ggplot(tab_df, aes(x = as.numeric(N), y = as.numeric(rejection), color = method)) +
  geom_line() +
  geom_point() +
  labs(x = "Number of respondents N (ratio)", y = "Rejection Rate r", color = "Method") +
  ylim(0, 0.1) + #To see better what is happening
  theme_minimal(base_size = 14) +
  legend_inside

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
colnames(tab_df) <- c('M', 'N', 'I', 'type', 'gamma', 'difSize', 'power', 'rejection')
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
  labs(x = "Number of respondents N (ratio)", y = "Power Rate w", color = "Method") +
  theme_minimal(base_size = 14) +
  legend_inside

PNUniformRejection <- ggplot(tab_df, aes(x = as.numeric(N), y = as.numeric(rejection), color = method)) +
  geom_line() +
  geom_point() +
  labs(x = "Number of respondents N (ratio)", y = "Rejection Rate r", color = "Method") +
  ylim(0, 0.1) + #To see better what is happening
  theme_minimal(base_size = 14) +
  legend_inside

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






# Itemwise plotting ------------------------------------------------------

# 1. get unique n’s (as characters, to index res_total)
unique_ns <- unique(tab_df$N)

# 2. for each n, run the function and grab only method + vote_ratio
results_extracted <- do.call(rbind, lapply(unique_ns, function(n_val) {
  # compute
  tmp <- calculate_itemwise_detection(res_total[[ as.character(n_val) ]])
  # select and add an 'n' column
  data.frame(
    N          = n_val,
    method     = tmp$method,
    vote_ratio = tmp$vote_ratio,
    stringsAsFactors = FALSE
  )
}))

# 3. inspect
head(results_extracted)

#same for one item
ggplot(results_extracted, aes(x = as.numeric(N), y = as.numeric(vote_ratio), color = method)) +
  geom_line() +
  geom_point() +
  labs(x = "Number of respondents N", y = "Vote Ratio", color = "Method") +
  theme_minimal(base_size = 14) +
  legend_inside



