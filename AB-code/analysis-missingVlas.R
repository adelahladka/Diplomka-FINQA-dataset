samples_n <- c(10, 20, 30, 40, 50, 60, 70, 80, 100, 200, 400, 500, 1000, 1200, 2000)
stats <- c("LogStat", "MantelStat", "SIBStat")
samples_n <- cbind(c(1,1), c(2,1), c(3,1), c(4,1), c(1,2), c(1,3), c(1,4))
samples_n <- c(0.4, 0.6, 0.8, 1.0, 1.2, 1.5)
samples_n <- c(20, 40, 50, 60, 70, 80, 100, 200)
res_total <- res_total2
#UNIFORM------------------------------------------------------------------------

# Count NAs
cat("NA counts:\n")
for (stat in stats) {
  cat(paste0("\n", stat, ":\n"))
  for (n in samples_n) {
    count_na <- sum(is.na(unlist(res_total[[as.character(n)]]$results[[stat]])))
    cat(paste0("  N = ", n, ": ", count_na, "\n"))
  }
}

# Count Infs
cat("\nInfinite value counts:\n")
for (stat in stats) {
  cat(paste0("\n", stat, ":\n"))
  for (n in samples_n) {
    count_inf <- sum(is.infinite(unlist(res_total[[as.character(n)]]$results[[stat]])))
    cat(paste0("  N = ", n, ": ", count_inf, "\n"))
  }
}
cat("\nZero value counts:\n")
for (stat in stats) {
  cat(paste0("\n", stat, ":\n"))
  for (n in samples_n) {
    count_zero <- sum(unlist(res_total[[as.character(n)]]$results[[stat]]) == 0, na.rm = TRUE)
    cat(paste0("  N = ", n, ": ", count_zero, "\n"))
  }
}


#NON UNIFORM---------------------------------------------------------------------

# Define your stats from the screenshot
stats <- c("MantelNormalStat", "MantelLowStat", "MantelHighStat", "LogStat", "cSIBStat")

# Count NAs
cat("NA counts:\n")
for (stat in stats) {
  cat(paste0("\n", stat, ":\n"))
  for (n in samples_n) {
    count_na <- sum(is.na(unlist(res_total[[as.character(n)]]$results[[stat]])))/40000*100
    cat(paste0("  N = ", n, ": ", count_na, "\n"))
  }
}

# Count Infs
cat("\nInfinite value counts:\n")
for (stat in stats) {
  cat(paste0("\n", stat, ":\n"))
  for (n in samples_n) {
    count_inf <- sum(is.infinite(unlist(res_total[[as.character(n)]]$results[[stat]])))
    cat(paste0("  N = ", n, ": ", count_inf, "\n"))
  }
}

# Count Zeros
cat("\nZero value counts:\n")
for (stat in stats) {
  cat(paste0("\n", stat, ":\n"))
  for (n in samples_n) {
    count_zero <- sum(unlist(res_total[[as.character(n)]]$results[[stat]]) == 0, na.rm = TRUE)
    cat(paste0("  N = ", n, ": ", count_zero, "\n"))
  }
}


########################################################################
item_params<- res_total$`10`$metadata$paramR
# Create the item parameters data frame
item_params <- as.data.frame(item_params)
item_params$Item <- rownames(item_params)
item_params <- item_params[, c("Item", setdiff(names(item_params), "Item"))]
# Sort by discrimination parameter 'a'
item_params_sorted <- item_params[order(item_params$a), ]
item_params_sorted$rank <- 1:nrow(item_params_sorted)
item_params_sorted$Item
# Create the plot
library(ggplot2)

# Basic plot
p1 <- ggplot(item_params_sorted, aes(x = reorder(gsub("Item", "", Item), rank), y = a)) +
  geom_point(aes(color = ifelse(Item == "Item1", "#E69F00", "#0072B2")), 
             size = 3, alpha = 0.7) +
  geom_line(aes(group = 1), color = "#0072B2", alpha = 0.5) +
  scale_color_identity() +
  labs(
    x = "Item Numbers (sorted by discrimination)",
    y = "Discrimination Parameter (a)"
  ) +
  theme_minimal()
