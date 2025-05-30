#Simulation from ground up
library(ltm)
library(difR)
#rm(list = ls())
###################################################################################
#Hyperparameters

N = 10 #number of simulations
m = 5 #number of items (should change to I)
n_tot = 1000 # number of respondents (700, 1000, 1200, 1500) 
rat_n = c(1:1) #(c(5,2), c(1,1), c(5,1), c(2,1))
set.seed(123)
#------------------------------------------------------------------------------
#Functions by piece
generate_group_sizes <- function(n_total, ratio = c(1, 1)) {#Reference: Focal
  if (length(ratio) != 2) stop("Ratio must be of length 2, like c(1, 3) or 1:3.")
  
  total_parts <- sum(ratio)
  n2 <- round(n_total * ratio[1] / total_parts)  # focal group (first in ratio)
  n1 <- n_total - n2                             # reference group (second in ratio)
  
  return(list(n = n_total, n1 = n1, n2 = n2))
}

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

simulation_abe <- function(N,n1,n2,skup, param1,param2,theta1, theta2){
  # Number of respondents per group
  #n1 <- length(theta1)
  #n2 <- length(theta2)
  
  # Initialize result storage
  Mantel     <- Bres <- Bres2 <- Log <- rep(NA, N)
  MantelStat <- BresStat <- Bres2Stat <- LogStat <- rep(NA, N)
  for(i in 1:N){
    data1=rmvlogis(n1, param1,IRT=TRUE, link="logit",z.vals=matrix(theta1)) #generates 0,1 asnwer for each responndet using 3PL model (b+c step)
    data2=rmvlogis(n2, param2,IRT=TRUE, link="logit",z.vals=matrix(theta2))
    data=rbind(data1,data2)
    
    MantelStat[i]=help=try(difMH(data, skup,focal.name=1, correct=TRUE,MHstat="MHChisq", exact=FALSE)[[1]][1],silent=TRUE)
    Mantel[i]=ifelse(is.numeric(help),1-pchisq(help,1),NA)
    BresStat[i]=help=try(difBD(data, skup,focal.name=1)[[1]][1,3],silent=TRUE)
    Bres[i]=ifelse(is.numeric(help),help,NA)
    Bres2Stat[i]=help=try(difBD(data, skup,focal.name=1,BDstat="trend")[[1]][1,3],silent=TRUE)
    Bres2[i]=ifelse(is.numeric(help),help,NA)
    LogStat[i]=L=try(difLogistic(data, skup,focal.name=1)[[1]][1],silent=TRUE)
    Log[i]=ifelse(is.numeric(L)==TRUE,1-pchisq(L,2),NA)
    print(paste(i,date()))
  }
  mat=rbind(Mantel,Bres,Bres2,Log)
  matStat=rbind(MantelStat,BresStat,Bres2Stat,LogStat)
  return(mat)
}


simul_total <- function(N, n_total, rat_n, I, mu_R, mu_F, type = "alpha", diffs = NULL) {
  # Group sizes
  list_ns <- generate_group_sizes(n_total, rat_n)
  n <- list_ns$n
  n1 <- list_ns$n1
  n2 <- list_ns$n2
  
  # Abilities
  theta1_ab <- rnorm(n1, mean = mu_R, sd = 1)
  theta2_ab <- rnorm(n2, mean = mu_F, sd = 1)
  
  # Item parameters (same initially)
  param1_ab <- generate_param(I)
  param2_ab <- param1_ab
  skup <- c(rep(0, n1), rep(1, n2))
  
  # Init result list
  results <- list()
  
  if (type == "alpha") {
    # Type I error — no DIF introduced
    result_alpha <- simulation_abe(N,n1,n2,skup, param1_ab, param2_ab, theta1_ab, theta2_ab)
    results[["alpha"]] <- result_alpha
    
  } else if (type == "unif") {
    # Uniform DIF — manipulate difficulty (b)
    for (i in seq_along(diffs)) {
      d <- diffs[i]
      param2_ab_mod <- param1_ab
      param2_ab_mod[1, 1] <- param1_ab[1, 1] + d # shift difficulty
      
      label <- paste0("Unif_DIF_d", d)
      results[[label]] <- simulation_abe(N,n1,n2,skup, param1_ab, param2_ab_mod, theta1_ab, theta2_ab)
    }
    
  } else if (type == "nonunif") {
    # Non-uniform DIF — manipulate discrimination (a)
    for (i in seq_along(diffs)) {
      delta <- diffs[i]
      param2_ab_mod <- param1_ab
      a_orig <- param1_ab[1, 2]
      c_orig <- param1_ab[1, 3]
      param2_ab_mod[1, 2] <- 2 * a_orig / (2 + delta * a_orig / ((1 - c_orig) * log(2)))
      
      label <- paste0("NonUnif_DIF_delta", delta)
      results[[label]] <- simulation_abe(N,n1,n2,skup, param1_ab, param2_ab_mod, theta1_ab, theta2_ab)
    }
  }
  
  return(results)
}
#Evaluating------------------------------------------------------------

# Type I error
res_alpha0 <- simul_total(N = 10, n_total = 1000, rat_n = c(1,1), I = 5, mu_R = 0, mu_F = -1, type = "alpha")

# Uniform DIF
res_unif0 <- simul_total(N = 10, n_total = 1000, rat_n = c(5,2), I = 5, mu_R = 0, mu_F = -1, 
                        type = "unif", diffs = c(0.5,1,2,4))

# Non-uniform DIF
res_nonunif0 <- simul_total(N = 10, n_total = 1000, rat_n = c(2,1), I = 5, mu_R = 0, mu_F = -1, 
                           type = "nonunif", diffs = c(0.4,0.6,0.8,1))


# to do --------------------------------------------------------------
# multiple items ???
# add sibtest and NU mantel haensel

# look into specific libraries 

?difMH
?difSIBTEST

#####################################################################
