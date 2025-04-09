#Simulation from ground up
library(ltm)
library(difR)
###################################################################################
N = 10 #number of simulations
m = 5 #number of items (should change to I)
set.seed(123)

#------------------------------------------------------------------------------
#1.choice
#n1 = 500 #number of  reference folks
#n2 = 200 #number of focal folks

#2.choice
#n1 = 500
#n2 = 500
#3.choice
#n1 = 1000
#n2 = 200
#4.choice
#n1 = 1000
#n2 = 500

#Automatic generation 
generate_group_sizes <- function(n_total, ratio = c(1, 1)) {#Reference: Focal
  if (length(ratio) != 2) stop("Ratio must be of length 2, like c(1, 3) or 1:3.")
  
  total_parts <- sum(ratio)
  n2 <- round(n_total * ratio[1] / total_parts)  # focal group (first in ratio)
  n1 <- n_total - n2                             # reference group (second in ratio)
  
  return(list(n = n_total, n1 = n1, n2 = n2))
}


list_number <- generate_group_sizes(800, c(1,3))
n = list_number$n
n1 = list_number$n1
n2 = list_number$n2

#n  = n1 + n2 # total number of respondents
theta1_ab <- rnorm(n1, mean = 0, sd = 1) # ability for reference group
theta2_ab <-rnorm(n2, mean = -1, sd = 1) #??? do I need to generate differently? -> Penfield gave -1
  


generate_param <- function(I) {
  # Generate parameters for the 3PL model
  b <- rnorm(I, mean = 0, sd = 1)
  z <- rnorm(I, mean = 0, sd = sqrt(0.1225))
  a <- exp(z)
  c <- rep(0.2, I)
  
  param <- cbind(b, a, c)
  rownames(param) <- paste0("Item", 1:I)
  colnames(param) <- c("b", "a", "c")
  return(param)
}


param1_ab <- generate_param(m)       # No DIF
param2_ab <-param1_ab
#param2_ab <- generate_param(5)            # or do I need to param2_ab <-param1_ab?

skup <- c(rep(0,n1),rep(1,n2))

#Generating 
simulace_alt <- function(N,param1,param2,theta1, theta2, text){
  for(i in 1:N){
    data1=rmvlogis(n1, param1,IRT=TRUE, link="logit",z.vals=matrix(theta1)) #generates prob of the asnwer for each responndet using 3PL model
    data2=rmvlogis(n2, param2,IRT=TRUE, link="logit",z.vals=matrix(theta2))
    data=rbind(data1,data2)
    
    
    #paramery pro 1pl model
    S11pl<-try(itemParEst(data1,model="1PL"),silent=TRUE)
    S21pl=try(itemParEst(data2,model="1PL"),silent=TRUE)
    if(is.matrix(S11pl)){
      beta1pl[1:m,i]=S11pl[,1]
      SDb1pl[1:m,i]=S11pl[,2]}
    if(is.matrix(S21pl)){
      beta1pl[(m+1):(2*m),i]=S21pl[,1]
      SDb1pl[(m+1):(2*m),i]=S21pl[,2]
      
    }
    #parametry pro 2pl model
    S12pl=try(itemParEst(data1,model="2PL"),silent=TRUE)
    S22pl=try(itemParEst(data2,model="2PL"),silent=TRUE)
    if(is.matrix(S12pl)){
      beta2pl[1:m,i]=S12pl[,2]
      alpha2pl[1:m,i]=S12pl[,1]
      SDb2pl[1:m,i]=S12pl[,4]
      SDa2pl[1:m,i]=S12pl[,3]
      cov2pl[1:m,i]=S12pl[,5]
    }
    if(is.matrix(S22pl)){
      beta2pl[(m+1):(2*m),i]=S22pl[,2]
      alpha2pl[(m+1):(2*m),i]=S22pl[,1]
      SDb2pl[(m+1):(2*m),i]=S22pl[,4]
      SDa2pl[(m+1):(2*m),i]=S22pl[,3]
      cov2pl[(m+1):(2*m),i]=S22pl[,5]
    }
    #parametry pro 3pl model
    S3pl=try(itemPar3PL(data),silent=TRUE)
    if(is.matrix(S3pl)){
      c=S3pl[,3]
      gamma3pl[,i]=c
      S13pl=try(itemParEst(data1,c=c,model="3PL"),silent=TRUE)
      S23pl=try(itemParEst(data2,c=c,model="3PL"),silent=TRUE)
      if(is.matrix(S13pl)){
        beta3pl[1:m,i]=S13pl[,2]
        alpha3pl[1:m,i]=S13pl[,1]
        SDb3pl[1:m,i]=S13pl[,4]
        SDa3pl[1:m,i]=S13pl[,3]
        cov3pl[1:m,i]=S13pl[,5]
      }
      if(is.matrix(S23pl)){
        beta3pl[(m+1):(2*m),i]=S23pl[,1]
        alpha3pl[(m+1):(2*m),i]=S23pl[,2]
        SDb3pl[(m+1):(2*m),i]=S23pl[,4]
        SDa3pl[(m+1):(2*m),i]=S23pl[,3]
        cov3pl[(m+1):(2*m),i]=S23pl[,5]
        
      }
    }
    
    if(is.matrix(S11pl)&is.matrix(S21pl)){
      Lord1Stat=L=try(difLord(irtParam=rbind(S11pl,S21pl),model="1PL",same.scale=FALSE)[[1]][1],silent=TRUE)#asi vyuzijeme vsetky p-hodnoty
      Lord1[i]=ifelse(is.numeric(L)==TRUE,1-pchisq(L,1),NA) 
      RajuU1Stat[i]=R=try(difRaju(irtParam=rbind(S11pl,S21pl),model="1PL",same.scale=FALSE,signed=FALSE)[[1]][1],silent=TRUE)
      
      RajuU1[i]=ifelse(is.numeric(R)==TRUE,2*pnorm(-abs(R),lower.tail=TRUE,FALSE),NA)
      RajuS1Stat[i]=R2=try(difRaju(irtParam=rbind(S11pl,S21pl),model="1PL",same.scale=FALSE,signed=TRUE)[[1]][1],silent=TRUE)
      RajuS1[i]=ifelse(is.numeric(R2)==TRUE,2*pnorm(-abs(R2),lower.tail=TRUE,FALSE),NA)   
    }
    if(is.matrix(S12pl)&is.matrix(S22pl)){
      Lord2Stat[i]=L=try(difLord(irtParam=rbind(S12pl,S22pl),model="2PL",same.scale=FALSE)[[1]][1],silent=TRUE)#asi vyuzijeme vsetky p-hodnoty
      Lord2[i]=ifelse(is.numeric(L)==TRUE,1-pchisq(L,2),NA) 
      RajuU2Stat[i]=R=try(difRaju(irtParam=rbind(S12pl,S22pl),model="2PL",same.scale=FALSE,signed=FALSE)[[1]][1],silent=TRUE)
      RajuU2[i]=ifelse(is.numeric(R)==TRUE,2*pnorm(-abs(R),lower.tail=TRUE,FALSE),NA)
      
      RajuS2Stat[i]=R2=try(difRaju(irtParam=rbind(S12pl,S22pl),model="2PL",same.scale=FALSE,signed=TRUE)[[1]][1],silent=TRUE)
      RajuS2[i]=ifelse(is.numeric(R2)==TRUE,2*pnorm(-abs(R2),lower.tail=TRUE,FALSE),NA)   
    }
    if(is.numeric(c)){
      if(is.matrix(S13pl)&is.matrix(S23pl)){
        Lord3Stat[i]=L=try(difLord(irtParam=rbind(S13pl,S23pl),model="3PL",c=c,same.scale=FALSE)[[1]][1],silent=TRUE)#asi vyuzijeme vsetky p-hodnoty
        Lord3[i]=ifelse(is.numeric(L)==TRUE,1-pchisq(L,2),NA) 
        RajuU3Stat[i]=R=try(difRaju(irtParam=rbind(S13pl,S23pl),model="3PL",c=c,same.scale=FALSE,signed=FALSE)[[1]][1],silent=TRUE)
        RajuU3[i]=ifelse(is.numeric(R)==TRUE,2*pnorm(-abs(R),lower.tail=TRUE,FALSE),NA)
        
        RajuS3Stat[i]=R2=try(difRaju(irtParam=rbind(S13pl,S23pl),model="3PL",c=c,same.scale=FALSE,signed=TRUE)[[1]][1],silent=TRUE)
        RajuS3[i]=ifelse(is.numeric(R2)==TRUE,2*pnorm(-abs(R2),lower.tail=TRUE,FALSE),NA)           
      }
    }
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
  mat=rbind(Lord1,Lord2,Lord3,RajuU1,RajuU2,RajuU3,RajuS1,RajuS2,RajuS3,Mantel,Bres,Bres2,Log)
  matStat=rbind(Lord1Stat,Lord2Stat,Lord3Stat,RajuU1Stat,RajuU2Stat,RajuU3Stat,RajuS1Stat,RajuS2Stat,
                RajuS3Stat,MantelStat,BresStat,Bres2Stat,LogStat)
  return(mat)
}

Lord1=Lord2=Lord3=RajuU1=RajuU2=RajuU3=RajuS1=RajuS2=RajuS3=Mantel=Stand=Bres=Bres2=Log=rep(NA,N)

Lord1Stat=Lord2Stat=Lord3Stat=RajuU1Stat=RajuU2Stat=RajuU3Stat=RajuS1Stat=RajuS2Stat=RajuS3Stat=MantelStat=
  StandStat=BresStat=Bres2Stat=LogStat=rep(NA,N)

beta1pl=beta2pl=beta3pl=alpha2pl=alpha3pl=matrix(NA,2*m,N)

SDb1pl=SDb2pl=SDb3pl=cov2pl=cov3pl=SDa2pl=SDa3pl=matrix(NA,2*m,N)

gamma3pl=matrix(NA,m,N)

text=paste("Hladina",sep="") #argument for exporting into csv
simul_alpha=simulace_alt(N,param1_ab,param2_ab,theta1_ab, theta2_ab, text) #delta = -> assesing type I error, why? (Penfield 2003)
View(simul_alpha)


#sila-uniformni DIF
rozdil=c(0.5,1,2,4)
param2_ab_unif <- param2_ab
for(i in 1:4){
 text=paste("Sila",rozdil[i],sep="")
  param2_ab_unif[1,1]=param1_ab[1,1]+rozdil[i] #shifting difficulty of the first item
  simul_power_unif=simulace_alt(N,param1_ab,param2_ab_unif,theta1_ab, theta2_ab, text) 
}

#sila-neuniformni DIF
delta=c(0.4,0.6,0.8,1)
param2_ab_nenif <- param2_ab

for(i in 1:5){
  param2_ab_nenif[1,2]=2*param1_ab[1,2]/(2+delta[i]*param1_ab[1,2]/((1-param1_ab[1,3])*log(2)))
  text=paste("NeuniD",delta[i],sep="")
  simul_power_nounif=simulace_alt(N,param1_ab,param2_ab_nenif, theta1_ab, theta2_ab, text)
}

# potential to do
# multiple items
# improve code
# look into specific libraries 



# param1_ab <- matrix(c(
#   0, 1, 0.2,   # item 1: b=0, a=1, c=0.2
#   0.5, 1.2, 0.25, # item 2
#   -1, 0.8, 0.1,   # item 3
#   1.5, 1.5, 0.3,  # item 4
#   0.2, 1.0, 0.15  # item 5
# ), nrow = 5, byrow = TRUE) #for not some reasonable params
# param2_ab <- param1_ab #no dif null hyphothesis
#Additional how to genereate params