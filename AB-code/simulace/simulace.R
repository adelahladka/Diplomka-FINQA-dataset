library(ltm)
library(difR)
simulace=function(N,param1,param2,text){
  for(i in 1:N){
    data1=rmvlogis(n1, param1,IRT=TRUE, link="logit",z.vals=matrix(theta1))
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
  write.csv2(matStat,paste("Statistika",m,"n",n1+n2,text,".csv",sep=""),row.names=FALSE)
  write.csv2(beta1pl,paste("beta1plm",m,"n",n1+n2,text,".csv",sep=""),row.names=FALSE)
  write.csv2(beta2pl,paste("beta2plm",m,"n",n1+n2,text,".csv",sep=""),row.names=FALSE)
  write.csv2(beta3pl,paste("beta3plm",m,"n",n1+n2,text,".csv",sep=""),row.names=FALSE)
  
  
  write.csv2(alpha2pl,paste("alpha2plm",m,"n",n1+n2,text,".csv",sep=""),row.names=FALSE)
  write.csv2(alpha3pl,paste("alpha3plm",m,"n",n1+n2,text,".csv",sep=""),row.names=FALSE)
  
  write.csv2(gamma3pl,paste("gamma3plm",m,"n",n1+n2,text,".csv",sep=""),row.names=FALSE)
  write.csv2(mat,paste(text,"m",m,"n",n1+n2,".csv",sep=""),row.names=FALSE)
  
  write.csv2(SDb1pl,paste("SDb1plm",m,"n",n1+n2,text,".csv",sep=""),row.names=FALSE)
  write.csv2(SDb2pl,paste("SDb2plm",m,"n",n1+n2,text,".csv",sep=""),row.names=FALSE)
  write.csv2(SDb3pl,paste("SDb3plm",m,"n",n1+n2,text,".csv",sep=""),row.names=FALSE)
  write.csv2(SDa2pl,paste("SDa2plm",m,"n",n1+n2,text,".csv",sep=""),row.names=FALSE)
  write.csv2(SDa3pl,paste("SDa3plm",m,"n",n1+n2,text,".csv",sep=""),row.names=FALSE)
  write.csv2(cov2pl,paste("cov2plm",m,"n",n1+n2,text,".csv",sep=""),row.names=FALSE)
  write.csv2(cov3pl,paste("cov3plm",m,"n",n1+n2,text,".csv",sep=""),row.names=FALSE)
  return(mat)
}


m=80#pocet polozek
N=1000#pocet simulaci
Lord1=Lord2=Lord3=RajuU1=RajuU2=RajuU3=RajuS1=RajuS2=RajuS3=Mantel=Stand=Bres=Bres2=Log=rep(NA,N)

Lord1Stat=Lord2Stat=Lord3Stat=RajuU1Stat=RajuU2Stat=RajuU3Stat=RajuS1Stat=RajuS2Stat=RajuS3Stat=MantelStat=
  StandStat=BresStat=Bres2Stat=LogStat=rep(NA,N)

beta1pl=beta2pl=beta3pl=alpha2pl=alpha3pl=matrix(NA,2*m,N)

SDb1pl=SDb2pl=SDb3pl=cov2pl=cov3pl=SDa2pl=SDa3pl=matrix(NA,2*m,N)

gamma3pl=matrix(NA,m,N)

Z3pl_alt <-  read.csv2("simulace/ParamZeny.csv",sep = ',', stringsAsFactors = FALSE)
Z3pl_alt[, 2:7] <- lapply(Z3pl_alt[, 2:7], function(x) as.numeric(gsub(",", ".", x)))
View(param1_alt)
str(param1_alt)
param1_alt=cbind(as.numeric(Z3pl_alt[,2]),Z3pl_alt[,1],as.numeric(Z3pl_alt[,6]))
param2_alt<-param1_alt


View(Z3pl)
Z3pl<-as.matrix(read.csv2("simulace/ParamZeny.csv",sep = ',', stringsAsFactors = FALSE))
param1=cbind(Z3pl[,2],Z3pl[,1],Z3pl[,6])
param2<-param1

theta1=c(as.matrix(read.csv2("simulace\\theta11000.csv"), sep = ','))#resp. theta1200 nebo theta1500
theta2=c(as.matrix(read.csv2("simulace\\theta21000.csv"), sep = ','))#resp. theta2200
n1=length(theta1)#pocet ludi pre prvu skupinu
n2=length(theta2)#velkost druhej skupiny
skup=c(rep(0,n1),rep(1,n2))

dim(Z3pl)
#hladina
text=paste("Hladina",sep="")
H=simulace(N,param1_alt,param2_alt,text)

#sila-uniformni DIF
rozdil=c(0.5,1,2,4)
for(i in 1:4){
    text=paste("Sila",rozdil[i],sep="")
    param2[1,1]=param1[1,1]+rozdil[i]
    H=simulace(N,param1,param2,text)
}
#sila-neuniformni DIF
delta=c(0.4,0.6,0.8,1)
for(i in 1:5){
  param2[1,2]=2*param1[1,2]/(2+delta[i]*param1[1,2]/((1-param1[1,3])*log(2)))
  text=paste("NeuniD",delta[i],sep="")
  Neuni=simulace(N,param1,param2,text)
}

