###################Odhad sily############################################
#Zacneme neuniformnim DIF

n1=200
n2=200
data1=read.csv2(paste("neun400/NeuniD0.4m",m,"n",n1+n2,".csv",sep=""))
data2=read.csv2(paste("neun400/NeuniD0.6m",m,"n",n1+n2,".csv",sep=""))
data3=read.csv2(paste("neun400/NeuniD0.8m",m,"n",n1+n2,".csv",sep=""))
data4=read.csv2(paste("neun400/NeuniD1m",m,"n",n1+n2,".csv",sep=""))

sila400=matrix(NA,13,4)
sila400[,1]=apply(data1,1,function(x){mean(x<=0.05,na.rm=TRUE)})
sila400[,2]=apply(data2,1,function(x){mean(x<=0.05,na.rm=TRUE)})
sila400[,3]=apply(data3,1,function(x){mean(x<=0.05,na.rm=TRUE)})
sila400[,4]=apply(data4,1,function(x){mean(x<=0.05,na.rm=TRUE)})
write.csv2(sila400,"OdhadSily/SilaNeun400.csv")

NA400=matrix(NA,13,4)
(NA400[,1]=apply(data1,1,function(x){sum(is.na(x))}))
(NA400[,2]=apply(data2,1,function(x){sum(is.na(x))}))
(NA400[,3]=apply(data3,1,function(x){sum(is.na(x))}))
(NA400[,4]=apply(data4,1,function(x){sum(is.na(x))}))
write.csv2(NA400,"OdhadSily/NaNeun400.csv")
#1500
n1=500
n2=1000

data1=read.csv2(paste("neun1500/NeuniD0.4m",m,"n",n1+n2,".csv",sep=""))
data2=read.csv2(paste("neun1500/NeuniD0.6m",m,"n",n1+n2,".csv",sep=""))
data3=read.csv2(paste("neun1500/NeuniD0.8m",m,"n",n1+n2,".csv",sep=""))
data4=read.csv2(paste("neun1500/NeuniD1m",m,"n",n1+n2,".csv",sep=""))

sila1500=matrix(NA,13,4)
sila1500[,1]=apply(data1,1,function(x){mean(x<=0.05,na.rm=TRUE)})
sila1500[,2]=apply(data2,1,function(x){mean(x<=0.05,na.rm=TRUE)})
sila1500[,3]=apply(data3,1,function(x){mean(x<=0.05,na.rm=TRUE)})
sila1500[,4]=apply(data4,1,function(x){mean(x<=0.05,na.rm=TRUE)})
write.csv2(sila1500,"OdhadSily/SilaNeun1500.csv")

NA1500=matrix(NA,13,4)
(NA1500[,1]=apply(data1,1,function(x){sum(is.na(x))}))
(NA1500[,2]=apply(data2,1,function(x){sum(is.na(x))}))
(NA1500[,3]=apply(data3,1,function(x){sum(is.na(x))}))
(NA1500[,4]=apply(data4,1,function(x){sum(is.na(x))}))

write.csv2(NA1500,"OdhadSily/NaNeun1500.csv")
#2000
n1=1000
n2=1000

data1=read.csv2(paste("neun2000/NeuniD0.4m",m,"n",n1+n2,".csv",sep=""))
data2=read.csv2(paste("neun2000/NeuniD0.6m",m,"n",n1+n2,".csv",sep=""))
data3=read.csv2(paste("neun2000/NeuniD0.8m",m,"n",n1+n2,".csv",sep=""))
data4=read.csv2(paste("neun2000/NeuniD1m",m,"n",n1+n2,".csv",sep=""))

sila2000=matrix(NA,13,4)
sila2000[,1]=apply(data1,1,function(x){mean(x<=0.05,na.rm=TRUE)})
sila2000[,2]=apply(data2,1,function(x){mean(x<=0.05,na.rm=TRUE)})
sila2000[,3]=apply(data3,1,function(x){mean(x<=0.05,na.rm=TRUE)})
sila2000[,4]=apply(data4,1,function(x){mean(x<=0.05,na.rm=TRUE)})
write.csv2(sila2000,"OdhadSily/SilaNeun2000.csv")

NA2000=matrix(NA,13,4)
(NA2000[,1]=apply(data1,1,function(x){sum(is.na(x))}))
(NA2000[,2]=apply(data2,1,function(x){sum(is.na(x))}))
(NA2000[,3]=apply(data3,1,function(x){sum(is.na(x))}))
(NA2000[,4]=apply(data4,1,function(x){sum(is.na(x))}))
write.csv2(NA2000,"OdhadSily/NaNeun2000.csv")
########################################################################
#Uniformne DIF
n1=200
n2=200
data1=read.csv2(paste("unif400/Sila0.5m",m,"n",n1+n2,".csv",sep=""))
data2=read.csv2(paste("unif400/Sila1m",m,"n",n1+n2,".csv",sep=""))
data3=read.csv2(paste("unif400/Sila2m",m,"n",n1+n2,".csv",sep=""))
data4=read.csv2(paste("unif400/Sila4m",m,"n",n1+n2,".csv",sep=""))

sila400U=matrix(NA,13,4)
sila400U[,1]=apply(data1,1,function(x){mean(x<=0.05,na.rm=TRUE)})
sila400U[,2]=apply(data2,1,function(x){mean(x<=0.05,na.rm=TRUE)})
sila400U[,3]=apply(data3,1,function(x){mean(x<=0.05,na.rm=TRUE)})
sila400U[,4]=apply(data4,1,function(x){mean(x<=0.05,na.rm=TRUE)})
write.csv2(sila400U,"OdhadSily/SilaUn400.csv")

NA400U=matrix(NA,13,4)
(NA400U[,1]=apply(data1,1,function(x){sum(is.na(x))}))
(NA400U[,2]=apply(data2,1,function(x){sum(is.na(x))}))
(NA400U[,3]=apply(data3,1,function(x){sum(is.na(x))}))
(NA400U[,4]=apply(data4,1,function(x){sum(is.na(x))}))
write.csv2(NA400U,"OdhadSily/NaUn400.csv")

#1500
n1=500
n2=1000
data1=read.csv2(paste("unif1500/Sila0.5m",m,"n",n1+n2,".csv",sep=""))
data2=read.csv2(paste("unif1500/Sila1m",m,"n",n1+n2,".csv",sep=""))
data3=read.csv2(paste("unif1500/Sila2m",m,"n",n1+n2,".csv",sep=""))
data4=read.csv2(paste("unif1500/Sila4m",m,"n",n1+n2,".csv",sep=""))

sila1500U=matrix(NA,13,4)
sila1500U[,1]=apply(data1,1,function(x){mean(x<=0.05,na.rm=TRUE)})
sila1500U[,2]=apply(data2,1,function(x){mean(x<=0.05,na.rm=TRUE)})
sila1500U[,3]=apply(data3,1,function(x){mean(x<=0.05,na.rm=TRUE)})
sila1500U[,4]=apply(data4,1,function(x){mean(x<=0.05,na.rm=TRUE)})
write.csv2(sila1500U,"OdhadSily/SilaUn1500.csv")
NA1500U=matrix(NA,13,4)
(NA1500U[,1]=apply(data1,1,function(x){sum(is.na(x))}))
(NA1500U[,2]=apply(data2,1,function(x){sum(is.na(x))}))
(NA1500U[,3]=apply(data3,1,function(x){sum(is.na(x))}))
(NA1500U[,4]=apply(data4,1,function(x){sum(is.na(x))}))
write.csv2(NA1500U,"OdhadSily/NaUn1500.csv")

#2000
n1=1000
n2=1000
data1=read.csv2(paste("unif2000/Sila0.5m",m,"n",n1+n2,".csv",sep=""))
data2=read.csv2(paste("unif2000/Sila1m",m,"n",n1+n2,".csv",sep=""))
data3=read.csv2(paste("unif2000/Sila2m",m,"n",n1+n2,".csv",sep=""))
data4=read.csv2(paste("unif2000/Sila4m",m,"n",n1+n2,".csv",sep=""))

sila2000U=matrix(NA,13,4)
sila2000U[,1]=apply(data1,1,function(x){mean(x<=0.05,na.rm=TRUE)})
sila2000U[,2]=apply(data2,1,function(x){mean(x<=0.05,na.rm=TRUE)})
sila2000U[,3]=apply(data3,1,function(x){mean(x<=0.05,na.rm=TRUE)})
sila2000U[,4]=apply(data4,1,function(x){mean(x<=0.05,na.rm=TRUE)})
sila2000U=round(sila2000U,3)
write.csv2(sila2000U,"OdhadSily/SilaUn2000.csv")

NA2000U=matrix(NA,13,4)
(NA2000U[,1]=apply(data1,1,function(x){sum(is.na(x))}))
(NA2000U[,2]=apply(data2,1,function(x){sum(is.na(x))}))
(NA2000U[,3]=apply(data3,1,function(x){sum(is.na(x))}))
(NA2000U[,4]=apply(data4,1,function(x){sum(is.na(x))}))
write.csv2(NA2000U,"OdhadSily/NaUn2000.csv")
