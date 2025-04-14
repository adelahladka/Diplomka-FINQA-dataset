n1=1000
n2=1000
data2000=read.csv2(paste("hladina2000/Hladinam",m,"n",n1+n2,".csv",sep=""))
odhad2000=apply(data2000,1,function(x){sum(x<=0.05,na.rm=TRUE)/sum(!is.na(x))})
odhad2000=round(odhad2000,3)
#pocet NA pro IRT modely
nNA2000=apply(data2000,1,function(x){sum(is.na(x))})
konfint2000=rep(NA,length(odhad2000))
for(i in 1:length(odhad2000)){
  konfint2000[i]=paste("(",round(binom.test(round((1000-nNA2000[i])*odhad2000[i]),n=(1000-nNA2000[i]))$conf.int[1],3),", ",
                     round(binom.test(round((1000-nNA2000[i])*odhad2000[i]),n=(1000-nNA2000[i]))$conf.int[2],3),")",sep="")                     
}

n1=500
n2=1000
data1500=read.csv2(paste("hladina1500/Hladinam",m,"n",n1+n2,".csv",sep=""))
odhad1500=apply(data1500,1,function(x){sum(x<=0.05,na.rm=TRUE)/sum(!is.na(x))})
odhad1500=round(odhad1500,3)
#pocet NA pro IRT modely
nNA1500=apply(data1500,1,function(x){sum(is.na(x))})
konfint1500=rep(NA,length(odhad1500))
for(i in 1:length(odhad1500)){
  konfint1500[i]=paste("(",round(binom.test(round((1000-nNA1500[i])*odhad1500[i]),n=(1000-nNA1500[i]))$conf.int[1],3),", ",
                       round(binom.test(round((1000-nNA1500[i])*odhad1500[i]),n=(1000-nNA1500[i]))$conf.int[2],3),")",sep="")                     
}



n1=200
n2=200
data400=read.csv2(paste("hladina400/Hladinam",m,"n",n1+n2,".csv",sep=""))
odhad400=apply(data400,1,function(x){sum(x<=0.05,na.rm=TRUE)/sum(!is.na(x))})
odhad400=round(odhad400,3)
#pocet NA pro IRT modely
nNA400=apply(data400,1,function(x){sum(is.na(x))})
konfint400=rep(NA,length(odhad400))
for(i in 1:length(odhad400)){
  konfint400[i]=paste("(",round(binom.test(round((1000-nNA400[i])*odhad400[i]),n=(1000-nNA400[i]))$conf.int[1],3),", ",
                       round(binom.test(round((1000-nNA400[i])*odhad400[i]),n=(1000-nNA400[i]))$conf.int[2],3),")",sep="")                     
}

odhad=cbind(odhad400,odhad1500,odhad2000)
nNA=cbind(nNA400,nNA1500,nNA2000)
konfint=cbind(konfint400,konfint1500,konfint2000)
write.csv2(odhad,"OdhadHladiny/hladina.csv")
write.csv2(nNA,"OdhadHladiny/nNA.csv")
write.csv2(konfint,"OdhadHladiny/konfint.csv")
