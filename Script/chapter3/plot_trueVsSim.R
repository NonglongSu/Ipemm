#Usage true vs sim: Rscript --vanilla ../Script/chapter3/plot_trueVsSim.R k12 1
#Usage true vs im:  Rscript --vanilla ../Script/chapter3/plot_trueVsSim.R k12 0

suppressWarnings(suppressMessages(library(tidyverse)))
suppressPackageStartupMessages(library(stats))
suppressPackageStartupMessages(library(jsonlite))
suppressPackageStartupMessages(library(ggplot2))

#read file
readFile = function(inD, pat){
  Files= list.files(inD,pattern=pat,full.names=T)
  index = as.numeric(str_extract(basename(Files),'[^.]+')) #rectify the order
  Files = Files[order(index)]
  return(Files)
}

#without wz
doMat17 = function(Files,n){
  dMat = matrix(0,n,nvar)
  colnames(dMat)=c('A','C','G','T','s1','s2','s3','s4','s5','s6','omega','tau','ext.I','ext.D','r0*t', 'r1*t','r2*t')
  for (i in 1:n) {
    Jmp      = fromJSON(Files[i])
    go       = Jmp$gap.openning
    go.comb  = c(go[1]+go[4],go[2]+go[5],go[3]+go[6])
    rt       = -log(1-go.comb)
    ge       = Jmp$gap.extension
    dMat[i,] = c(Jmp$nuc.freq, Jmp$sigmas, Jmp$omega, Jmp$branch.length, ge, rt) 
  }
  return(dMat)
}

#with wz
doMat16 = function(Files,n){
  dMat  = matrix(0,n,nvar)
  colnames(dMat)=c('A','C','G','T','s1','s2','s3','s4','s5','s6','omega','tau','ext.I','ext.D','r*t','omegaZ')
  for (i in 1:n) {
    Jmp       = fromJSON(Files[i])
    g         = Jmp$gap.openning
    rt        = -log(1-2*g)/2
    ge        = Jmp$gap.extension
    dMat[i,]  = c(Jmp$nuc.freq, Jmp$sigmas, Jmp$omega, Jmp$branch.length, ge, rt, Jmp$omega.z) 
  }
  return(dMat)
}


#qqplot 
qqplt = function(tMat,dMat,p,q){
  colorblindP = "#CC6677"
  for(j in p:q){
    qqplot(tMat[,j],dMat[,j],xlab=colnames(tMat)[j],ylab=paste0(colnames(tMat)[j],'_est'),main=NULL)
    abline(0,1,col=colorblindP,lwd=2)
  }
}

#error perc
hiplt = function(diff,name,p,q){
  colorblindP = "#CC6677"
  xmax = round(max(diff[,p:q]),2)
  # dd   = density(diff[,j])
  # ymax = max(dd$y)
  for (j in p:q){
    hist(diff[,j],xlab=name[j],main=NULL,freq=FALSE,breaks=20,xlim=c(0,xmax))
    lines(density(sort(diff[,j])),col=colorblindP,lwd=2)
  }
}

####################################
readMatrix = function(inD){
  Files = readFile(inD,'est')
  n     = length(Files)
 
  if(nvar==17){
    res.mat = doMat17(Files,n)
  }else{
    res.mat = doMat16(Files,n)
  }
  return(res.mat)
}


##true vs sim/im plot
plot_true_vs = function(tMat,dMat,ouD){
  ouF1=paste0(ouD,"/","qqplot.pdf")
  ouF2=paste0(ouD,"/","error_perc.pdf")
  
  
  #qqplot
  pdf(ouF1,onefile=TRUE)
  
  ##sigmas
  par(mfrow=c(2,3),mar=c(5.1, 6.1, 4.1, 2.1))
  qqplt(tMat,dMat,5,10)
  
  ##omega,tau
  par(mfrow=c(1,2),mar=c(10.1, 6.1, 10.1, 2.1))
  qqplt(tMat,dMat,11,12)
  
  ##gap extensions
  par(mfrow=c(1,2))
  qqplt(tMat,dMat,13,14)
  
  if(nvar==16){
    par(mfrow=c(1,2))
    qqplt(tMat,dMat,15,16)
  }else{
    par(mfrow=c(1,3),mar=c(22.1, 6.1, 10, 2.1))
    qqplt(tMat,dMat,15,17)
  }
  dev.off()
  
  #Error percentage distribution plot
  name = colnames(tMat)
  diff = abs(dMat-tMat)/tMat
  print(colMeans(diff))
  
  pdf(ouF2,onefile=TRUE)
  par(mar=rep(4,4))
  
  #sigmas
  par(mfrow=c(2,3))
  hiplt(diff,name,5,10)
  
  #omega & tau
  par(mfrow=c(2,1))
  hiplt(diff,name,11,12)

  ##gap extensions
  par(mfrow=c(2,1))
  hiplt(diff,name,13,14)
  
  if(nvar==16){
    par(mfrow=c(2,1))
    hiplt(diff,name,15,16)
  }else{
    par(mfrow=c(3,1))
    hiplt(diff,name,15,17)
  }
  dev.off()
  
  # summary(diff[,1:17])
  # var(diff[,1:17])     
}



############################################################
#setwd("~/Dropbox (ASU)/Indel_project/chapter3")
#setwd("~/Dropbox (ASU)/Indel_project/chapter4")


#sim: num=1
#im:  num=0
#targ ='k12'
#targ = "NULL"

main = function(targ,num){
  n = as.numeric(num)
  
  #read true input
  inF  = "trueP.100.txt"
  tP   = read.table(inF,header=T)
  tMat = as.matrix(tP)
  nvar <<- ncol(tMat)
  
  #read simulated input
  if(targ!="NULL"){
    inD1 = paste0(targ,"/Results/Gse")
    inD2 = paste0(targ,"/Results/Pise")
    ouD1 = paste0(targ,"/Figure/500_sim")
    ouD2 = paste0(targ,"/Figure/500_im")
  }else{
    inD1 = "Results/Gse"
    inD2 = "Results/Pise"
    ouD1 = "Figure/500_sim"
    ouD2 = "Figure/500_im"
  }
  
  
  if(n==1){
    dMat    = readMatrix(inD1)
    plot_true_vs(tMat,dMat,ouD1)
  }else{
    dMat    = readMatrix(inD2)
    plot_true_vs(tMat,dMat,ouD2)
  }
 
}



########################################
args = commandArgs(trailingOnly=TRUE)
main(args[1],args[2])


