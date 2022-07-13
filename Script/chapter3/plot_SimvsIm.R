#Usage: Rscript --vanilla ../Script/chapter3/plot_alignVsSim.R  k12  16

suppressPackageStartupMessages(library(stats))
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(jsonlite))
suppressPackageStartupMessages(library(ggplot2))

#read file
readFile = function(inD, pat){
  Files= list.files(inD,pattern=pat,full.names=T)
  index = as.numeric(str_extract(basename(Files),'[^.]+')) #rectify the order
  Files = Files[order(index)]
  n    = length(Files)
  res  = list(Files,n)
  return(res)
}

#without wz
doMat17 = function(Files1,Files2,n){
  dMat1 = matrix(0,n,nvar)
  colnames(dMat1)=c('A','C','G','T','s1','s2','s3','s4','s5','s6','omega','tau','ext.I','ext.D','r0*t', 'r1*t','r2*t')
  for (i in 1:n) {
    Jmp      = fromJSON(Files1[i])
    go       = Jmp$gap.openning
    go.comb  = c(go[1]+go[4],go[2]+go[5],go[3]+go[6])
    rt       = -log(1-go.comb)
    ge       = Jmp$gap.extension
    dMat1[i,] = c(Jmp$nuc.freq, Jmp$sigmas, Jmp$omega, Jmp$branch.length, ge, rt) 
  }
  
  #summary statistics
  dMat2 = matrix(0,n,5)
  colnames(dMat2)=c('avg.gap.size.I','avg.gap.size.D','num.phase0','num.phase1','num.phase2')
  for (i in 1:n) {
    Jmp2      = fromJSON(Files2[i])
    numP      = Jmp2$gap.phases
    num.phase = c(numP[1]+numP[4],numP[2]+numP[5],numP[3]+numP[6])
    dMat2[i,] = c(Jmp2$avg.gap.size,num.phase)  
  }
  return(list(dMat1,dMat2))
}

#with wz
doMat16 = function(Files1,Files2,n){
  dMat1 = matrix(0,n,nvar)
  colnames(dMat1)=c('A','C','G','T','s1','s2','s3','s4','s5','s6','omega','tau','ext.I','ext.D','r*t','omegaZ')
  for (i in 1:n) {
    Jmp       = fromJSON(Files1[i])
    g         = Jmp$gap.openning
    rt        = -log(1-2*g)/2
    ge        = Jmp$gap.extension
    dMat1[i,] = c(Jmp$nuc.freq, Jmp$sigmas, Jmp$omega, Jmp$branch.length, ge, rt, Jmp$omega.z) 
  }
  
  #summary statistics
  dMat2 = matrix(0,n,8)
  colnames(dMat2)=c('avg.gap.size.I','avg.gap.size.D','num.zn0','num.zn1','num.zn2','num.zs0','num.zs1','num.zs2')
  for (i in 1:n) {
    Jmp2      = fromJSON(Files2[i])
    numN      = Jmp2$zn.phases
    num.zn    = c(numN[1]+numN[4],numN[2]+numN[5],numN[3]+numN[6])
    numS      = Jmp2$zs.phases
    num.zs    = c(numS[1]+numS[4],numS[2]+numS[5],numS[3]+numS[6])
    dMat2[i,] = c(Jmp2$avg.gap.size,num.zn,num.zs)  
  }
  return(list(dMat1,dMat2))
}

#qqplot
qqplt = function(tMat,dMat,p,q){
  colorblindP = "#CC6677"
  for(j in p:q){
    qqplot(tMat[,j],dMat[,j],xlab=colnames(tMat)[j],ylab=paste0(colnames(tMat)[j],'_est'),main=NULL)
    abline(0,1,col=colorblindP,lwd=2)
  }
}

#histogram
hiplt = function(diff,name,p,q){
  colorblindP = "#CC6677"
  xmax = round(max(diff[,p:q]),2)
  for (j in p:q){
    hist(diff[,j],xlab=name[j],main=NULL,freq=FALSE,breaks=20,xlim=c(0,xmax))
    lines(density(sort(diff[,j])),col=colorblindP,lwd=2)
  }
}

######################################
readMatrix = function(inD){
  Files1 = readFile(inD,'est')[[1]]
  Files2 = readFile(inD,'sum')[[1]]
  n      = readFile(inD,'est')[[2]]
  #parameter estimate
  
  if(nvar==17){
    res.mat = doMat17(Files1,Files2,n)
  }else{
    res.mat = doMat16(Files1,Files2,n)
  }
  
  return(res.mat)
}


##true align vs sim plot
plot_sim_im = function(dMat1,dMat2,ouD){
  ouF1=paste0(ouD,"/","SimvsIm.qqplot.pdf")
  ouF2=paste0(ouD,"/","SimvsIm.error_perc.pdf")
  
  #qqplot
  pdf(ouF1,onefile=TRUE)
  
  ##sigmas:par(mar=c(5.1,4.1,4.1,2.1))
  par(mar=rep(4,4))
  par(mfrow=c(2,3))
  qqplt(dMat1,dMat2,5,10)
  
  ##omega,tau
  par(mar=c(5.1,4.1,4.1,2.1))
  par(mfrow=c(2,1))
  qqplt(dMat1,dMat2,11,12)
  
  ##gap extensions
  par(mfrow=c(2,1))
  qqplt(dMat1,dMat2,13,14)
  
  if(nvar==16){
    par(mfrow=c(2,1))
    qqplt(dMat1,dMat2,15,16)
  }else{
    par(mfrow=c(3,1))
    qqplt(dMat1,dMat2,15,17)
  }
  dev.off()
  
  #Error percentage distribution plot
  name = colnames(dMat1)
  diff = abs(dMat2-dMat1)/dMat1
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

##summary stats
plot_sum_stats = function(dMat1.1,dMat2.1,ouD){
  ouF=paste0(ouD,"/","SimvsIm.sum.pdf")
  pdf(ouF,onefile=TRUE)
  #summary stat: num of gap phases, avergage gap sizes
  
  #qqplot
  ##extensions
  par(mfrow=c(2,1))
  qqplt(dMat1.1,dMat2.1,1,2)
  
  if(nvar==16){
    par(mfrow=c(2,3))
    qqplt(dMat1.1,dMat2.1,3,8)
  }else{
    par(mfrow=c(3,1))
    qqplt(dMat1.1,dMat2.1,3,5)
  }
  
  #error perc
  name = colnames(dMat1.1)
  diff = abs(dMat1.1-dMat2.1)/dMat1.1
  print(colMeans(diff))
  
  par(mfrow=c(2,1))
  hiplt(diff,name,1,2)
  
  if(nvar==16){
    par(mfrow=c(2,3))
    hiplt(diff,name,3,8)
  }else{
    par(mfrow=c(3,1))
    hiplt(diff,name,3,5)
  }
  dev.off()
}



#setwd("~/Dropbox (ASU)/Indel_project/chapter4")
#nv  = "16"
#targ="k12"
main = function(targ,nv){
  
  #read true alignment & simulated input
  nvar <<- as.numeric(nv)
  
  inD1 = paste0(targ,"/Results/Gse")
  inD2 = paste0(targ,'/Results/Pise')
  
  ouD = paste0(targ,"/Figure")
  
  dMat1 = readMatrix(inD1)[[1]]
  dMat2 = readMatrix(inD2)[[1]]
  
  dMat1.1 = readMatrix(inD1)[[2]]
  dMat2.1 = readMatrix(inD2)[[2]]
  
  
  plot_sim_im(dMat1,dMat2,ouD)
  plot_sum_stats(dMat1.1,dMat2.1,ouD)
  
}

###############################
args = commandArgs(trailingOnly=T)
main(args[1],args[2])
