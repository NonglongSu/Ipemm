#ML estimate of gillespie simulations

suppressPackageStartupMessages(library(Matrix))
suppressPackageStartupMessages(library(stringr))
suppressPackageStartupMessages(library(seqinr))
suppressPackageStartupMessages(library(stringi))
suppressPackageStartupMessages(library(Biostrings))
suppressPackageStartupMessages(library(jsonlite))
suppressPackageStartupMessages(library(dfoptim))
suppressPackageStartupMessages(library(tidyverse))

#setwd("~/Dropbox (ASU)/Indel_project/chapter3")

######################################
#Count the nucleotide freq
count_freq = function(input){
  nuc.count = 0
  for (i in input){
    dna       = readDNAStringSet(i)
    nuc.count = nuc.count + oligonucleotideFrequency(dna,width=1)
  }
  nuc.freq  = colSums(nuc.count)/sum(nuc.count)
  if(sum(nuc.freq)==1){
    return(nuc.freq)
  }else{
    print("The nucleotide frequency sums up not 1!")
  }
}

# -Log-likelihood of
LL_min  = function(theta){
  nuc_codon_freq = init_f(cod,DD,64)
  fw  = nuc_codon_freq[[1]]
  cfw = nuc_codon_freq[[2]]
  
  rmat = GTR(theta[1:6], fw)
  Rmat = MG94(rmat,theta[7],cod,codonstrs,syn,64)
  -sum(log(expm(Rmat)*cfw)*DD)
}

######################################
main = function(inD,ouF){
  
  sub1 = "../Script/sources/"
  sub2 = "../Script/chapter3/"
  source(paste0(sub1,"codon_call.R"))
  source(paste0(sub1,"gtr.R"))
  source(paste0(sub1,"mg94.R"))
  source(paste0(sub1,"init_f.R"))
  source(paste0(sub1,"LL.R"))
  source(paste0(sub2,"phase_indel_prob3.R"))
  
  #construct codons and its degeneracy
  co.res    = codon_call()
  codons    = co.res[[1]]
  codonstrs = co.res[[2]]
  syn       = co.res[[3]]
  
  #nmkb method
  cod      <<- codons
  codonstrs<<- codonstrs
  syn      <<- syn
  
  #read input 
  #inD  = "Gs/3"
  Files = list.files(inD, full.names=T)
  index = as.numeric(str_extract(basename(Files),'[^.]+'))
  Files = Files[order(index)]
  K     = length(Files)
  Alist = list()
  for (i in 1:K) {Alist[[i]] = readBStringSet(Files[i])}
  
  #summary statistic
  Gm           = matrix(0,K,6)
  Mm           = matrix(0,K,6)
  codon_array  = array(0,c(64,64,K))
  llz          = rep(0,K)
  gl           = matrix(0,K,2)
  
  #Init parameters
  set.seed(8088)
  f0   = count_freq(Files)
  p0   = rep(0.1,7)
  g0   = rep(0.01,6)
  e0   = rep(0.4,2)
  
  rmat = GTR(p0[1:6],f0)
  t0   = -sum(diag(rmat)*f0) 
  
  #>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  iter = 1
  repeat{
    if(iter>20){
      print("Pass the 20 iterations limits!")
      break}
    
    #pre-cal avg gap size
    avg.gap = 1/(1-e0)
    if(any(avg.gap<1)){
      print("Average gap length less than 3!")
      break
    }
    
    #prellocate gtr/mg94 mat, codon freq.
    cf0   = sapply(seq(64), function(x){prod(f0[match(cod[x,],DNA_BASES)])})
    Rmat  = MG94(rmat,p0[7],cod,codonstrs,syn)
    Pmat  = log(expm(Rmat))
    
    
    for(i in 1:K){#E-step
      A      = Alist[[i]]
      res1   = ziqi_prob(A,g0,e0,Pmat,codonstrs,f0)
      Gm[i,] = res1[[1]]
      Mm[i,] = res1[[2]]
      codon_array[,,i] = res1[[3]]
      llz[i] = res1[[4]]
      gl[i,] = res1[[5]]
    }
    
    #M-step: parameter estimates
    
    #pseudo weight
    Wi = rep(1/K,K)
    
    ##gap openning
    nnew = colSums(Wi*Gm)         
    mnew = colSums(Wi*Mm) 
    gnew = nnew/(nnew+mnew) 
  
    
    ##gap extension
    w.avg.gap = c(mean(gl[which(gl[,1]>0),1]),mean(gl[which(gl[,2]>0),2]))
    enew      = 1-1/w.avg.gap
    
    ##codon
    datw = matrix(0,64,64)
    for (j in 1:K) {
      dat.wei = (1/K)*codon_array[,,j]
      datw    = datw + dat.wei
    }
    DD <<- datw
    pb = nmkb(fn=LL_min, par=p0, lower=0, upper=1.5, control=list(tol=1e-5,trace=F,maxfeval=5000))  #change the tolerance
    if(pb$convergence != 0){
      cat("Warning: failed convergence!")
    }else{
      pnew = pb$par
    }
    
    #cal. the tau
    nuc_f = init_f(cod,DD,64)    
    fnew  = nuc_f[[1]]
    rmat  = GTR(pnew[1:6],fnew)
    tnew  = -sum(diag(rmat)*fnew)
    
    print(gnew[1]+gnew[4])
    print(gnew[2]+gnew[5])
    print(gnew[3]+gnew[6])
    print(enew)
    print(w.avg.gap)
    print(pnew[1:6]/tnew)
    print(pnew[7])
    print(tnew)
    
    #rmse tolerance
    p    = c(g0,e0,p0)
    q    = c(gnew,enew,pnew)
    delta= (q-p)^2
    rmse = sqrt(mean(delta))
    cat(sprintf("iter:%i, rmse:%.6f\n",iter, rmse))
    if(rmse<=1e-4){
      break
    }else{
      iter=iter+1
      f0  = fnew
      p0  = pnew
      g0  = gnew
      e0  = enew
    }
  }
  
  ##output Jsons
  par_lst = list('nuc.freq'=fnew, 'sigmas'=pnew[1:6]/tnew, 'omega'=pnew[7], 'branch.length'=tnew, 
                 'gap.openning'=gnew,'gap.extension'=enew)
  
  sum_lst = list('gap.phases'=nnew,'avg.gap.size'=w.avg.gap)
  
  par_est = toJSON(par_lst)
  sum_stat= toJSON(sum_lst) 
  #ouF    = "Results/Gse/1.est.json"
  dname = dirname(ouF)
  fname = str_extract(basename(ouF),'[^.]+')
  ouF2  = paste0(dname,'/',fname,'.sum.json')
  write(par_est,ouF)
  write(sum_stat,ouF2)
}

###########################
args = commandArgs(trailingOnly=T)
main(args[1],args[2])