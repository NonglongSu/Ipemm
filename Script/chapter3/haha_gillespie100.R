suppressPackageStartupMessages(library(Biostrings))
suppressPackageStartupMessages(library(stringr))
suppressPackageStartupMessages(library(seqinr))
suppressPackageStartupMessages(library(plyr))
suppressPackageStartupMessages(library(R.utils))

#setwd("~/Dropbox (ASU)/Indel_project/chapter3")

#deletion starts to the right(pos included)
sim_del = function(s, dnaB, pos, posB, k){
  po = pos %% 3
  flag = 0
  if(po == 2){#phase 1
    ref    = paste0(s[pos-1],  s[pos+k], s[pos+k+1],collapse="")
    sub1   = paste0(s[pos-1],  s[pos],   s[pos+1],  collapse="")
    sub2   = paste0(s[pos+k-1],s[pos+k], s[pos+k+1],collapse="")
    sec    = syn[[ref]]
    if(sub1 %in% sec || sub2 %in% sec){
      flag = 1
    }
  }else if(po == 0){#phase2
    ref    = paste0(s[pos-2],  s[pos-1],  s[pos+k],collapse="")
    sub1   = paste0(s[pos-2],  s[pos-1],  s[pos],  collapse="")
    sub2   = paste0(s[pos+k-2],s[pos+k-1],s[pos+k],collapse="")
    sec    = syn[[ref]]
    if(sub1 %in% sec || sub2 %in% sec){
      flag = 1
    }
  }else{#phase0
    flag = 1
  }
  
  #AAA A-- -AA
  #AAA ANN NAA
  #--- --- NAA
  if((flag==1) || (runif(1L)<=Wz)){#selection check
    i=0
    while(i<k){
      #AAA --- A-- -AA
      #AA- --- --- --A
      if(dnaB[posB+i]=='-'){
        k=k+1
      }else{
        dnaB[posB+i]='-'
      }
      i=i+1
    }
  }
  
  return(dnaB)
}

#substition starts on the choosing position.
sim_sub = function(s, dnaB, pos, posB, gtr){
  mu = sample(DNA_BASES, 1L, prob = gtr[which(DNA_BASES == s[pos]),])  #If nuc A, it must change to C,G,T.
  
  po = pos %% 3
  flag = 0
  if(po == 2){#phase 1
    flag  = 0
  }else if(po == 0){#phase2
    ref = paste0(s[(pos-2):pos],collapse="")
    sub = paste0(s[pos-2],s[pos-1],mu,collapse="")
    sec = syn[[ref]]
    if(sub %in% sec){
      flag = 1
    }
  }else{#phase0
    ref = paste0(s[pos:(pos+2)],collapse="")
    sub = paste0(mu, s[pos+1],s[pos+2],collapse="")
    sec = syn[[ref]]
    if(sub %in% sec){
      flag = 1
    }
  }
  
  if((flag==1) || (runif(1L)<W)){
    dnaB[posB] = mu
  }
  return(dnaB)
}

#insertion starts to the left(including immortal link)
sim_ins = function(s, dnaB, pos, posB, k){
  ins = sample(DNA_BASES,k,prob=Pi,replace=TRUE)  #throw a k-sided dice for inserted nucs.
  
  #type-N:0, type-S:1
  po   = pos %% 3
  flag = 0
  
  if(po == 2){#phase2
    ref    = paste0(s[pos-1], s[pos], s[pos+1], collapse="")
    sub1   = paste0(s[pos-1], s[pos], ins[1],   collapse="")
    sub2   = paste0(ins[k-1], ins[k], s[pos+1], collapse="")
    sec    = syn[[ref]]
    if(sub1 %in% sec || sub2 %in% sec){
      flag = 1
    }
  }else if(po == 1){#phase1
    ref    = paste0(s[pos],  s[pos+1], s[pos+2], collapse = "")
    sub1   = paste0(s[pos],  ins[1],   ins[2],   collapse = "")
    sub2   = paste0(ins[k],  s[pos+1], s[pos+2], collapse = "")
    sec    = syn[[ref]]
    if(sub1 %in% sec || sub2 %in% sec){
      flag = 1
    }
  }else{#phase0
    flag = 1
  }
  
  #selection check
  if( (flag==1) || (runif(1L)<=Wz) ){
    newB = append(dnaB,ins,posB)                     
    stamp = 1
  }else{
    newB  = dnaB
    stamp = 0
  }
  res = list(newB,stamp)
  return(res)
}

#sequential simulation with same brlen#######
D_onestep = function(r1,ext,dna0){
  dnaB = dna0
  tau  = 0
  #total rates
  phase   = (seq_along(dnaB)+2) %% 3 + 1
  phase.r = cbind(r1[phase],rep(0,length(phase)))
  
  for (iter in 1:500) {
    #no events should occur again on the deletion box
    phase.r[which(dnaB=='-')] = 0
    
    #sample position and mutation types
    Pos  = sample(nrow(phase.r), 1L, prob=phase.r[,1])
    
    #Poisson waiting time
    total_rate = sum(phase.r[,1])
    tau = tau + rexp(1,total_rate)
    if((total_rate==0) || (tau>=brlen)){
      break
    }
    
    #count number of gaps before current position
    ngap = length(which(phase.r[(1:Pos),1]==0))
    
    #update dna
    dna  = dnaB[which(dnaB!='-')]
    pos  = Pos-ngap
    
    #draw geometric dist of indel length  
    k = rgeom(1,1-ext)
    k = 3*(k+1)        
    if(length(dna)-pos<k){#delete the end 
      next
    }
    
    dnaB = sim_del(dna, dnaB, pos, Pos, k)
  }
  cat(sprintf("iter:%d  tau:%.6f\n", iter,tau))
  
  if((length(dnaB)!=length(dna0)) || (length(dnaB) %%3 !=0)){
    print("Warning:after deletion, alignment length should be the same!")
  }
  return(dnaB)
}

S_onestep = function(r2,gtr,dnaB){
  tau = 0
  #total rates
  phase = (seq_along(dnaB)+2) %% 3 + 1
  
  for (iter in 1:1000) {
    #update the mu
    phase.r = cbind(r2[dnaB],rep(0,length(phase)))
    
    #sample position and mutation types
    Pos  = sample(nrow(phase.r), 1L, prob=phase.r[,1])
    
    #Poisson waiting time
    total_rate = sum(phase.r[,1])
    tau = tau + rexp(1,total_rate)
    if((total_rate==0) || (tau>=brlen)){
      break
    }
    
    #count number of gaps before current position
    ngap = length(which(phase.r[(1:Pos),1]==0))
    
    #update dna
    dna  = dnaB[which(dnaB!='-')]
    pos  = Pos-ngap
    dnaB = sim_sub(dna, dnaB, pos, Pos, gtr)
  }
  
  cat(sprintf("iter:%d  tau:%.6f\n", iter,tau))
  return(dnaB)
}

I_onestep= function(r1,ext,dnaB1){
  tau = 0
  
  #total rates
  phase.r = matrix(0,length(dnaB1)+1,2)
  phase   = (seq_along(dnaB1)) %% 3 + 1
  phase.r[1,1] = r1[1]
  phase.r[2:nrow(phase.r),] = cbind(r1[phase],rep(0,length(phase)))
  #no events should occur again on the deletion box
  phase.r[which(dnaB1=='-')+1] = 0
  
  
  for (iter in 1:500) {
    #sample position and mutation types
    Pos  = sample(nrow(phase.r), 1L, prob=phase.r[,1])
    
    #Poisson waiting time
    total_rate = sum(phase.r[,1])
    tau = tau + rexp(1,total_rate)
    if((total_rate==0) || (tau>=brlen)){
      break
    }
    
    #draw geometric dist of indel length
    k = rgeom(1,1-ext)
    k = 3*(k+1)
    
    #count number of gaps before current position
    posB = Pos - 1
    if(posB == 0){#immortal link
      ins      = sample(DNA_BASES,k,prob=Pi,replace=TRUE)
      dnaB1    = append(dnaB1,ins,posB) 
      Isinsert = 1
    }else{
      ngap     = length(which(phase.r[(1:posB),1]==0))
      #update dna
      dna      = dnaB1[which(dnaB1!='-')]
      pos      = posB-ngap
      dnaBs    = sim_ins(dna, dnaB1, pos, posB, k)
      dnaB1    = dnaBs[[1]]
      Isinsert = dnaBs[[2]]
    }
    
    if(Isinsert==1){
      Bi = append(phase.r[,2],rep(1,k),Pos)     
      if(Pos %% 3 == 1){#phase0
        Rate = append(phase.r[,1],rep(r1[c(2,3,1)],k/3),Pos)
      }else if(Pos %% 3 == 2){#phase1
        Rate = append(phase.r[,1],rep(r1[c(3,1,2)],k/3),Pos)
      }else{#phase2
        Rate = append(phase.r[,1],rep(r1,k/3),Pos)
      }
      phase.r = cbind(Rate,Bi)
    }
  }
  
  cat(sprintf("iter:%d  tau:%.6f\n", iter,tau))
  if(nrow(phase.r)!=length(dnaB1)+1){
    print("Warning: immortal link does not exist!")
  }
  res = list(dnaB1,phase.r)
  return(res)
}

#Align it back
align_back = function(dna0,dnaB2){
  A   = dna0
  A   = insert(A,1,'#')
  B   = dnaB2[[1]]
  bi  = dnaB2[[2]][,2]
  
  bipos = which(bi==1)
  if(length(bipos)!=0){
    for(i in 1:length(bipos)){
      A = append(A,'-',bipos[i]-1)  
    }
  }
  A = A[-1]
  if(length(A)!=length(B)){
    print("Warning: someting wrong with the back alignment!")
  }
  
  DnaA  = paste0(A,collapse="")
  DnaB  = paste0(B,collapse="")
  Align = BStringSet(c(DnaA,DnaB))
  return(Align)
}



#setwd("~/Dropbox (ASU)/Indel_project/chapter3")
#setwd("~/Dropbox (ASU)/Indel_project/chapter4")
#######################################################
main = function(ouD,inF,l,omega_z,ss){
  
  #read input
  #ouD = "Gs/98"
  #inF = "trueP.100.txt"
  cmd = paste0("mkdir -p ",ouD, sep=' ')
  system(cmd)
  
  tag = as.numeric(str_extract(basename(ouD), "[^.]+"))
  tP  = read.table(inF,header=T,sep='') 
  tp  = unlist(tP[tag,])
  
  Pi    <<- tp[1:4]
  Sigma <<- tp[5:10]
  W     <<- tp[11]
  brlen <<- tp[12]
  ext    = tp[13:14]
  r1     = tp[15:17]/(2*brlen)
  
  
  
  
  omegaz = as.numeric(omega_z)  #znzs: omega_z=1 (default)
  Wz     <<- 1
  
  
  # construct codons and its degeneracy
  codons = cbind(rep(DNA_BASES, each = 16),
                 rep(DNA_BASES, times = 4, each = 4),
                 rep(DNA_BASES, 16))
  codonstrs  = apply(codons, 1, stringr::str_c, collapse = "")
  syn        = syncodons(codonstrs)
  names(syn) = toupper(names(syn))
  syn        <<- lapply(syn, toupper)
  
  #construct GTR model
  #Sigma = c(0.7135135,2.9447236,0.3246753,1.1608040,2.9123377,0.6428571)
  #Sigma  = 1:6
  r      = matrix(0,4,4)
  r[lower.tri(r)] = Sigma
  r       = r + t(r)
  r       = t(r*Pi)
  diag(r) = -rowSums(r)
  T   = -sum(diag(r)*Pi)
  gtr = r / T
  print(-sum(diag(gtr)*Pi))

  r2        = -diag(gtr)
  r2[5]     = 0
  names(r2) = c(DNA_BASES,'-')
  diag(gtr) = 0                #avoid negative prob.
  
  
  #########################################PART II Run simulation
  set.seed(8088)
  len   = as.numeric(l)
  dna0  = sample(DNA_BASES, len, prob=Pi, replace=TRUE)
  dna00 = paste0(dna0,collapse='')
  
  sim.IDS = list()
  ssize   = as.numeric(ss)
  for(i in 1:ssize){
    dnaB  = D_onestep(r1,ext[2],dna0)
    dnaB1 = S_onestep(r2,gtr,dnaB)
    dnaB2 = I_onestep(r1,ext[1],dnaB1)
    Align = align_back(dna0,dnaB2)
    
    names(Align) = c('Seq1','Seq2')
    writeXStringSet(Align, paste0(ouD,'/',i,".fa"))
    sim.IDS[[i]] = Align
  }
  
}


########################################
args = commandArgs(trailingOnly=T)
main(args[1],args[2],args[3],args[4],args[5])