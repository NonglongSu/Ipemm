#1: output the prob of ziqi's model
#2: output the summuary stats

#Cal. the inserted freq prob
prob_insert = function(As2,g1,f0){
  if(length(g1)==0){
    return (0)
  }else{
      pos1 = start(g1)
      pos2 = end(g1)
      insN = c()
      for(j in 1:length(g1)){
        insN = c(insN,As2[pos1[j]:pos2[j]])
      }
  }
  f.rank = match(insN,DNA_BASES)
  f.prob = sum(log(f0[f.rank]))
  return(f.prob)
}


# 0 1 1; 0 1 0; 0 1 1; 0 2 0; 0 0 1; 0 0 2
#Count number of gap edges.
#Format: ins: del: 
count_gap = function(As,g){
  phase = matrix(0,3,2)
  lG = c(0,0)
  lT = c(0,0)
  if(length(unlist(g)) == 0){
    return(list(phase,lG,lT))
  }else{#separate the ins from del
    for (j in seq(2)) {
      pos  = start(g[[j]])
      wid  = width(g[[j]])
      
      rem = pos %% 3
      phase[1,j] = length(which(rem == 1))
      phase[2,j] = length(which(rem == 2))
      phase[3,j] = length(which(rem == 0))
    }
  }
  #gap number; gap length
  lG = c(length(g[[1]]),length(g[[2]]))
  lT = c(sum(width(g[[1]])), sum(width(g[[2]])))/3
  
  res = list(phase,lG,lT)
  return(res)
}

# 4 3 3; 6 5 6; 4 3 3; 4 2 4; 6 6 5; 4 4 2
#1--- AAA AAA    2AAA --- AAA      3---AAA
#1AAA --- AAA    2AAA AAA ---      3AAAAAA

#Count number of no-gap edges
count_nogap = function(As,g,gCount){
  lenA  = length(As[[1]])
  ug    = unlist(g)
  lenT  = sum(width(ug))
  
  last0 = 0 #last site
  if(lenA %in% end(ug)){
    lastg0 = ug[which(end(ug) %in% lenA)]
    last0  = last0 + 1
    if((start(lastg0)-1) %in% end(ug)){#2
      seclastg0 = ug[which((start(lastg0)-1) %in% end(ug))]
      last0     = last0 + 1
    }
  }
  
  endflag = 0 #ins-1, del-2, match-0
  if(last0==1){
    end1 = end(g)[[1]]
    if(end(lastg0) %in% end1){#ins>end
      gCount[1,1] = gCount[1,1]-1
      endflag = 1
    }else{#del>end
      gCount[1,2] = gCount[1,2]-1
      endflag = 2
    }
  }else if(last0==2){#ins>del>end
    gCount[1,] = gCount[1,]-1
    endflag    = 2
  }
  
  exp.edge = (lenA-lenT)/3
  M        = matrix(0,3,2)
  for (j in 1:3) {
    M[j,]   = exp.edge - gCount[j,]
  }
  res = list(M,endflag)
  return(res)
}


#Remove all gap-positioned string
rmgap = function(A,As,g){
  g    = unlist(g)
  rm.A = stri_sub_replace_all(A,from=sort(start(g)),to=sort(end(g)),replacement='')
  rm.A = DNAStringSet(rm.A)
  if(all(width(rm.A)%%3 == 0)){
    return(rm.A)
  }else{
    print("Warning:removed-gap sequences are not multiple of three")
    break
  }
}

#Generate obs-codon matrix
countN = function(rA, codonstrs, ncd){
  seqs = str_split(rA,'')
  len  = length(seqs[[1]])
  nmat = matrix(0,ncd,ncd)
  i=1
  while(i<len) {
    c1 = paste0(seqs[[1]][i:(i+2)], collapse = '')
    c2 = paste0(seqs[[2]][i:(i+2)], collapse = '')
    coor1 = which(codonstrs %in% c1)
    coor2 = which(codonstrs %in% c2)
    nmat[coor1,coor2] = nmat[coor1,coor2] + 1
    i=i+3
  }
  
  return(nmat)
}




ziqi_prob = function(A,g0,e3,Pmat,codonstrs,f0,ncd){
  
  As = str_split(A,'')
  g  = IRangesList(lapply(As, function(x){IRanges(x=='-')}))
  
  #cal the inserted freq prob
  insP = prob_insert(As[[2]],g[[1]],f0)
  
  
  #count number of gap edges
  Cg    = count_gap(As,g)
  N.012 = Cg[[1]]
  num.g = Cg[[2]]
  len.g = Cg[[3]]
  
  Cng   = count_nogap(As,g,N.012)
  M.012 = Cng[[1]]
  Eflag = Cng[[2]]
  
  
  ##end state prob
  if(Eflag==0){#match
    endP = log(1-g0[1])
    endP = unname(endP)
  }else if(Eflag==1){#ins
    endP = log(1-e3[1])
  }else{#del
    endP = log(1)
  }
  
  #extension prob of ins/del
  avg.g = len.g/num.g
  if(all(len.g==0)){
    avg.g = c(0,0)
  }else if(len.g[1]==0){
    avg.g[1] = 0
  }else if(len.g[2]==0){
    avg.g[2] = 0
  }
  scoreE = (len.g[1]-num.g[1])*log(e3[1])+num.g[1]*log(1-e3[1]) + (len.g[2]-num.g[2])*log(e3[2])+num.g[2]*log(1-e3[2])
  
  
  ##transition prob
  scoreT = log(prod(g0^N.012)) + sum(log((1-g0)^M.012))
  #print(scoreT+scoreE+endP)
  
  rA     = rmgap(A,As,g)
  dat    = countN(rA,codonstrs,ncd)
  scoreP = sum(Pmat*dat)
  
  ##sum LL
  score_ziqi = scoreE + scoreT + endP + scoreP + insP
  if(score_ziqi==-Inf){
    print("minus infinity z score!")
  }
  #summary stat
  res = list(c(N.012),c(M.012),dat,score_ziqi,avg.g)
  
  #print(score_ziqi)
  return(res)
}
