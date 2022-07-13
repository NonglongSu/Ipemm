#EM method obtain the est of init nuc freq.

init_f = function(cod,dat,num){
  #empirical f
  f1     = rowSums(dat) + colSums(dat)  
  f1.id  = sapply(seq(num), function(x){match(cod[x, ], DNA_BASES)})
  base.f = c()
  for (i in seq(4)) {
    base.id = sapply(seq(num), function(x){length(which(f1.id[, x] == i))})
    base.f  = c(base.f, sum(base.id*f1))
  }
  f  = base.f/sum(base.f)
  if(round(sum(f),3)!=1){print("Warning:sum of nuc freq is not 1!")}
  cf = sapply(seq(num), function(x){prod(f[match(cod[x, ], DNA_BASES)])})
  
  if(num==61){
    f0=f
    sum_61 = sum(f1)
    for (i in 1:20) {
      Pi_stp = c(f0[1]^2*f0[4],f0[1]*f0[3]*f0[4],f0[3]*f0[1]*f0[4]) 
      Pi_61  = 1 - sum(Pi_stp)
      
      sum_stp  = sum_61/Pi_61 - sum_61      
      stp_norm = Pi_stp/sum(Pi_stp) *sum_stp
      
      cf0   = sapply(seq(num), function(x){prod(f0[match(cod[x, ], DNA_BASES)])})
      ll    = sum(f1*log(cf0))- sum_61*log(Pi_61)                         
      cat(sprintf("%i: %.6f\n", i, ll))
      
      denom = 3*(sum_61+sum_stp)
      f0    = c(base.f[1]+2*stp_norm[1]+stp_norm[2]+stp_norm[3], 
                base.f[2],
                base.f[3]+stp_norm[2]+stp_norm[3],
                base.f[4]+sum_stp)/denom
    }
    cf0 = cf0/sum(cf0)
    res = list(f0,cf0)
  }else{
    if(round(sum(cf),3)!=1){print("Warning:sum of codon freq is not 1!")}
    res = list(f,cf)
  }
  

  return(res)
}

