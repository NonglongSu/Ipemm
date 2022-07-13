#init sigma
#create nuc transition matrix [4-fold-degeneracy-codons only]
init_sigma = function(syn,dat){
  fourD_id = which(sapply(syn, function(x){length(x) == 4}))
  k=1
  ndat = matrix(0,4,4)
  while (k<length(fourD_id)) {
    i = j = fourD_id[k:(k+3)]
    ndat  = ndat+dat[i,j]
    k = k + 4
  }
  fn = rowSums(ndat) + colSums(ndat)        #neutral freq
  fn = fn/sum(fn)
  #print(sum(fn))
  
  #symmetric average of dat
  sdat  = matrix(0,4,4)
  diag(sdat) = diag(ndat)
  sdat  = (ndat + t(ndat))/2 
  sf    = colSums(sdat)/sum(sdat)
  #print(sum(sf))
  Dhat  = diag(sf)
  phat  = t(sdat/colSums(sdat))
  #logm(phat)
  
  eigP  = eigen(phat)
  Athat = eigP$vectors %*% diag(log(eigP$values)) %*%  inv(eigP$vectors)
  that  = -sum(diag(Athat*Dhat))
  Ahat  = t(Athat)/that
  #print(that)  
  #print(Ahat)
  s     = Ahat[lower.tri(t(Ahat))]
  
  if(any(s<0)){#HKY85 
    mean(diag(Ahat))
    i = c(1,2,3,4)
    j = c(3,4,1,2)
    ts = mean(Ahat[cbind(i,j)])   #transitions
    diag(Ahat)=0 
    tv = (sum(Ahat)-4*ts)/8       #transversions
    
    hky = matrix(0,4,4)
    hky[lower.tri(hky)] = c(tv,ts,tv,tv,ts,tv)
    hky = hky + t(hky)
    hky = t(hky*sf)
    diag(hky) = -rowSums(hky)
    s   = hky[lower.tri(t(hky))]
  }
  return(s)
}
