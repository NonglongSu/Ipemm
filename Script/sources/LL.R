#Log-likelihood
LL = function(cod,codonstrs,syn,theta,f0,cf0,dat,num=64){
  rmat  = GTR(theta[1:6],f0)
  Rmat  = MG94(rmat,theta[7],cod,codonstrs,syn,num)
  llt   = sum(log(expm(Rmat)*cf0)*dat)
  return(llt)
}