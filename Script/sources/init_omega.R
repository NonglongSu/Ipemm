#find the init omega
init_omega = function(cod,codonstrs,syn,dat,f0,s,omega.id,num){
  #obs non-syn change. 
  obs.non = sum(sapply(omega.id, function(x){dat[x[1],x[2]]}))
  f1      = rowSums(dat) + colSums(dat)  
  #exp non-syn change
  gtr.2  = GTR(s, f0)
  mg94.2 = MG94(gtr.2,1,cod,codonstrs,syn,num)
  P94.2  = expm(mg94.2)
  print(sum(P94.2))
  
  exp.non = 0
  for (i in 1:num) {
    ith.non = omega.id[which(sapply(omega.id,'[[',1)==i)]
    exp.non = exp.non + f1[i]*sum(sapply(ith.non, function(x){P94.2[x[1],x[2]]}))
  }
  w = obs.non/(exp.non/2)
  return(w)
}
