#Construct codons and its degeneracy
#61 or 64
codon_call = function(num=64){
  stp   = c(49,51,57)
  cod64 = cbind(rep(DNA_BASES, each = 16),
                rep(DNA_BASES, times = 4, each = 4),
                rep(DNA_BASES, 16))
  if(num==61){
    cod = cod64[-stp,]
  }else{
    cod = cod64 
  }
  
  codonstrs  = apply(cod, 1, stringr::str_c, collapse = "")            
  syn        = syncodons(codonstrs)
  names(syn) = toupper(names(syn))
  syn        = lapply(syn, toupper)
  
  
  res = list(cod,codonstrs,syn)
  return(res)
}