#Build sigma list for M-step
##Locating all 6-sigma locations in 64*64/61*61 R matrix.
sigma_coor = function(cod, num){
  ii = GTR(1:6, rep(1,4))
  diag(ii) = 0
  I = matrix(0, num, num)
  for (i in 1:num) {
    for (j in 1:num) {
      if(i == j){
        I[i, j] = 0
      }else if(sum(cod[i, ] != cod[j, ]) > 1){
        I[i, j] = 0
      }else{
        pos = which(cod[i, ] != cod[j, ])
        x   = which(DNA_BASES == cod[i, pos])
        y   = which(DNA_BASES == cod[j, pos])
        I[i, j] = ii[x, y]
      }
    }
  }
  sigma.id = sapply(1:6, function(x){which(I == x)})
  if(num==64){
    sigma.id = split(sigma.id,col(sigma.id))
  }
  return(sigma.id)
}

  
