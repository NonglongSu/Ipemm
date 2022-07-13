#Build omega list for M-step.
##Locating all non-syn locations in 64*64 R matrix.

omega_coor = function(cod, codonstrs, syn, num){
  omega.id = c()
  for (i in 1:num) {
    for (j in 1:num) {
      if((i != j) && 
         (sum(cod[i, ] != cod[j, ]) == 1) && 
         (!(codonstrs[j] %in% syn[[codonstrs[i]]])) ){
        omega.id = c(omega.id, i, j)          
      }
    }
  }
  omega.id = split(omega.id, ceiling(seq_along(omega.id) / 2))
  return(omega.id)
}


