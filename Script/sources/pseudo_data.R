#Generate pseudo data
pseudo_data = function(ssize,pmat,num){
  dat  = sample(num*num, ssize, replace=TRUE, prob=pmat)
  dat  = table(dat)
  dat  = as.data.frame(dat)
  id1  = as.numeric(as.vector(dat[[1]]))
  id2  = as.numeric(as.vector(dat[[2]]))
  dat1 = matrix(0, num, num)
  for (i in 1:length(id1)) {
    dat1[id1[i]] = id2[i]
  }
  return(dat1)
}

