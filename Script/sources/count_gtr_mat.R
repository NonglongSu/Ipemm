#Create a 4*4 nuc subs matrix
#input:a pairwise alignment

#remove the gap positions
rm.gap = function(A,As,g){
  g    = unlist(g)
  if(length(g)==0){
    return(A)
  }else{
    rm.A = stri_sub_replace_all(A,from=sort(start(g)),to=sort(end(g)),replacement='')
    rm.A = DNAStringSet(rm.A)
    if(all(width(rm.A)%%3 == 0)){
      return(rm.A)
    }else{
      print("Warning:removed-gap sequences are not multiple of three")
      break
    }
  }
}

#generate a 4*4 gtr count
count.nsub = function(A){
  As = str_split(A,'')
  g  = IRangesList(lapply(As, function(x){IRanges(x=='-')}))
  rA = rmgap(A,As,g)
  
  nmat = matrix(0,4,4)
  colnames(nmat)=DNA_BASES
  rownames(nmat)=DNA_BASES
  
  seqs = str_split(rA,'')
  s1   = seqs[[1]]
  s2   = seqs[[2]]
  len  = length(seqs[[1]])
  for(i in 1:len){
    nmat[s1[i],s2[i]] = nmat[s1[i],s2[i]]+1
  }
  return(nmat)
}




