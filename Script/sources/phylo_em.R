
#phylo-em algorithm to solve the parameter estimate
phylo_em = function(p){
   source("sources/gtr.R")
   source("sources/mg94.R")
  
    s   = p[1:6]
    w   = p[7]
    r   = GTR(s, f0)
    R   = MG94(r, w, cod, codonstrs, syn, num)
    #print(sum(R) %>% round(13))
    
    #Simulate branch length
    T1  = -sum(diag(r)*f0)          
    #print(T1)
    
    ##make r, R symmetric
    fmat = outer(sqrt(cf0), 1/sqrt(cf0))
    S    = R * fmat
    
    ##calculate eigenvectors and values.
    eig = eigen(S)
    D   = eig$values
    V   = eig$vectors
    
    ##calculate Prob(b|a)
    pab = V %*% diag(exp(D)) %*% t(V)           
    Pab = pab * t(fmat)
    #print(rowSums(Pab))
    
    
    ##Log likelihood "the smaller, the better"
    ll = sum(log(Pab*cf0)*dat)
    cat(sprintf("%.6f\n", ll))
    
    
    
    ##construct the Jkl matrix       
    J = outer(D/T1, D/T1, function(x,y) {
      ifelse(x-y == 0,
             T1*exp(x*T1),
             exp(y*T1)*(expm1((x-y)*T1))/(x-y))
    })
    
    ##calculate the expected values
    # W[a,b,i,i] is the expected time spent state i on a branch going from a -> b
    # U[a,b,i,j] is the expected number of events going from i -> j on a branch going from a->b
    W = array(0, c(num,num,num,num))     
    U = array(0, c(num,num,num,num))
    
    tm = system.time(
      for(a in 1:num) {
        for(b in 1:num) {
          for(i in 1:num) {
            for(j in 1:num) {
              ff = sqrt(cf0[i]*cf0[b]/cf0[a]/cf0[j])
              o  = outer(V[a,]*V[i,], V[j,]*V[b,])   ##cheat: V[i,] = t(V)[,i]
              W[a,b,i,j] = ff * sum(o*J)
            }
          }
          W[a,b,,] = W[a,b,,] / Pab[a,b]
          U[a,b,,] = R * W[a,b,,]
        }
      }
    ) 
    
    ##calculate expected values by summing over observations --a,b is sumable. 
    Wh = array(0, c(num,num))
    Uh = array(0, c(num,num))
    for(i in 1:num) {
      for(j in 1:num) {
        Wh[i,j] = sum(W[,,i,j] * dat)
        Uh[i,j] = sum(U[,,i,j] * dat)
      }
    }
    Wh = diag(Wh)
    
    ##M-Step maximize sigmas.
    sigma.Cij = sapply(seq(6), function(x){sum(Uh[sigma.id[[x]]])})
    sigma.Wij = c()
    for (k in 1:6) {
      ichunks   = ceiling(sigma.id[[k]]/num)
      sigma.Wij = c(sigma.Wij, sum(Wh[ichunks]* t(R)[sigma.id[[k]]])/s[k])
    }
    s = sigma.Cij/sigma.Wij
    
    
    ##M-Step maximize omega
    w.Cij = sum(sapply(omega.id, function(x){Uh[x[1], x[2]]}))
    Rii   = c()
    for (i in 1:num) {
      ith.non = omega.id[which(sapply(omega.id, "[[", 1) == i)]
      ith.sum = sum(sapply(ith.non, function(x){R[x[1], x[2]]})) / w    
      Rii     = c(Rii, ith.sum)
    }
    w.Wij = sum(Wh*Rii)
    w     = w.Cij/w.Wij
    
    ##reconstruct gtr and mg94.
    pNew = c(s,w)
    return(pNew)
}