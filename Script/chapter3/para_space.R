library(expm)
library(stats)
library(Biostrings)

#GTR matrix
GTR = function(si, pai){
  r1 = matrix(0, 4, 4)
  r1[lower.tri(r1)] = si
  r  = r1 + t(r1)
  r  = t(r*pai)
  diag(r) = -rowSums(r)
  return(r)
}


#######################################
set.seed(8088)
#set up the nuc freq. 
#narrow the freq range so the repeats would be less.
nd  = 100
PiA = runif(nd,0.1,0.4)   
PiT = PiA
PiC = 0.50 - PiA
PiG = PiC
Pi.all = list(PiA, PiC, PiG, PiT)
Pi.all = sapply(1:nd, function(x){sapply(Pi.all, "[[", x)})
#print(colSums(Pi.all)) 


#set up omega from zhou's paper
wv = runif(nd,0.02,0.5)

#keep the brlen between [0,0.1]
#set up 6 sigmas, choose the mean rate as 0.1, lower the cv
cv = 0.5
a  = 1/(cv^2)
b  = 0.1/a                   #>>adjustable     
Sigmas = matrix(0,6,nd)
tv = c()
for (i in 1:nd) {
  pai        = Pi.all[,i]
  si         = rgamma(6, shape=a, scale=b)
  gtr        = GTR(si, pai)
  Sigmas[,i] = si 
  tv[i]      = -sum(diag(gtr)*pai)
}

#summary(tv)
#hist(tv, prob = TRUE, xlim = c(0,2))
#stv = sort(tv)
#lines(stv,dgamma(stv,shape=a,scale=b),col='magenta',lwd=2,lty='dotted')
#plot(density(tv))
#dev.off()

#normalize sigma. 
norm.Sig = sapply(1:nd, function(x){Sigmas[,x]/tv[x]})
#print(norm.Sig)

#set up the indel rate ([12,16]/100)
#assume P_ins = P_del
r0  = runif(nd,0.05,0.15)
r1  = runif(nd,0.05,0.15)
r2  = runif(nd,0.05,0.15)

r       = matrix(c(r0,r1,r2),nd,3)
rt      = r*tv
g.open  = 1-exp(-r*tv)


#setup the gap extension 
#pick a range for avg gap size and convert to extention prob 
#T'=T/3, e=(T'-G)/T' -> T'/G=1/(1-e))
lower =3.5/3
higher=4
avg.gap1 = runif(nd,lower,higher)   
avg.gap2 = runif(nd,lower,higher)  
avg.gap = matrix(c(avg.gap1,avg.gap2),nd)
ext     = 1-1/avg.gap

##output the 100 parameters
tPar = matrix(0,100,17)
for (i in 1:nd) {
  tPar[i,1:4]  = Pi.all[,i]
  tPar[i,5:10] = norm.Sig[,i]
  tPar[i,11]   = wv[i]
  tPar[i,12]   = tv[i]
  tPar[i,13:14]= ext[i,]
  tPar[i,15:17]= rt[i,]
  
}
colnames(tPar)=c('A','C','G','T','s1','s2','s3','s4','s5','s6','omega','tau','ext.I','ext.D','r0*t','r1*t','r2*t')
write.table(tPar,"trueP.100.txt",quote=F,sep="\t",
            row.names=F)

