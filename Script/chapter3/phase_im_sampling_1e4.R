# Take the coatiM (viterbi algo) and weight score as input
# Apply the importance sampling Wi = g/f
# g(x)--ziqi's phase_coati model; f(x)--juan's coatiM
#
# A -- ancestor, B -- descendent
# insertion and deletion length distribution combined.
# inserting bases are ignored since they are treated as deletions.

suppressWarnings(suppressMessages(library(tidyverse)))
suppressPackageStartupMessages(library(Matrix))
suppressPackageStartupMessages(library(stringr))
suppressPackageStartupMessages(library(seqinr))
suppressPackageStartupMessages(library(stringi))
suppressPackageStartupMessages(library(Biostrings))
suppressPackageStartupMessages(library(jsonlite))
suppressPackageStartupMessages(library(dfoptim))
suppressPackageStartupMessages(library(doParallel))
suppressPackageStartupMessages(library(profvis))

#

######################################
#Count the nucleotide freq
count_freq = function(input){
  nuc.count = 0
  for (i in input){
    dna       = readDNAStringSet(i)
    nuc.count = nuc.count + oligonucleotideFrequency(dna,width=1)
  }
  nuc.freq  = colSums(nuc.count)/sum(nuc.count)
  if(sum(nuc.freq)==1){
    return(nuc.freq)
  }else{
    print("The nucleotide frequency sums up not 1!")
  }
}

#Read sample
read_sample = function(D,K){
  Alist = list()
  for (i in 1:K) {
    Alist[[i]] = DNAStringSet(c(D$Seq1[i],D$Seq2[i]))
  }
  return(Alist)
}

#-Log-likelihood of
LL_min  = function(theta){
  nuc_codon_freq = init_f(cod,DD,ncd)
  fw  = nuc_codon_freq[[1]]
  cfw = nuc_codon_freq[[2]]
  
  rmat = GTR(theta[1:6], fw)
  Rmat = MG94(rmat,theta[7],cod,codonstrs,syn,ncd)
  -sum(log(expm(Rmat)*cfw)*DD)
}

#remove effect of scaling factor [ancestor/descendant]
test_pab = function(ab,f0){
  dnaAB    = DNAStringSet(gsub('-','',ab))
  ab.count = oligonucleotideFrequency(dnaAB,width=1)
  sumPab   = sum(ab.count[1,]*log(f0)) + sum(ab.count[2,]*log(f0))
  return(sumPab)
}




##############################################
# ouD   = "JsonD/58/"
# inD   = "Gs_trim/58"
# ouD   = "hmr/JsonD/3/"
# inD   = "hmr/Data/3"

# inD   ="../test_90_species/Raw_data/cds/06_Nematode_aligned_cds"
# ouD   ="90/JsonD/06_Nematode_aligned_cds/"
# ouF   ="90/Results/PISE/06_Nematode_aligned_cds.est.json"
# ssize ="100"
# ncdon ='61'
# spList='90/species/06_Nematode_aligned_cds.txt'
main = function(inD,ouD,spList,ouF,ssize,ncdon){
  
  sub1 = "../Script/sources/"
  sub2 = "../Script/chapter3/"
  source(paste0(sub1,"codon_call.R"))
  source(paste0(sub1,"gtr.R"))
  source(paste0(sub1,"mg94.R"))
  source(paste0(sub1,"init_f.R"))
  source(paste0(sub1,"LL.R"))
  source(paste0(sub2,"phase_indel_prob3.R"))
  
  #construct codons and its degeneracy
  ncd      <<- as.numeric(ncdon)
  co.res    = codon_call(ncd)
  codons    = co.res[[1]]
  codonstrs = co.res[[2]]
  syn       = co.res[[3]]
  
  #nmkb method
  cod      <<- codons
  codonstrs<<- codonstrs
  syn      <<- syn
  
  ######################################
  Files = list.files(inD, full.names=T)
  namev = str_extract(basename(Files),'[^.]+')
  
  if(namev[1]=='1'){#simulation
    index = as.numeric(namev)
    Files = Files[order(index)]
    tag   = as.numeric(basename(inD))
  }else{#90 species
    Fbname = list.files(inD, full.names=F)
    kept   = read.table(spList,header=F)[,1]
    index  = which(Fbname %in% kept)
    Files  = Files[index]
  }
  n     = length(Files)

   
  #Initial parameters
  set.seed(8088)
  f0  = count_freq(Files)
  p0  = rep(0.1,7)
  
  g0  = rep(0.01,6)
  e0  = rep(0.8,2)      #single codon
  
  rmat= GTR(p0[1:6],f0)
  t0  = -sum(diag(rmat)*f0) 
  
  #Iterate through the black box
  ##default:ssize=100
  max.it = 20
  K      = as.numeric(ssize)
  lljv   = matrix(NA,n,max.it)
  llzv   = matrix(NA,n,max.it)
  
  N      = matrix(NA,n,6)
  M      = matrix(NA,n,6)
  E      = matrix(NA,n,2)
  colnames(N) = c('i0','i1','i2','d0','d1','d2')
  colnames(M) = c('no_i0','no_i1','no_i2','no_d0','no_d1','no_d2')
  colnames(E) = c('avg.gap.len.I','avg.gap.len.D')
  
  Wv     = list()          
  
  
  #>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>RUNNING
  iter=1
  tm = system.time({
    repeat{
      D = matrix(0,ncd,ncd)
      if(iter>max.it){
        cat("Pass the", max.it, "iterations limits!")
        break}
      
      #pre-cal avg gap size
      avg.gap = 1/(1-e0)
      if(any(avg.gap<3)){
        print("Average gap length less than 3!")
        break
        }
      e3    = 1-1/(avg.gap/3)
      
      #prellocate gtr/mg94 mat, codon freq.
      cf0   = sapply(seq(ncd), function(x){prod(f0[match(cod[x,],DNA_BASES)])})   
      cf0   = cf0/sum(cf0)
      Rmat  = MG94(rmat,p0[7],cod,codonstrs,syn,ncd)
      Pmat  = log(expm(Rmat))
      
      for (j in 1:n) {#E-step
        ##coati-sampler
        input = Files[j]
        ouJ   = paste0(ouD,j,'.json')
        cmd   = paste("bash ../Script/chapter3/coati_sampler.sh", input, ouJ, K,
                      mean(g0),mean(e0),f0[1],f0[2],f0[3],f0[4],p0[1],p0[2],p0[3],p0[4],p0[5],p0[6],p0[7],t0,sep=' ')
        system(cmd)
        Input.tmp    = fromJSON(ouJ)
        Dat.tmp      = Input.tmp %>% tidyr::unpack(aln)
        
        colnames(Dat.tmp)[1:2] = c('Seq1','Seq2')
        aln      = Dat.tmp %>% dplyr::group_by(Seq1,Seq2)
        aln2     = aln %>% dplyr::group_keys(Seq1,Seq2,weight,log_weight)
        gpsize   = group_size(aln)
        ngroups  = n_groups(aln)
        
        print(ngroups)
      
        llj          = aln2$log_weight
        lljv[j,iter] = sum(llj*gpsize)/K
        Alist        = read_sample(aln2,ngroups)
        
        #llj2 = exp(llj)/sum(exp(llj)*gpsize)
        #print(sum(llj2*gpsize))
        
        ##preset summary stats
        Gm        = matrix(0,ngroups,6)
        Mm        = matrix(0,ngroups,6)
        codon_arr = array(0,c(ncd,ncd,ngroups))
        llz       = rep(0,ngroups)
        gl        = matrix(0,ngroups,2)
        
        ##adding scalling factor 
        scal.pab = test_pab(Alist[[1]],f0)
        
        
        #profvis({
          for(i in 1:ngroups){
            A      = Alist[[i]]
            res1   = ziqi_prob(A,g0,e3,Pmat,codonstrs,f0,ncd)
            Gm[i,] = res1[[1]]
            Mm[i,] = res1[[2]]
            codon_arr[,,i] = res1[[3]]
            llz[i] = res1[[4]]-scal.pab
            gl[i,] = res1[[5]]
          } 
        #})
         
        
        llzv[j,iter] = sum(llz*gpsize)/K
        if(any(llz==-Inf)){#low-quality alignment
          next
        }
        
        ############
        #cal. weight
        uwi          = llz-llj
        uwi.max      = max(llz-llj)
        uwi.dif      = uwi - uwi.max
        Wi           = exp(uwi.dif)/sum(exp(uwi.dif)*gpsize)
        print(sum(Wi*gpsize))
        Wv[[j]]      = Wi*gpsize
        
      
        ##weighted gap phases
        N[j,] = colSums(Wi*gpsize*Gm)
        M[j,] = colSums(Wi*gpsize*Mm)

        ##weighted average gap length
        #E[j,] = colSums(Wi*gl*gpsize,na.rm=T) #>
        E[j,] = colSums(Wi*gl*gpsize) 
        
        datw = matrix(0,ncd,ncd)    
        for (j in 1:ngroups) {
          dat.wei = Wi[j]*codon_arr[,,j]
          datw    = datw + dat.wei*gpsize[j]
        }
        D = D+datw/n
      }
      
      
      
      
      #######################
      #M step: para estimates
      ##gap opening 
      nnew = colMeans(N,na.rm=T)
      mnew = colMeans(M,na.rm=T)
      gnew = nnew/(nnew+mnew)
      
      
      ##gap extension
      #w.avg.gap = c(mean(E[which(E[,1]!=0),1]),mean(E[which(E[,2]!=0),2]))
      w.avg.gap = c(mean(E[which(E[,1]>=1),1]),mean(E[which(E[,2]>=1),2]))
      enew1     = 1 - 1/(w.avg.gap*3)
      enew      = 1 - 1/w.avg.gap
      if(any(enew<0)){
        print("Warning: Updated gap extension prob (unit of 3) < 0!")
        break
        }
      
      DD <<- D
      pb = nmkb(fn=LL_min, par=p0, lower=0, control=list(tol=1e-5,trace=F,maxfeval=5000))   #change the tolerance
      if(pb$convergence != 0){
        cat("Warning: failed convergence!")
      }else{
        pnew = pb$par
      }
      
      ##cal. the tau
      nuc_f = init_f(cod,DD,ncd)     
      fnew  = nuc_f[[1]]
      rmat  = GTR(pnew[1:6],fnew)
      tnew  = -sum(diag(rmat)*fnew)
      
      #print the results
      print(gnew[1]+gnew[4])
      print(gnew[2]+gnew[5])
      print(gnew[3]+gnew[6])
      print(mean(enew))
      print(pnew[1:6]/tnew)
      print(pnew[7])
      print(tnew)
      
      #record the pars
      # gv   = rbind(gv,gnew)
      # pv   = rbind(pv,pnew)
      # ev   = c(ev,enew)
      # tv   = c(tv,tnew)
      
      #rmse tolerance
      p    = c(g0,e0,p0)
      q    = c(gnew,enew1,pnew)
      delta= (q-p)^2
      rmse = sqrt(mean(delta))
      cat(sprintf("iter:%i, rmse:%.6f\n",iter, rmse))
      if(rmse<=1e-4){
        break
      }else{
        # if(iter>10){
        #   K=1e+4
        #   if(tag==98){
        #     max.it=11
        #   }
        # }
        iter= iter+1
        f0  = fnew
        p0  = pnew
        g0  = gnew
        e0  = enew1
        t0  = tnew
      }
    }
  })
  
  cat(sprintf("Running loops: %i\n  Running time:%.3f mins", iter, tm[3]/60))
  ################################################################################
  
  #output the parameter estimates summary stats
  par_lst = list('nuc.freq'=fnew, 'sigmas'=pnew[1:6]/tnew, 'omega'=pnew[7], 'branch.length'=tnew, 
                 'gap.openning'=gnew,'gap.extension'=enew)
  sum_lst = list('gap.phases'=nnew,'avg.gap.size'=w.avg.gap)
  
  par_est = toJSON(par_lst)
  sum_stat= toJSON(sum_lst)
  
  dname = dirname(ouF)
  fname = str_extract(basename(ouF),'[^.]+')
  ouF2  = paste0(dname,'/',fname,'.sum.json')
  write(par_est,ouF)
  write(sum_stat,ouF2)
  
  #output the weight matrix (last weight)
  dname2 = gsub('(.*)/\\w+', '\\1', dname)
  ouF3   = paste0(dname2,'/weightM/',fname,'.txt')
  lapply(Wv,write,ouF3,append=T,ncolumns=1e+4)
  
  #output rmse
  ouF3   = paste0(dname2,'/RMSE/',fname,'.txt')
  write(rmse,ouF3,ncolumns=1)
  
  #output llzv matrix
  ouF4   = paste0(dname2,'/LLZ/',fname,'.tsv')
  write.table(llzv,ouF4,col.names=F,row.names=F,sep='\t')
  
}


################################
args = commandArgs(trailingOnly=T)
main(args[1],args[2],args[3],args[4],args[5],args[6])





# source(paste0(sub1,"count_gtr_mat.R"))
# 
# 
# aFiles = list.files( "../test_90_species/Raw_data/align_max/06_Nematode_aligned_cds",full.names = T)
# Cnt = matrix(0,4,4)
# for(i in 1:length(aFiles)){
#   A   = readDNAStringSet(aFiles[i])
#   Cnt = Cnt + count.nsub(A)
#   print(i)
# }

