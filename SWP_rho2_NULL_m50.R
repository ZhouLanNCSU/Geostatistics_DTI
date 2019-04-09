library("foreach")
library("doParallel")
registerDoParallel(cores=16)

SWP_rho2_NULL_m50<-foreach(st=1:50,.errorhandling=c('pass') )%dopar%{
set.seed(st)
  source("/mnt/home/zlan/SWP_all/New_Simulation/CODEV2.R")
  library("bayess",lib.loc="/mnt/home/zlan/R-Library")
 correlation=2
  NR=400
  size=20
  seqx=seq(1, size, length.out = size)
  seqy=seq(1, size, length.out = size)
  s <- expand.grid(x = seqx, y = seqy)
  
  c=size
  g_full <- graph.lattice(length=c,dim=2)
  n=c^2
  #V(g_full)$index=1:n
  net=get.adjacency(g_full,attr=NULL)
  K=5
  rho=5
  label_T0<-pottshm(ncol=K,niter=1,c,m=c,beta=rho)
  label_T1<-pottshm(ncol=K,niter=1,c,m=c,beta=rho)
  cc=c/3
  lower=0.25
  label_T0=label_T0*0
  label_T0[10:15,10:15]=1
  
 
  cov=varcov.spatial(coords =s,cov.model = "exponential",cov.pars=c(1,correlation))$`varcov`
  
  nonzero=which(label_T0==1)
  p_ind=rep(0,c^2)
  p_ind[nonzero]=0.5
  p_age=rep(0.25,c^2)
  
  beta11=t(cbind(rmvn(n = 1, mu=rep(0,size^2),  cov*0.01),
                 rmvn(n = 1, mu=p_ind,  cov*0.01),
                 rmvn(n = 1, mu=p_age,  cov*0.01)))
  
  beta22=t(cbind(rmvn(n = 1, mu=rep(0,size^2),  cov*0.01),
                 rmvn(n = 1, mu=p_ind,  cov*0.01),
                 rmvn(n = 1, mu=p_age,  cov*0.01)))
  
  beta33=t(cbind(rmvn(n = 1, mu=rep(0,size^2),  cov*0.01),
                 rmvn(n = 1, mu=p_ind,  cov*0.01),
                 rmvn(n = 1, mu=p_age,  cov*0.01)))
  
  beta21=t(cbind(rmvn(n = 1, mu=rep(0,size^2),  cov*0.01),
                 rmvn(n = 1, mu=rep(0,size^2),  cov*0.01),
                 rmvn(n = 1, mu=p_age,  cov*0.01)))
  
  beta31=t(cbind(rmvn(n = 1, mu=rep(0,size^2),  cov*0.01),
                 rmvn(n = 1, mu=rep(0,size^2),  cov*0.01),
                 rmvn(n = 1, mu=p_age,  cov*0.01)))
  
  beta32=t(cbind(rmvn(n = 1, mu=rep(0,size^2),  cov*0.01),
                 rmvn(n = 1, mu=rep(0,size^2),  cov*0.01),
                 rmvn(n = 1, mu=p_age,  cov*0.01)))
  N=10
  X=lapply(1:N,function(sub) c(1, ifelse(sub>5, 0, 1), abs(rnorm(1,0,1)) )  )
  
  
  Data=lapply(1:N, function(i) NULL)
  for (sub in 1:N){
    
    D1=rmvn(n = 50, mu=rep(0,c^2), cov)
    D2=rmvn(n = 50, mu=rep(0,c^2),  cov)
    D3=rmvn(n = 50, mu=rep(0,c^2),  cov)
    
    for(v in 1:c^2){
      
      L=matrix(0,3,3)
      L[1,1]=exp(X[[sub]]%*%beta11[,v])
      L[2,2]=exp(X[[sub]]%*%beta22[,v])
      L[3,3]=exp(X[[sub]]%*%beta33[,v])
      L[2,1]=X[[sub]]%*%beta21[,v]
      L[3,1]=X[[sub]]%*%beta31[,v]
      L[3,2]=X[[sub]]%*%beta32[,v]
      
      mm=cbind(D1[,v],D2[,v],D3[,v])
      Data[[sub]][[v]]<-L%*%t(mm)%*%mm%*%t(L)/50
    }
    
  }
  
  rho=mean_rho=correlation=2
  result=Working_Model(Data,s,X,
                        mean_nu=-1,sd_nu=1,
                        mean_range=0,sd_range=1,
                        a_var=.01,b_var=.01,
                        iters=5000,burn=1000,
                        NR=NR) 
  result
}


save.image("SWP_rho2_NULL_m50.Rdata")