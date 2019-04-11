library("rlang")
library('igraph')
library('rlist')
library('MASS')
library('matrixcalc')
library('MCMCpack')
library('inline')
library('RcppArmadillo')
library('Rcpp')
library("rstan")
library("Boom")
#library("profvis")
library("mclust")
library("plot3D")
library("geoR")
library("Matrix")
library("foreach")
library("doParallel")
library("ramcmc")


CPP_inv.body<-"
RNGScope scope;
Rcpp::List mylist(list_input);
int pp=as<int>(parts_input);

for(int i=0; i<pp; i++) {
SEXP ll = mylist[i];
arma :: mat ll_c=Rcpp :: as < arma :: mat >(ll);
mylist[i]=inv(ll_c);
}

return(wrap(mylist));
"
CPP_inv<-cxxfunction(signature(list_input="List",parts_input="integer"), body = CPP_inv.body,plugin="RcppArmadillo")

CPP_inv<-function(cov_list,parts){
  lapply(1:parts,function(i) solve(cov_list[[i]]))
}


residual_inv<-function(s_list,rho,smooth,parts){
  
  cov_list=lapply(1:parts, function(pp) varcov.spatial(coords=s_list[[pp]],cov.model = "matern",kappa=smooth,cov.pars=c(1,rho))$`varcov`^2)
  III=CPP_inv(cov_list,parts)
  INV=  duplicate(III,shallow = FALSE)
  #INV=lapply(cov_list, solve)
  return(INV)
}

mean_inv<-function(s_list,rho,smooth,parts){
  
  cov_list=lapply(1:parts, function(pp) varcov.spatial(coords=s_list[[pp]],cov.model = "matern",kappa=smooth,cov.pars=c(1,rho))$`varcov`)
  III=CPP_inv(cov_list,parts)
  INV=  duplicate(III,shallow = FALSE)
  #INV=lapply(cov_list, solve)
  return(INV)
}


logdet=function(X){
  (log(det(X)))
}



mysum<-function(X){
  
  sum(X,na.rm = TRUE)
}


Working_Model<-function(Data,s,X,
                        mean_nu=-1,sd_nu=1,
                        mean_range=0,sd_range=1,
                        a_var=.01,b_var=.01,
                        iters=5000,burn=1000,
                        NR=50){
  
  ##########################################
  ####### OTHERS ###########################
  ##########################################
  cat("OTHERS", "\n")
  N        <- length(Data)
  p        <- length(X[[1]])
  n=nrow(s)
  parts=(n/NR)
  partition_list=lapply(1:parts,function(pp) ((pp-1)*NR+1):(NR*pp) )
  s_list=lapply(1:parts, function(pp) s[partition_list[[pp]],] )
  ddd=as(diag(NR),"sparseMatrix")
  DDD=as(diag(n),"sparseMatrix")
  blockX_list=lapply(1:N, function(sub) ddd%x%matrix(X[[sub]],nrow=1) )
  FullX_list=lapply(1:N, function(sub) DDD%x%matrix(X[[sub]],nrow=1) )
  
  ##########################################
  ####### Data Preparation#################
  ##########################################
  cat("Data Preparation", "\n")
  
  Chol=lapply(1:N,function(sub) lapply(Data[[sub]],chol))
  t11=sapply(1:N, function(sub) sapply(1:n, function(v) Chol[[sub]][[v]][1,1]))
  t22=sapply(1:N, function(sub) sapply(1:n, function(v) Chol[[sub]][[v]][2,2]))
  t33=sapply(1:N, function(sub) sapply(1:n, function(v) Chol[[sub]][[v]][3,3]))
  t21=sapply(1:N, function(sub) sapply(1:n, function(v) Chol[[sub]][[v]][1,2]))
  t31=sapply(1:N, function(sub) sapply(1:n, function(v) Chol[[sub]][[v]][1,3]))
  t32=sapply(1:N, function(sub) sapply(1:n, function(v) Chol[[sub]][[v]][2,3]))
  
  logt11=log(t11)
  logt22=log(t22)
  logt33=log(t33)
  logtii_list=list(logt11,logt22,logt33)
  
  q2logt11=sqrt(2)*logt11
  q2logt22=sqrt(2)*logt22
  q2logt33=sqrt(2)*logt33
  
  q2logtii_list=list(q2logt11,q2logt22,q2logt33)
  
  ##########################################
  ####### Initial Values####################
  ##########################################
  cat("Initial Values", "\n")
  ##Coefficents a p*n vector
  plugin_beta11=beta11= sapply(1:n,function(v) coef(lm(log((t11[v,]))~0+list.rbind(X))))
  plugin_beta22=beta22= sapply(1:n,function(v) coef(lm(log((t22[v,]))~0+list.rbind(X))))
  plugin_beta33=beta33= sapply(1:n,function(v) coef(lm(log((t33[v,]))~0+list.rbind(X))))
  plugin_beta21=beta21= sapply(1:n,function(v) coef(lm(t21[v,]~0+list.rbind(X))))
  plugin_beta31=beta31= sapply(1:n,function(v) coef(lm(t31[v,]~0+list.rbind(X))))
  plugin_beta32=beta32= sapply(1:n,function(v) coef(lm(t32[v,]~0+list.rbind(X))))
  
  plugin_list=list(plugin_beta11,plugin_beta22,plugin_beta33,plugin_beta21,plugin_beta31,plugin_beta32)
  
  
  sd_beta11= sapply(1:n,function(v) summary(lm(log((t11[v,]))~0+list.rbind(X)))$sigma)
  sd_beta22= sapply(1:n,function(v) summary(lm(log((t22[v,]))~0+list.rbind(X)))$sigma)
  sd_beta33= sapply(1:n,function(v) summary(lm(log((t33[v,]))~0+list.rbind(X)))$sigma)
  sd_beta21= sapply(1:n,function(v) summary(lm(t21[v,]~0+list.rbind(X)))$sigma)
  sd_beta31= sapply(1:n,function(v) summary(lm(t31[v,]~0+list.rbind(X)))$sigma)
  sd_beta32= sapply(1:n,function(v) summary(lm(t32[v,]~0+list.rbind(X)))$sigma)
  
  
  
  ##residual 
  var_m=var(c(q2logt11,q2logt22,q2logt33,t21,t32,t31))
  rho=2
  smooth=0.5
  INV=residual_inv(s_list,rho,smooth,parts)
  INV_LOGDET=lapply(INV, logdet)
  
  ##mean
  s_beta=s_beta0=mean(c(sd_beta11,sd_beta22,sd_beta33,sd_beta31,sd_beta32,sd_beta21))^2
  mean_rho=2
  mean_smooth=0.5
  mean_INV=residual_inv(s_list,mean_rho,mean_smooth,parts)
  mean_INV_LOGDET=lapply(mean_INV, logdet)
  
  
  
  ###other important things
  residual22=lapply(1:N, function(sub)  sapply(1:n, function(v) exp(-X[[sub]]%*%plugin_beta22[,v]) ) )
  residual33=lapply(1:N, function(sub)  sapply(1:n, function(v) exp(-X[[sub]]%*%plugin_beta33[,v]) ) )
  
  d11=lapply(1:N, function(sub) sapply(1:n, function(v) t11[v,sub]*exp(-X[[sub]]%*%plugin_beta11[,v]) ) )
  d22=lapply(1:N, function(sub) sapply(1:n, function(v) t22[v,sub]*exp(-X[[sub]]%*%plugin_beta22[,v]) ) )
  
  ##########################################
  ####### MCMC##############################
  ##########################################
  MCMC_rho=rep(0,iters)
  MCMC_smooth=rep(0,iters)
  MCMC_var_m=rep(0,iters)
  
  MCMC_mean_rho=rep(0,iters)
  MCMC_mean_smooth=rep(0,iters)
  MCMC_s_beta=rep(0,iters)
  
  MCMC_beta=list("beta11"=lapply(1:iters,function(it) NULL),
                 "beta22"=lapply(1:iters,function(it) NULL),
                 "beta33"=lapply(1:iters,function(it) NULL),
                 "beta21"=lapply(1:iters,function(it) NULL),
                 "beta31"=lapply(1:iters,function(it) NULL),
                 "beta32"=lapply(1:iters,function(it) NULL))
  
  
  
  flush.console()       
  pb <- txtProgressBar(min = 0, max = iters, style = 3)
  
  AccepRateSpatial=0
  AccepRateMean=0
  S=diag(2)
  
  for(it in 1:iters){
    
    setTxtProgressBar(pb, it)
    flush.console()
    cat("Residual Spatial Acceptance Rate", AccepRateSpatial/it, ", the range is", rho,", the smoothness is" , smooth,"\n")
    cat("Mean Spatial Acceptance Rate", AccepRateMean/it, ", the range is", mean_rho,", the smoothness is" , mean_smooth,"\n")
    cat("Residual Variance is", var_m,"\n")
    cat("Mean Variance is", s_beta,"\n")
    
    ##########################################
    #######Mean Spatial Prameters ############
    ##########################################
    if(it==1|sample(c(0,1),1,prob=c(0,1))==1){
    ####variance
    a<- n*N/2+a_var
    b<- b_var
    beta_all_list=list(beta11,beta22,beta33,beta21,beta31,beta32)
    DDD=bdiag(mean_INV)%x% Diagonal(p)
    for(ii in 1:6){
      for(sub in 1:N){
        R      <- as.numeric(beta_all_list[[ii]])-as.numeric(plugin_list[[ii]])
        b      <- t(R)%*% (DDD)%*%R/2+b
      }
      
    }
    s_beta=1/rgamma(1,a,as.numeric(b))
    }
    s_beta=rgamma(1,10,10)*s_beta0
    ####covariance
    u <- rnorm(2)*0.05
    theta_prop <- log(c(mean_rho,mean_smooth)) + S %*% u 
    can_mean_rho= exp(theta_prop[1])
    can_mean_smooth=exp(theta_prop[2])
    
    can_mean_INV=residual_inv(s_list,can_mean_rho,can_mean_smooth,parts)
    can_mean_INV_LOGDET=lapply(can_mean_INV, logdet)
    
    can=Reduce("+",can_mean_INV_LOGDET)*0.5*6*p+
      dnorm(log(can_mean_rho),mean_range,sd_range,TRUE)+
      dnorm(log(can_mean_smooth),mean_nu,sd_nu,TRUE)
    old=Reduce("+",mean_INV_LOGDET)*0.5*6*p+
      dnorm(log(mean_rho),mean_range,sd_range,TRUE)+
      dnorm(log(mean_smooth),mean_nu,sd_nu,TRUE)
    
    can_DDD=bdiag(can_mean_INV)%x% Diagonal(p)
    for(ii in 1:6){
      for(sub in 1:N){
        R <- as.numeric(beta_all_list[[ii]])-as.numeric(plugin_list[[ii]])
        can=can-0.5*t(R)%*%(can_DDD)%*%R/s_beta
        old=old-0.5*t(R)%*%(DDD)%*%R/s_beta
      }
      
    }
    can=as.numeric(can)
    old=as.numeric(old)
    prob=min(1,exp(can-old))
    if(is.infinite(prob) & can>old){
      prob=0
      Accept=1
    }else{
      Accept=sample(c(1,0),1,prob = c(prob,1-prob))
    }
    
    if(Accept==1){
      mean_rho=can_mean_rho
      mean_smooth=can_mean_smooth
      mean_INV=can_mean_INV
      mean_INV_LOGDET=can_mean_INV_LOGDET
      AccepRateMean=AccepRateMean+1
      
    }
    
    
    
    
    ##########################################
    #######Residual Spatial Prameters ########
    ##########################################
    u <- rnorm(2)*0.05
    theta_prop <- log(c(rho,smooth)) + S %*% u 
    can_rho= exp(theta_prop[1])
    can_smooth=exp(theta_prop[2])
    
    can_INV=residual_inv(s_list,can_rho,can_smooth,parts)
    can_INV_LOGDET=lapply(can_INV, logdet)
    
    can=Reduce("+",can_INV_LOGDET)*0.5*N*6+
      dnorm(log(can_rho),mean_range,sd_range,TRUE)+
      dnorm(log(can_smooth),mean_nu,sd_nu,TRUE)
    old=Reduce("+",INV_LOGDET)*0.5*N*6+
      dnorm(log(rho),mean_range,sd_range,TRUE)+
      dnorm(log(smooth),mean_nu,sd_nu,TRUE)
    
    ####### diagonals ################
    beta_ii_list=list(beta11,beta22,beta33)
    for (ii in 1:3){
      logtii=logtii_list[[ii]]
      betaii=beta_ii_list[[ii]]
      mean=logtii-sapply(1:N,function(sub) sapply(1:n,function(v) X[[sub]]%*%betaii[,v]))
      can=can-
        0.5*sum(sapply(1:N,function(sub) mysum(sapply(1:parts,function(pp) (mean[partition_list[[pp]],sub]%*%can_INV[[pp]])%*%mean[partition_list[[pp]],sub]/(var_m/2)))))
      old=old-
        0.5*sum(sapply(1:N,function(sub) mysum(sapply(1:parts,function(pp) (mean[partition_list[[pp]],sub]%*%INV[[pp]])%*%mean[partition_list[[pp]],sub]/(var_m/2)))))
    }
    
    ####### off-diagonals ################
    ###t21
    central_t21=t21-t11*sapply(1:N, function(sub) sapply(1:n, function(v) exp(-X[[sub]]%*%plugin_beta11[,v])) )*
      sapply(1:N, function(sub) sapply(1:n, function(v) X[[sub]]%*%beta21[,v]) )
    can=can-
      0.5*mysum(sapply(1:N,function(sub) mysum(sapply(1:parts,function(pp) (central_t21[partition_list[[pp]],sub]%*%(diag(sapply(partition_list[[pp]], function(v) exp(-X[[sub]]%*%plugin_beta22[,v])))
                                                                                                                     %*%can_INV[[pp]]%*%
                                                                                                                       diag(sapply(partition_list[[pp]], function(v) exp(-X[[sub]]%*%plugin_beta22[,v])))))%*%
                                                        central_t21[partition_list[[pp]],sub]/(var_m))
      )
      )
      )
    old=old-
      0.5*mysum(sapply(1:N,function(sub) mysum(sapply(1:parts,function(pp) (central_t21[partition_list[[pp]],sub]%*%(diag(sapply(partition_list[[pp]], function(v) exp(-X[[sub]]%*%plugin_beta22[,v])))
                                                                                                                     %*%INV[[pp]]%*%
                                                                                                                       diag(sapply(partition_list[[pp]], function(v) exp(-X[[sub]]%*%plugin_beta22[,v])))))%*%
                                                        central_t21[partition_list[[pp]],sub]/(var_m))
      )
      )
      )
    ### t31
    central_t31=t31-t11*sapply(1:N, function(sub) sapply(1:n, function(v) exp(-X[[sub]]%*%plugin_beta11[,v])) )*
      sapply(1:N, function(sub) sapply(1:n, function(v) X[[sub]]%*%beta31[,v]) )
    can=can-
      0.5*mysum(sapply(1:N,function(sub) mysum(sapply(1:parts,function(pp) (central_t31[partition_list[[pp]],sub]%*%(diag(sapply(partition_list[[pp]], function(v) exp(-X[[sub]]%*%plugin_beta33[,v])))
                                                                                                                     %*%can_INV[[pp]]%*%
                                                                                                                       diag(sapply(partition_list[[pp]], function(v) exp(-X[[sub]]%*%plugin_beta33[,v])))))%*%
                                                        central_t31[partition_list[[pp]],sub]/(var_m))
      )
      )
      )
    old=old-
      0.5*mysum(sapply(1:N,function(sub) mysum(sapply(1:parts,function(pp) (central_t31[partition_list[[pp]],sub]%*%(diag(sapply(partition_list[[pp]], function(v) exp(-X[[sub]]%*%plugin_beta33[,v])))
                                                                                                                     %*%INV[[pp]]%*%
                                                                                                                       diag(sapply(partition_list[[pp]], function(v) exp(-X[[sub]]%*%plugin_beta33[,v])))))%*%
                                                        central_t31[partition_list[[pp]],sub]/(var_m))
      )
      )
      )
    
    ### t32
    central_t32=t32-t22*sapply(1:N, function(sub) sapply(1:n, function(v) exp(-X[[sub]]%*%plugin_beta22[,v])) )*
      sapply(1:N, function(sub) sapply(1:n, function(v) X[[sub]]%*%beta32[,v]) )
    can=can-
      0.5*mysum(sapply(1:N,function(sub) mysum(sapply(1:parts,function(pp) (central_t32[partition_list[[pp]],sub]%*%(diag(sapply(partition_list[[pp]], function(v) exp(-X[[sub]]%*%plugin_beta33[,v])))
                                                                                                                     %*%can_INV[[pp]]%*%
                                                                                                                       diag(sapply(partition_list[[pp]], function(v) exp(-X[[sub]]%*%plugin_beta33[,v])))))%*%
                                                        central_t32[partition_list[[pp]],sub]/(var_m))
      )
      )
      )
    old=old-
      0.5*mysum(sapply(1:N,function(sub) mysum(sapply(1:parts,function(pp) (central_t32[partition_list[[pp]],sub]%*%(diag(sapply(partition_list[[pp]], function(v) exp(-X[[sub]]%*%plugin_beta33[,v])))
                                                                                                                     %*%INV[[pp]]%*%
                                                                                                                       diag(sapply(partition_list[[pp]], function(v) exp(-X[[sub]]%*%plugin_beta33[,v])))))%*%
                                                        central_t32[partition_list[[pp]],sub]/(var_m))
      )
      )
      )
    prob=min(1,exp(can-old))
    if(is.infinite(prob) & can>old){
      prob=0
      Accept=1
    }else{
      Accept=sample(c(1,0),1,prob = c(prob,1-prob))
    }
    
    if(Accept==1){
      rho=can_rho
      smooth=can_smooth
      INV=can_INV
      INV_LOGDET=can_INV_LOGDET
      AccepRateSpatial=AccepRateSpatial+1
      
    }
    
    #S <- ramcmc::adapt_S(S, u, prob, it - 1)
    
    
    
    
    ##############################################:
    #####          VARIANCE (Gibbs)        #######:
    ##############################################:
    
    #var_m=1/50*rgamma(1,10,10)
    if(it==1|sample(c(0,1),1,prob=c(0,1))==1){
    a<- 6*N*n/2+a_var
    b<- b_var
    ###diagonals
    beta_ii_list=list(beta11,beta22,beta33)
    for(ii in 1:3){
      y=q2logtii_list[[ii]]
      for(sub in 1:N){
        R      <- y[,sub]-sqrt(2)*FullX_list[[sub]]%*%as.numeric(beta_ii_list[[ii]])
        b      <- t(R)%*%bdiag(INV)%*%R/2+b
      }
      
    }
    ###off-diagonals
    #t21
    R_full=central_t21
    for(sub in 1:N){
      R=R_full[,sub]
      PPP=Diagonal(n)*(residual22[[sub]])
      b<- t(R)%*% ( PPP%*% bdiag(INV) %*% PPP )  %*%R/2+b
    }
    
    #t31
    R_full=central_t31
    for(sub in 1:N){
      R=R_full[,sub]
      PPP=Diagonal(n)*(residual33[[sub]])
      b<- t(R)%*% ( PPP%*% bdiag(INV) %*% PPP )  %*%R/2+b
    }
    
    #t32
    R_full=central_t32
    for(sub in 1:N){
      R=R_full[,sub]
      PPP=Diagonal(n)*(residual33[[sub]])
      b<- t(R)%*% ( PPP%*% bdiag(INV) %*% PPP )  %*%R/2+b
    }
    
    var_m=1/rgamma(1,a,as.numeric(b))
    }
    
    
    ##########################################
    ####### Digonal Coeffient ################
    ##########################################
    beta_ii_list=list(beta11,beta22,beta33)
    for(ii in 1:3){
      logtii=logtii_list[[ii]]
      betaii_all=NULL
      for(pp in 1:parts){
        if(sample(c(0,1),1,prob=c(0,1))==1){
          XVX=lapply(1:N,function(sub) +t(blockX_list[[sub]])%*%INV[[pp]]%*%blockX_list[[sub]]/(var_m/2))
          VAR=solve(Reduce("+",XVX)+mean_INV[[pp]]%x%diag(p)/s_beta)
          MEAN=Reduce("+",lapply(1:N,function(sub) VAR%*%t(blockX_list[[sub]])%*%(INV[[pp]]/(var_m/2))%*%logtii[partition_list[[pp]],sub]))
          betaii=matrix(rmvn(1,as.numeric(MEAN),matrix(VAR,p*NR,p*NR)),nrow=p)
          betaii_all=cbind(betaii_all,betaii)
        }else{
          betaii_all=cbind(betaii_all,beta_ii_list[[ii]][,partition_list[[pp]]])
        }
        
      }
      beta_ii_list[[ii]]=betaii_all
    }
    beta11=beta_ii_list[[1]]
    beta22=beta_ii_list[[2]]
    beta33=beta_ii_list[[3]]
    
    ##########################################
    ####### Off-Digonal Coeffient ############
    ##########################################
    ### beta21
    betakl_all=NULL
    for(pp in 1:parts){
      if(sample(c(0,1),1,prob=c(0,1))==1){
        INV_beta21=lapply(1:N, function(sub)  diag(residual22[[sub]][partition_list[[pp]]])%*% 
                            INV[[pp]]%*%
                            diag(residual22[[sub]][partition_list[[pp]]]) )
        X21=lapply(1:N, function(sub)  d11[[sub]][partition_list[[pp]]]*blockX_list[[sub]])
        XVX=lapply(1:N,function(sub) t(X21[[sub]])%*%INV_beta21[[sub]]%*%X21[[sub]]/(var_m))
        VAR=solve(Reduce("+",XVX)+mean_INV[[pp]]%x%diag(p)/s_beta)
        MEAN=Reduce("+",lapply(1:N,function(sub) VAR%*%t(X21[[sub]])%*%(INV_beta21[[sub]]/(var_m))%*%t21[partition_list[[pp]],sub]))
        betakl=matrix(rmvn(1,as.numeric(MEAN),matrix(VAR,p*NR,p*NR)),nrow=p)
        betakl_all=cbind(betakl_all,betakl)
      }else{
        betakl_all=cbind(betakl_all,beta21[,partition_list[[pp]]])
      }
    }
    beta21=betakl_all
    
    
    ### beta31
    betakl_all=NULL
    for(pp in 1:parts){
      if(sample(c(0,1),1,prob=c(0,1))==1){
        INV_beta31=lapply(1:N, function(sub)  diag(residual33[[sub]][partition_list[[pp]]])%*% 
                            INV[[pp]]%*%
                            diag(residual33[[sub]][partition_list[[pp]]]) )
        X31=lapply(1:N, function(sub)  d11[[sub]][partition_list[[pp]]]*blockX_list[[sub]])
        XVX=lapply(1:N,function(sub) t(X31[[sub]])%*%INV_beta31[[sub]]%*%X31[[sub]]/(var_m))
        VAR=solve(Reduce("+",XVX)++mean_INV[[pp]]%x%diag(p)/s_beta)
        MEAN=Reduce("+",lapply(1:N,function(sub) VAR%*%t(X31[[sub]])%*%(INV_beta31[[sub]]/(var_m))%*%t31[partition_list[[pp]],sub]))
        betakl=matrix(rmvn(1,as.numeric(MEAN),matrix(VAR,p*NR,p*NR)),nrow=p)
        betakl_all=cbind(betakl_all,betakl)
      }else{
        betakl_all=cbind(betakl_all,beta31[,partition_list[[pp]]])
      }
    }
    beta31=betakl_all
    
    ### beta32
    betakl_all=NULL
    for(pp in 1:parts){
      if(sample(c(0,1),1,prob=c(0,1))==1){
        INV_beta32=lapply(1:N, function(sub)  diag(residual33[[sub]][partition_list[[pp]]])%*% 
                            INV[[pp]]%*%
                            diag(residual33[[sub]][partition_list[[pp]]]) )
        X32=lapply(1:N, function(sub)  d22[[sub]][partition_list[[pp]]]*blockX_list[[sub]])
        XVX=lapply(1:N,function(sub) t(X32[[sub]])%*%INV_beta32[[sub]]%*%X32[[sub]]/(var_m))
        VAR=solve(Reduce("+",XVX)+mean_INV[[pp]]%x%diag(p)/s_beta)
        MEAN=Reduce("+",lapply(1:N,function(sub) VAR%*%t(X32[[sub]])%*%(INV_beta32[[sub]]/(var_m))%*%t32[partition_list[[pp]],sub]))
        betakl=matrix(rmvn(1,as.numeric(MEAN),matrix(VAR,p*NR,p*NR)),nrow=p)
        betakl_all=cbind(betakl_all,betakl)
      }else{
        betakl_all=cbind(betakl_all,beta32[,partition_list[[pp]]])
      }
    }
    beta32=betakl_all
    
    
    
    
    ##########################################
    #######Finalizing#########################
    ##########################################
    
    MCMC_rho[it]=rho
    MCMC_smooth[it]=smooth
    MCMC_var_m[it]=var_m
    
    MCMC_mean_rho[it]=mean_rho
    MCMC_mean_smooth[it]=mean_smooth
    MCMC_s_beta[it]=s_beta
    
    MCMC_beta$"beta11"[[it]]=beta11
    MCMC_beta$"beta22"[[it]]=beta22
    MCMC_beta$"beta33"[[it]]=beta33
    MCMC_beta$"beta21"[[it]]=beta21
    MCMC_beta$"beta31"[[it]]=beta31
    MCMC_beta$"beta32"[[it]]=beta32
    
    
    
    
    
  }
  
  return(list("MCMC_rho"=MCMC_rho,"MCMC_smooth"=MCMC_smooth,
              "MCMC_var_m"=MCMC_var_m,
              "MCMC_mean_rho"=MCMC_mean_rho,
              "MCMC_mean_smooth"=MCMC_mean_smooth,
              "MCMC_s_beta"=MCMC_s_beta,
              "MCMC_beta"=MCMC_beta)
         )
  
  
}


