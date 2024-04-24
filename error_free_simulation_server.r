## approximate likelihood for extra error-free covariates in server

rm(list=ls())
#library(mvtnorm, parallel)
# library(mvtnorm, lib.loc = "~/R_local_lib") ## run on serverpd
# library(parallel, lib.loc = "~/R_local_lib")
# library(nlme, lib.loc = "~/R_local_lib")
# library(statmod, lib.loc = "~/R_local_lib")
library(mvtnorm) ## run on serverpd
library(parallel)
library(nlme)
library(statmod)
set.seed(1)
#setwd("C:/Users/trant/Desktop/University of Padova/PhD Thesis/software")

source('functions.r')

n <- c(10, 20, 50)
beta0 <- 0
beta1 <- 1
tau2 <- c(0.1, 0.5, 1)
n.rep <- 1000

## skewnorm
# loc.xi <- 0 ## location
# scale.xi <- 1 ## scale
# shape.xi <- 5 ## shape
# mu.xi<-loc.xi+scale.xi*(shape.xi/sqrt(1+shape.xi^2))*sqrt(2/pi)
# sigma2.xi <- scale.xi^2*(1-(2/pi)*(shape.xi/sqrt(1+shape.xi^2))^2)

for(j in 1:2){
  for(k in 1:3){
    #load(paste('latitude_xinorm_n',n[j],'_tau',tau2[k],'.RData', sep=''))
    #load(paste('latitude_xiskewnorm_n',n[j],'_tau',tau2[k],'.RData', sep=''))
    #load(paste('severity_xinorm_n',n[j],'_tau',tau2[k],'.RData', sep=''))
    #load(paste('severity_xiskewnorm_n',n[j],'_tau',tau2[k],'.RData', sep=''))
    
    # ris.lik <- mclapply(1:n.rep, function(b, dati.all.=dati.all) {
    #   dati <- dati.all.[[b]]
    #   nj <- nrow(dati)
    #   m.hat <- lm(y.obs ~ x.obs + zi, data=dati) ## LS
    #   coef.naive <- c(coef(m.hat), mean(resid(m.hat)^2))
    #   se.naive <- c(sqrt(diag(vcov(m.hat))), (sqrt( 2*(nj-3)*coef.naive[4]^2/nj^2 )))
    #   
    #   # theta.start <- c(coef(m.hat),                 ## beta0, beta1, beta2
    #   #                  mean(dati[,2]),               ## mu.xi
    #   #                  mean(resid(m.hat)^2),   ## (tau^2)
    #   #                  sd(dati[,2])^2) ## (sigma.xi^2)
    #   theta.start <- c(coef(m.hat),                 ## beta0, beta1, beta2
    #                    mean(dati[,2]),               ## mu.xi
    #                    log(mean(resid(m.hat)^2)),   ## log.tau2
    #                    log(sd(dati[,2])^2)) ## log.sigma2.xi
    #   
    #   est <- se <- sand.se <- rep(NA, 6)
    #   conv <- NA
    #   
    #   ## approximate normal model
    #   model <- try(optim(theta.start, error.free.lik.approx.repa, dati=dati, control=list(fnscale=-1, maxit=5000), hessian=TRUE), silent=TRUE)
    #   # if (class(model)=='try-error'){
    #   #   model <- try(optim(theta.start, error.free.lik.approx, dati=dati, control=list(fnscale=-1, maxit=5000), hessian=TRUE, method="BFGS"), silent=TRUE)
    #   # }
    #   # if (class(model)=='try-error'){
    #   #   model <- try(optim(c(theta.start[-c(5,6)],theta.start[5]*0.5,theta.start[6]), error.free.lik.approx, dati=dati, control=list(fnscale=-1, maxit=5000), hessian=TRUE), silent=TRUE)
    #   # }
    #   
    #   if (class(model)!='try-error'){
    #     est <- model$par
    #     #se <- sqrt(diag(solve(-model$hessian)))
    #     conv <- model$convergence
    #     
    #     #H <- 
    #     G <- matrix(0,length(est),length(est))
    #     for(i in 1:nrow(dati)){
    #       a <- fdHess(est, error.free.lik.approx.repa, dati=dati[i,]) ## compute Hessian matrix
    #       #H <- H + a$Hessian
    #       values.gradient <- a$gradient
    #       G <- G + values.gradient%*%t(values.gradient) ## compute J matrix
    #     }
    #     invH <- try(solve(model$hessian), silent=TRUE)
    #     if (length(class(invH))>1){
    #       # se <- sqrt(diag(-invH))
    #       # sand.se <- sqrt(diag(invH%*%G%*%invH))
    #       se <- sqrt(diag(-diag(c(rep(1,4),exp(est[5:6])))%*%invH%*%diag(c(rep(1,4),exp(est[5:6])))))
    #       sand.se <- sqrt(diag(diag(c(rep(1,4),exp(est[5:6])))%*%invH%*%G%*%invH%*%diag(c(rep(1,4),exp(est[5:6])))))
    #     }
    #     est[5:6]<-exp(est[5:6])
    #   }
    #   
    #   # sand.cov <- try(sand.lik(est, error.affect.lik.approx, dati=dati), silent=TRUE)
    #   # if (class(sand.cov)!='try-error'){
    #   #   sand.se <- sqrt(diag(sand.cov))
    #   # }
    #   
    #   return(list(coef.naive=coef.naive,se.naive=se.naive,est=est, se=se, conv=conv, sand.se=sand.se))
    #   
    #   ## approximate normal model
    #   # model.normal <- try(optim(theta.start, error.affect.lik.approx, dati=dati, control=list(fnscale=-1, maxit=5000), hessian=TRUE), silent=TRUE)
    #   # #model.normal <- try(optim(theta.start, error.affect.lik.exact, dati=dati, n.node=10, control=list(fnscale=-1, maxit=5000), hessian=TRUE), silent=TRUE)
    #   # 
    #   # if (class(model.normal)!='try-error'){
    #   #   return(list(coef.naive=coef.naive,se.naive=se.naive,estimate=model.normal$par,se=sqrt(diag(solve(-model.normal$hessian))),convergence=model.normal$convergence))
    #   # }
    # }, mc.cores = 20) ## mc.cores specifies the number of cores to be used
    
    coef.naive <- se.naive <- matrix(NA, ncol=4, nrow=n.rep) ## coefficients of LS model
    
    # ris.lik.est <- matrix(NA, ncol=6, nrow=n.rep) ## estimates
    # ris.lik.se <- matrix(NA, ncol=6, nrow=n.rep)
    # ris.lik.sand.se <- matrix(NA, ncol=6, nrow=n.rep)
    # ris.lik.convergence <- rep(NA, length=n.rep)
    
    for(b in 1:n.rep){
      # if (is.null(ris.lik[[b]])==FALSE){
      #   ris.lik.est[b,] <- ris.lik[[b]]$est
      #   ris.lik.se[b,] <- ris.lik[[b]]$se
      #   ris.lik.sand.se[b,] <- ris.lik[[b]]$sand.se
      #   ris.lik.convergence[b] <- ris.lik[[b]]$conv
      # }
      
      dati <- dati.all[[b]]
      nj <- nrow(dati)
      w <- 1/var.eta
      
      ## error-free covariate
      m.hat <- lm(y.obs ~ x.obs+zi, data=dati, weights=w) ## weighted naive
      coef.naive[b,] <- c(coef(m.hat),                 ## beta0, beta1, beta2
                          mean(resid(m.hat)^2)) ## sigma2.eta
      se.naive[b,] <- c(sqrt(diag(vcov(m.hat))), (sqrt( 2*(nj-3)*coef.naive[b,4]^2/nj^2 )))
      
      
      # coef.naive[b,] <- c(coef(m.hat),                 ## beta0, beta1, beta2
      #                     mean(dati[,2]),               ## mu.xi
      #                     mean(resid(m.hat)^2)) ## sigma2.xi
      # se.naive[b,] <- c(sqrt(diag(vcov(m.hat))), sd(dati[,2])/sqrt(nj), (sqrt( 2*(nj-3)*coef.naive[b,5]^2/nj^2 )))
      
    }
    #save.image(paste('latitude_xinorm_n',n[j],'_tau',tau2[k],'_LIK.RData', sep='')) ##within-study varz is available
    #save.image(paste('latitude_xinorm_n',n[j],'_tau',tau2[k],'_NAIVE.RData', sep='')) ##within-study varz is available
    
    #save.image(paste('latitude_xiskewnorm_n',n[j],'_tau',tau2[k],'_LIK.RData', sep='')) ##within-study varz is available
    #save.image(paste('latitude_xiskewnorm_n',n[j],'_tau',tau2[k],'_NAIVE.RData', sep='')) ##within-study varz is available
    
    #save.image(paste('severity_xinorm_n',n[j],'_tau',tau2[k],'_LIK.RData', sep='')) ##within-study varz is available
    #save.image(paste('severity_xinorm_n',n[j],'_tau',tau2[k],'_NAIVE.RData', sep='')) ##within-study varz is available

    #save.image(paste('severity_xiskewnorm_n',n[j],'_tau',tau2[k],'_LIK.RData', sep='')) ##within-study varz is available
    #save.image(paste('severity_xiskewnorm_n',n[j],'_tau',tau2[k],'_NAIVE.RData', sep='')) ##within-study varz is available
    print(b)
    
  } 
}
