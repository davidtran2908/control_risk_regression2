## approximate likelihood for extra error-affected covariates in server

rm(list=ls())
# library(mvtnorm, lib.loc = "~/R_local_lib") ## run on serverpd
# library(parallel, lib.loc = "~/R_local_lib")
# library(nlme, lib.loc = "~/R_local_lib")
# library(statmod, lib.loc = "~/R_local_lib")
library(mvtnorm) ## run on serverpd
library(parallel)
library(nlme)
library(statmod)
set.seed(1)
#setwd("C:/Users/trant/Dropbox/Thien_Phuc_Tran/software")

source('functions.r')

n <- c(10, 20, 50) #15
beta0 <- 0
beta1 <- 1
tau2 <- c(0.1, 0.5, 1)
n.rep <- 1000

## skewnorm
# loc.xi <- 0 ## location
# scale.xi <- 1 ## scale
# shape.xi <- 3 ## shape
# mu.xi<-loc.xi+scale.xi*(shape.xi/sqrt(1+shape.xi^2))*sqrt(2/pi)
# sigma2.xi <- scale.xi^2*(1-(2/pi)*(shape.xi/sqrt(1+shape.xi^2))^2)

for(j in 1:2){
  for(k in 1:3){
    #load(paste('mean_dep_xinorm_n',n[j],'_tau',tau2[k],'.RData', sep=''))
    #load(paste('mean_dep_xiskewnorm_n',n[j],'_tau',tau2[k],'.RData', sep=''))
    
    #load(paste('gender_dep_xinorm_n',n[j],'_tau',tau2[k],'.RData', sep=''))
    #load(paste('gender_dep_xiskewnorm_n',n[j],'_tau',tau2[k],'.RData', sep=''))
    
    # ris.lik <- mclapply(1:n.rep, function(b, dati.all.=dati.all) {
    #   dati <- dati.all.[[b]]
    #   nj <- nrow(dati)
    #   # dati$x.obs.cen<-scale(dati$x.obs, scale=FALSE)
    #   # dati$zi.obs.cen<-scale(dati$zi.obs, scale=FALSE)
    #   m.hat <- lm(y.obs ~ x.obs + zi.obs, data=dati) ## LS
    #   #m.hat <- lm(y.obs ~ x.obs.cen + zi.obs.cen, data=dati) ## LS
    #   coef.naive <- c(coef(m.hat), mean(resid(m.hat)^2))
    #   se.naive <- c(sqrt(diag(vcov(m.hat))), (sqrt( 2*(nj-3)*coef.naive[4]^2/nj^2 )))
    # 
    #   # theta.start <- c(coef(m.hat),                 ## beta0, beta1, beta2
    #   #                  mean(dati[,2]),               ## mu.xi
    #   #                  mean(dati[,3]),               ## mu.z
    #   #                  mean(resid(m.hat)^2),   ## (tau^2)
    #   #                  sd(dati[,2])^2, ## (sigma.xi^2)
    #   #                  sd(dati[,3])^2)       ## (sigma.z^2)
    #   theta.start <- c(coef(m.hat),                 ## beta0, beta1, beta2
    #                    mean(dati[,2]),               ## mu.xi
    #                    mean(dati[,3]),               ## mu.z
    #                    log(mean(resid(m.hat)^2)),   ## (tau^2)
    #                    log(sd(dati[,2])^2), ## (sigma.xi^2)
    #                    log(sd(dati[,3])^2))       ## (sigma.z^2)
    # 
    #   est <- se <- sand.se <- rep(NA, 8)
    #   conv <- NA
    # 
    #   ## approximate normal model
    #   model <- try(optim(theta.start, error.affect.lik.approx.repa, dati=dati, pseudo=FALSE, control=list(fnscale=-1, maxit=5000), hessian=TRUE), silent=TRUE)
    #   if (class(model)!='try-error'){
    #     est <- model$par
    #     conv <- model$convergence
    # 
    #     G <- matrix(0,length(est),length(est))
    #     for(i in 1:nrow(dati)){
    #       a <- fdHess(est, error.affect.lik.approx.repa, dati=dati[i,], pseudo=FALSE) ## compute Hessian matrix
    #       #H <- H + a$Hessian
    #       values.gradient <- a$gradient
    #       G <- G + values.gradient%*%t(values.gradient) ## compute J matrix
    #     }
    #     invH <- try(solve(model$hessian), silent=TRUE)
    #     if (length(class(invH))>1){
    #       # se <- sqrt(diag(-invH))
    #       # sand.se <- sqrt(diag(invH%*%G%*%invH))
    #       se <- sqrt(diag(-diag(c(rep(1,5),exp(est[6:8])))%*%invH%*%diag(c(rep(1,5),exp(est[6:8])))))
    #       sand.se <- sqrt(diag(diag(c(rep(1,5),exp(est[6:8])))%*%invH%*%G%*%invH%*%diag(c(rep(1,5),exp(est[6:8])))))
    #     }
    #     est[6:8]<-exp(est[6:8])
    #   }
    # 
    #   # sand.cov <- try(sand.lik(est, error.affect.lik.approx, dati=dati), silent=TRUE)
    #   # if (class(sand.cov)!='try-error'){
    #   #   sand.se <- sqrt(diag(sand.cov))
    #   # }
    # 
    #   return(list(coef.naive=coef.naive,se.naive=se.naive,est=est, se=se, conv=conv, sand.se=sand.se))
    # }, mc.cores = 10)
    
    # ris.lik.pseudo.approx <- mclapply(1:n.rep, function(b, dati.all.=dati.all) {
    #   dati <- dati.all.[[b]]
    #   nj <- nrow(dati)
    #   # dati$x.obs.cen<-scale(dati$x.obs, scale=FALSE)
    #   # dati$zi.obs.cen<-scale(dati$zi.obs, scale=FALSE)
    #   m.hat <- lm(y.obs ~ x.obs + zi.obs, data=dati) ## LS
    #   #m.hat <- lm(y.obs ~ x.obs.cen + zi.obs.cen, data=dati) ## LS
    #   coef.naive <- c(coef(m.hat), mean(resid(m.hat)^2))
    #   se.naive <- c(sqrt(diag(vcov(m.hat))), (sqrt( 2*(nj-3)*coef.naive[4]^2/nj^2 )))
    # 
    #   # theta.start <- c(coef(m.hat),                 ## beta0, beta1, beta2
    #   #                  mean(dati[,2]),               ## mu.xi
    #   #                  mean(dati[,3]),               ## mu.z
    #   #                  mean(resid(m.hat)^2),   ## (tau^2)
    #   #                  sd(dati[,2])^2, ## (sigma.xi^2)
    #   #                  sd(dati[,3])^2)       ## (sigma.z^2)
    #   theta.start <- c(coef(m.hat),                 ## beta0, beta1, beta2
    #                    mean(dati[,2]),               ## mu.xi
    #                    mean(dati[,3]),               ## mu.z
    #                    log(mean(resid(m.hat)^2)),   ## (tau^2)
    #                    log(sd(dati[,2])^2), ## (sigma.xi^2)
    #                    log(sd(dati[,3])^2))       ## (sigma.z^2)
    # 
    #   est <- se <- sand.se <- rep(NA, 8)
    #   conv <- NA
    # 
    #   ## approximate normal model
    #   model <- try(optim(theta.start, error.affect.lik.approx.repa, dati=dati, pseudo=TRUE, control=list(fnscale=-1, maxit=5000), hessian=TRUE), silent=TRUE)
    #   if (class(model)!='try-error'){
    #     est <- model$par
    #     conv <- model$convergence
    # 
    #     G <- matrix(0,length(est),length(est))
    #     for(i in 1:nrow(dati)){
    #       a <- fdHess(est, error.affect.lik.approx.repa, dati=dati[i,], pseudo=TRUE)
    #       #H <- H + a$Hessian
    #       values.gradient <- a$gradient
    #       G <- G + values.gradient%*%t(values.gradient) ## compute J matrix
    #     }
    #     invH <- try(solve(model$hessian), silent=TRUE)
    #     if (length(class(invH))>1){
    #       # se <- sqrt(diag(-invH))
    #       # sand.se <- sqrt(diag(invH%*%G%*%invH))
    #       se <- sqrt(diag(-diag(c(rep(1,5),exp(est[6:8])))%*%invH%*%diag(c(rep(1,5),exp(est[6:8])))))
    #       sand.se <- sqrt(diag(diag(c(rep(1,5),exp(est[6:8])))%*%invH%*%G%*%invH%*%diag(c(rep(1,5),exp(est[6:8])))))
    #     }
    #     est[6:8]<-exp(est[6:8])
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

    # ris.lik.pseudo.exact <- mclapply(1:n.rep, function(b, dati.all.=dati.all) {
    #   dati <- dati.all.[[b]]
    #   nj <- nrow(dati)
    #   m.hat <- lm(y.obs ~ x.obs + zi.obs, data=dati) ## LS
    #   coef.naive <- c(coef(m.hat), mean(resid(m.hat)^2))
    #   se.naive <- c(sqrt(diag(vcov(m.hat))), (sqrt( 2*(nj-3)*coef.naive[4]^2/nj^2 )))
    # 
    #   # theta.start <- c(coef(m.hat),                 ## beta0, beta1, beta2
    #   #                  mean(dati[,2]),               ## mu.xi
    #   #                  mean(dati[,3]),               ## mu.z
    #   #                  mean(resid(m.hat)^2),   ## (tau^2)
    #   #                  sd(dati[,2])^2, ## (sigma.xi^2)
    #   #                  sd(dati[,3])^2)       ## (sigma.z^2)
    #   theta.start <- c(coef(m.hat),                 ## beta0, beta1, beta2
    #                    mean(dati[,2]),               ## mu.xi
    #                    mean(dati[,3]),               ## mu.z
    #                    log(mean(resid(m.hat)^2)),   ## (tau^2)
    #                    log(sd(dati[,2])^2), ## (sigma.xi^2)
    #                    log(sd(dati[,3])^2))       ## (sigma.z^2)
    # 
    #   est <- se <- sand.se <- rep(NA, 8)
    #   conv <- NA
    # 
    #   ## approximate normal model
    #   # model2 <- try(optim(theta.start, error.affect.lik.approx.repa, dati=dati, pseudo=TRUE, control=list(fnscale=-1, maxit=5000), hessian=TRUE), silent=TRUE)
    #   # if (class(model2)!='try-error'){
    #   #   theta.start <- model2$par
    #   # }
    #   model <- try(optim(theta.start, error.affect.lik.exact.repa, dati=dati, n.node=20, cov.type='mean', control=list(fnscale=-1, maxit=5000), hessian=TRUE), silent=TRUE)
    #   if (class(model)!='try-error'){
    #     est <- model$par
    #     conv <- model$convergence
    # 
    #     G <- matrix(0,length(est),length(est))
    #     for(i in 1:nrow(dati)){
    #       a <- fdHess(est, error.affect.lik.exact.repa, dati=dati[i,], n.node=20, cov.type='mean') ## compute Hessian matrix
    #       #H <- H + a$Hessian
    #       values.gradient <- a$gradient
    #       G <- G + values.gradient%*%t(values.gradient) ## compute J matrix
    #     }
    #     invH <- try(solve(model$hessian), silent=TRUE)
    #     if (length(class(invH))>1){
    #       # se <- sqrt(diag(-invH))
    #       # sand.se <- sqrt(diag(invH%*%G%*%invH))
    #       se <- sqrt(diag(-diag(c(rep(1,5),exp(est[6:8])))%*%invH%*%diag(c(rep(1,5),exp(est[6:8])))))
    #       sand.se <- sqrt(diag(diag(c(rep(1,5),exp(est[6:8])))%*%invH%*%G%*%invH%*%diag(c(rep(1,5),exp(est[6:8])))))
    #     }
    #     est[6:8]<-exp(est[6:8])
    #   }
    # 
    #   # sand.cov <- try(sand.lik(est, error.affect.lik.approx, dati=dati), silent=TRUE)
    #   # if (class(sand.cov)!='try-error'){
    #   #   sand.se <- sqrt(diag(sand.cov))
    #   # }
    # 
    #   return(list(coef.naive=coef.naive,se.naive=se.naive,est=est, se=se, conv=conv, sand.se=sand.se))
    # 
    # }, mc.cores = 20) ## mc.cores specifies the number of cores to be used
    
    coef.naive <- se.naive <- matrix(NA, ncol=4, nrow=n.rep) ## coefficients of LS model
    
    # ris.lik.est <- matrix(NA, ncol=8, nrow=n.rep) ## estimates
    # ris.lik.se <- matrix(NA, ncol=8, nrow=n.rep)
    # ris.lik.sand.se <- matrix(NA, ncol=8, nrow=n.rep)
    # ris.lik.convergence <- rep(NA, length=n.rep)
    
    # ris.lik.pseudo.approx.est <- matrix(NA, ncol=8, nrow=n.rep) ## estimates
    # ris.lik.pseudo.approx.se <- matrix(NA, ncol=8, nrow=n.rep)
    # ris.lik.pseudo.approx.sand.se <- matrix(NA, ncol=8, nrow=n.rep)
    # ris.lik.pseudo.approx.convergence <- rep(NA, length=n.rep)
    
    # ris.lik.pseudo.exact.est <- matrix(NA, ncol=8, nrow=n.rep) ## estimates
    # ris.lik.pseudo.exact.se <- matrix(NA, ncol=8, nrow=n.rep)
    # ris.lik.pseudo.exact.sand.se <- matrix(NA, ncol=8, nrow=n.rep)
    # ris.lik.pseudo.exact.convergence <- rep(NA, length=n.rep)
    
    for(b in 1:n.rep){
      # if (is.null(ris.lik[[b]])==FALSE){
      #   ris.lik.se[b,] <- ris.lik[[b]]$se
      #   ris.lik.sand.se[b,] <- ris.lik[[b]]$sand.se
      #   ris.lik.convergence[b] <- ris.lik[[b]]$conv
      # }
      
      # if (is.null(ris.lik.pseudo.approx[[b]])==FALSE){
      #   ris.lik.pseudo.approx.est[b,] <- ris.lik.pseudo.approx[[b]]$est
      #   ris.lik.pseudo.approx.se[b,] <- ris.lik.pseudo.approx[[b]]$se
      #   ris.lik.pseudo.approx.sand.se[b,] <- ris.lik.pseudo.approx[[b]]$sand.se
      #   ris.lik.pseudo.approx.convergence[b] <- ris.lik.pseudo.approx[[b]]$conv
      # }
      
      # if (is.null(ris.lik.pseudo.exact[[b]])==FALSE){
      #   ris.lik.pseudo.exact.est[b,] <- ris.lik.pseudo.exact[[b]]$est
      #   ris.lik.pseudo.exact.se[b,] <- ris.lik.pseudo.exact[[b]]$se
      #   ris.lik.pseudo.exact.sand.se[b,] <- ris.lik.pseudo.exact[[b]]$sand.se
      #   ris.lik.pseudo.exact.convergence[b] <- ris.lik.pseudo.exact[[b]]$conv
      # }
      
      dati <- dati.all[[b]]
      nj <- nrow(dati)
      w <- 1/var.eta
      
      ## error-prone covariate
      m.hat <- lm(y.obs ~ x.obs+zi.obs, data=dati, weights=w) ## weighted naive
      coef.naive[b,] <- c(coef(m.hat),                 ## beta0, beta1, beta2
                          mean(resid(m.hat)^2)) ## sigma2.eta
      se.naive[b,] <- c(sqrt(diag(vcov(m.hat))), (sqrt( 2*(nj-3)*coef.naive[b,4]^2/nj^2 )))
    }
    #save.image(paste('mean_dep_xinorm_n',n[j],'_tau',tau2[k],'_peLIK.RData', sep='')) ##subgroup summary are available 
    
    #save.image(paste('mean_dep_xinorm_n',n[j],'_tau',tau2[k],'_NAIVE.RData', sep='')) ##subgroup summary are available 
    #save.image(paste('mean_dep_xiskewnorm_n',n[j],'_tau',tau2[k],'_NAIVE.RData', sep='')) ##subgroup summary are available 
    
    #save.image(paste('gender_dep_xinorm_n',n[j],'_tau',tau2[k],'_NAIVE.RData', sep='')) ##subgroup summary are available 
    #save.image(paste('gender_dep_xiskewnorm_n',n[j],'_tau',tau2[k],'_NAIVE.RData', sep='')) ##subgroup summary are available 
    
    print(b)
  } 
}
