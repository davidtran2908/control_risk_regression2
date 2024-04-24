## generation of data with gender  

rm(list=ls())
library(mvtnorm,EnvStats)
library(sn)
set.seed(1)
setwd("C:/Users/trant/Desktop/University of Padova/PhD Thesis/software")
#setwd("C:/Users/trant/Dropbox/Thien_Phuc_Tran/software")

n <- c(10, 20) 
beta0 <- 0
beta1 <- 1
beta2 <- 0.8 ## associated to Z
tau2 <- c(0.1, 0.5, 1) ## tau2 is the variance of epsilon, sigma2 in our notes

# mu.xi <- 0 ## norm
# sigma2.xi <- 1

loc.xi <- 0 ## skewnorm
scale.xi <- 1
skew.xi <- -5

mu.xi <- loc.xi+scale.xi*(skew.xi/sqrt(1+skew.xi^2))*sqrt(2/pi) ## centered skew-norm
sigma2.xi <- scale.xi^2*(1-(2/pi)*(skew.xi/sqrt(1+skew.xi^2))^2)

tau <- sqrt(tau2)

mu.z <- 0
sigma2.z <- 1 
#theta.true <- c(beta0, beta1, beta2, tau2, mu.xi, sigma2.xi, mu.z, sigma2.z)
#names(theta.true) <- c("beta0", "beta1", "beta2", "tau2", "mu.xi", "sigma2.xi", "mu.z", "sigma2.z")

n.rep <- 1000 ## B

for(j in 1:length(n)){
  ni.t <- round(runif(n[j], 15, 200)) ## number of treated 
  ni.c <- round(runif(n[j], 15, 200)) ## number of controls
  ni <- ni.t+ni.c
  
  # rho.yz<-runif(n[j], -1, 1)
  # rho.xz<-runif(n[j], -1, 1)
  
  for(k in 1:length(tau2)){
    dati.all <- list(NULL)
    for(b in 1:n.rep){
      zi <- rnorm(n[j], mu.z, sqrt(sigma2.z))
      p.male <- plogis(zi)
      
      #xi <- rnorm(n[j], mu.xi, sqrt(sigma2.xi)) ## normal dist
      xi <- rsn(n[j], loc.xi, scale.xi, skew.xi) ## skewnorm
      
      eta <- beta0 + beta1 * xi + beta2*zi + rnorm(n[j], 0, sqrt(tau2[k])) ## insert Z here 
      p.mort.t <- plogis(eta) ## true.treatmentrisk = exp(eta)/(1+exp(eta))
      
      ## generate joint probability 
      p11.t <- apply(matrix(cbind(p.mort.t, p.male), ncol=2), 1, function(x) runif(1, 0, min(x[1], x[2])))
      p10.t <- p.mort.t - p11.t
      p01.t <- p.male - p11.t
      p00.t <- 1 - (p11.t + p10.t + p01.t)
      repeat{ ## correct negative probability
        p11.t[which(p00.t<0)] <- apply(matrix(cbind(p.mort.t[which(p00.t<0)], p.male[which(p00.t<0)]), ncol=2), 1, function(x) runif(1, 0, min(x[1], x[2])))
        p10.t <- p.mort.t - p11.t
        p01.t <- p.male - p11.t
        p00.t <- 1 - (p11.t + p10.t + p01.t)
        if(min(p00.t)>=0){
          break
        }
      }
      # p00.t<-(1-p.mort.t)*(1-p.male)+rho.yz*sqrt(p.mort.t*(1-p.mort.t)*p.male*(1-p.male))
      # p10.t<-1-p.male-p00.t
      # p01.t<-1-p.mort.t-p00.t
      # p11.t<-1-(p00.t+p10.t+p01.t)
      
      p.mort.c <- plogis(xi) ## true.controlrisk = exp(xi)/(1+exp(xi))
      p11.c <- apply(matrix(cbind(p.mort.c, p.male), ncol=2), 1, function(x) runif(1, 0, min(x[1], x[2])))
      p10.c <- p.mort.c - p11.c
      p01.c <- p.male - p11.c
      p00.c <- 1 - (p11.c + p10.c + p01.c)
      repeat{
        p11.c[which(p00.c<0)] <- apply(matrix(cbind(p.mort.c[which(p00.c<0)], p.male[which(p00.c<0)]), ncol=2), 1, function(x) runif(1, 0, min(x[1], x[2])))
        p10.c <- p.mort.c - p11.c
        p01.c <- p.male - p11.c
        p00.c <- 1 - (p11.c + p10.c + p01.c)
        if(min(p00.c)>=0){
          break
        }
      }
      # p00.c<-(1-p.mort.c)*(1-p.male)+rho.xz*sqrt(p.mort.c*(1-p.mort.c)*p.male*(1-p.male))
      # p10.c<-1-p.male-p00.c
      # p01.c<-1-p.mort.c-p00.c
      # p11.c<-1-(p00.c+p10.c+p01.c)
      
      ## subgroup summary and group summary
      subgr.t <- t(apply(matrix(cbind(ni.t,p11.t,p10.t,p01.t,p00.t), ncol=5), 1, function(x) rmultinom(1, x[1], x[-1]))) ## treatment group
      subgr.c <- t(apply(matrix(cbind(ni.c,p11.c,p10.c,p01.c,p00.c), ncol=5), 1, function(x) rmultinom(1, x[1], x[-1]))) ## control group
      
      ## number of patients in subgroups (1: events, 0: no events)
      male.i.treated1 <- subgr.t[,1]
      male.i.treated0 <- subgr.t[,3] ## treated + male + no events
      female.i.treated1 <- subgr.t[,2] ## treated + female + events
      female.i.treated0 <- subgr.t[,4] ## treated + female + no events
      male.i.treated <- male.i.treated1+male.i.treated0 ## number of males in treatment group
      female.i.treated <- ni.t-male.i.treated ## number of females in treatment group
      
      male.i.control1 <- subgr.c[,1] ## control + male + events
      male.i.control0 <- subgr.c[,3] ## control + male + no events
      female.i.control1 <- subgr.c[,2] ## control + female + events
      female.i.control0 <- subgr.c[,4] ## control + female + no events
      male.i.control <- male.i.control1+male.i.control0 ## number of males in control group
      female.i.control <- ni.c-male.i.control ## number of females in the control group
      
      y.i <- male.i.treated1+female.i.treated1 ## number of events in treatment group
      x.i <- male.i.control1+female.i.control1 ## number of events in control group
      
      male.i <- male.i.treated+male.i.control ## number of males in the i-th study
      female.i <- ni-male.i ## number of females in the i-th study
      
      ## observed measures
      eta.hat <- log( (y.i)/(ni.t-y.i) )
      #w.new <- 1/( 1/(y.i) + 1/(ni.t-y.i) ) ## inverse sd of eta.hat
      id.eta <- which(eta.hat==Inf | eta.hat==-Inf) 
      if(length(id.eta)>0){ ## check for infinite eta.hat and correct them
        eta.hat[id.eta] <- log( (y.i[id.eta]+0.5)/(ni.t[id.eta]-y.i[id.eta]+0.5) )
        #w.new[id.eta] <- 1/( 1/(y.i[id.eta]+0.5) + 1/(ni.t[id.eta]-y.i[id.eta]+0.5) )
      }
      xi.hat <- log( (x.i)/(ni.c-x.i) )
      id.xi <- which(xi.hat==Inf | xi.hat==-Inf)
      if(length(id.xi)>0) ## check for infinite xi.hat and correct them
        xi.hat[id.xi] <- log( (x.i[id.xi]+0.5)/(ni.c[id.xi]-x.i[id.xi]+0.5) )
      
      zi.obs <- log( (male.i)/(ni-male.i) )
      id.zi <- which(zi.obs==Inf | zi.obs==-Inf)
      if(length(id.zi)>0) ## check for infinite zi.obs and correct them
        zi.obs[id.zi] <- log( (male.i[id.zi]+0.5)/(ni[id.zi]-male.i[id.zi]+0.5) )
      
      ## within-study variance
      var.eta <- 1/y.i+1/(ni.t-y.i) 
      var.xi <- 1/x.i+1/(ni.c-x.i)
      
      var.z <- 1/male.i+1/(ni-male.i)
      
      id.eta <- which(var.eta==Inf)
      if(length(id.eta)>0) ## check for infinite eta.hat variance and correct them
        var.eta[id.eta] <- 1/(y.i[id.eta]+0.5)+1/(ni.t[id.eta]-y.i[id.eta]+0.5)
      id.xi <- which(var.xi==Inf)
      if(length(id.xi)>0) ## check for infinite xi.hat variance and correct them
        var.xi[id.xi] <- 1/(x.i[id.xi]+0.5)+1/(ni.c[id.xi]-x.i[id.xi]+0.5)
      
      id.zi <- which(var.z==Inf)
      if(length(id.zi)>0) ## check for infinite zi.obs variance and correct them
        var.z[id.zi] <- 1/(male.i[id.zi]+0.5)+1/(ni[id.zi]-male.i[id.zi]+0.5)
      
      ## now compute covariances at the within-study level
      cov.obs.etazi <- (male.i.treated1/y.i-male.i.treated0/(ni.t-y.i))/male.i-(female.i.treated1/y.i-female.i.treated0/(ni.t-y.i))/female.i
      cov.obs.etazi[is.na(cov.obs.etazi)] <- 0.0 ##when the presence of one single obs does not allow to compute the variance
      cov.obs.xizi <- (male.i.control1/x.i-male.i.control0/(ni.c-x.i))/male.i-(female.i.control1/x.i-female.i.control0/(ni.c-x.i))/female.i
      cov.obs.xizi[is.na(cov.obs.xizi)] <- 0.0 ##when the presence of one single obs does not allow to compute the variance
      
      dati <- data.frame(eta.hat=eta.hat, xi.hat=xi.hat, zi.obs=zi.obs, var.eta=var.eta,
                         cov.etaxi=0, cov.etaz=cov.obs.etazi, cov.xieta=0, var.xi=var.xi, cov.xiz=cov.obs.xizi,
                         cov.zeta=cov.obs.etazi, cov.zxi=cov.obs.xizi, var.z=var.z, y.i=y.i, x.i=x.i, zi=zi, ni.t=ni.t, ni.c=ni.c, male.i.treated=male.i.treated, male.i.control=male.i.control)
      
      # dati <- data.frame(eta.hat=eta.hat, xi.hat=xi.hat, zi.obs=zi.obs, var.eta=var.eta,
      #                    cov.etaxi=0, cov.etaz=0, cov.xieta=0, var.xi=var.xi, cov.xiz=0,
      #                    cov.zeta=0, cov.zxi=0, var.z=var.z, y.i=y.i, x.i=x.i, zi=zi, male.i.treated=male.i.treated, male.i.control=male.i.control)
       
      colnames(dati) <- c('y.obs', 'x.obs', 'zi.obs', 'var.y', 'cov.yx', 'cov.yz',
                          'cov.yx', 'var.x', 'cov.xz', 'cov.yz', 'cov.xz', 'var.z', 'y', 'x', 'zi', 'ni.t', 'ni.c', 'male.treated', 'male.control')
      
      dati.all[[b]] <- dati
      
      #save.image(paste('gender_dep_xinorm_n',n[j],'_tau',tau2[k],'.RData', sep='')) 
      save.image(paste('gender_dep_xiskewnorm_n',n[j],'_tau',tau2[k],'.RData', sep='')) 
      
      print(b) ## keep track of the loop
    } 
  }
}