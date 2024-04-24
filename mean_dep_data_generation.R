## generation of data with Z measured with error and aggregation bias

rm(list=ls())
library(mvtnorm,EnvStats,sn)
library(sn)
set.seed(1)
#set.seed(123) ## for n=15
setwd("C:/Users/trant/Desktop/University of Padova/PhD Thesis/software")
#setwd("C:/Users/trant/Dropbox/Thien_Phuc_Tran/software")

n <- c(10, 20) #15
beta0 <- 0
beta1 <- 1
beta2 <- 0.8 ## associated to Z
tau2 <- c(0.1, 0.5, 1) ## tau2 is the variance of epsilon, sigma2 in our notes

# mu.xi <- 0 ## norm
# sigma2.xi <- 1

loc.xi <- 0 ## norm
scale.xi <- 1
skew.xi <- -5

mu.xi <- loc.xi+scale.xi*(skew.xi/sqrt(1+skew.xi^2))*sqrt(2/pi) ## skewnorm
sigma2.xi <- scale.xi^2*(1-(2/pi)*(skew.xi/sqrt(1+skew.xi^2))^2)

tau <- sqrt(tau2)

mu.z <- 0
sigma2.z <- 1

#theta.true <- c(beta0, beta1, beta2, tau2, mu.xi, sigma2.xi, mu.z, sigma2.z)
#names(theta.true) <- c("beta0", "beta1", "beta2", "tau2", "mu.xi", "sigma2.xi", "mu.z", "sigma2.z")
var.zij <- 1 ## within-study variance
sd.zij <- sqrt(var.zij)

n.rep <- 1000 ## B

for(j in 1:length(n)){
  ni.t <- round(runif(n[j], 15, 200)) ## number of treated
  ni.c <- round(runif(n[j], 15, 200)) ## number of controls
  ni <- ni.t+ni.c
  
  rho.t <- runif(n[j], -1, 1) ## within-study correlation between ziT1 and ziT0
  rho.c <- runif(n[j], -1, 1) ## within-study correlation between ziC1 and ziC0
  for(k in 1:length(tau2)){
    dati.all <- list(NULL)
    
    for(b in 1:n.rep){
      zi <- rnorm(n[j], mu.z, sqrt(sigma2.z))
      
      #xi <- rnorm(n[j], mu.xi, sqrt(sigma2.xi)) ## norm
      xi <- rsn(n[j], loc.xi, scale.xi, skew.xi) ## skewnorm
      
      eta <- beta0 + beta1 * xi + beta2 * zi  + rnorm(n[j], 0, sqrt(tau2[k]))
      
      p.mort.t <- plogis(eta)
      p.mort.c <- plogis(xi)
      
      y.i <- apply(matrix(cbind(ni.t, p.mort.t), ncol=2), 1, function(x) rbinom(1, x[1], x[2])) ## events for the treated--> eta.hat
      x.i <- apply(matrix(cbind(ni.c, p.mort.c), ncol=2), 1, function(x) rbinom(1, x[1], x[2])) ## events for the conrols --> xi.hat
      
      ## subgroup summary and group summary
      zi.obs.mean.T1 <- zi.obs.mean.T0 <- zi.obs.mean.T <- rep(NA, length(y.i))
      for(i in 1:length(y.i)){
        if(y.i[i]==ni.t[i]) 
          zi.obs.mean.T1[i] <- zi.obs.mean.T[i] <- rnorm(1, zi[i], sd.zij/sqrt(ni.t[i]))
        
        if(y.i[i]==0) 
          zi.obs.mean.T0[i] <- zi.obs.mean.T[i] <- rnorm(1, zi[i], sd.zij/sqrt(ni.t[i]))
        
        x <- NA
        if(y.i[i]>0 & y.i[i]<ni.t[i]){
          x <- rmvnorm(1, c(zi[i],zi[i]), matrix(c(var.zij/y.i[i],rho.t[i]*sqrt((var.zij/y.i[i])*(var.zij/(ni.t[i]-y.i[i]))),rho.t[i]*sqrt((var.zij/y.i[i])*(var.zij/(ni.t[i]-y.i[i]))),var.zij/(ni.t[i]-y.i[i])), ncol=2))
          zi.obs.mean.T1[i] <- x[1]
          zi.obs.mean.T0[i] <- x[2]
          zi.obs.mean.T[i] <- (y.i[i]*zi.obs.mean.T1[i]+(ni.t[i]-y.i[i])*zi.obs.mean.T0[i])/ni.t[i]
        }
      }
      
      zi.obs.mean.C1 <- zi.obs.mean.C0 <- zi.obs.mean.C <- rep(NA, length(x.i))
      for(i in 1:length(x.i)){
        if(x.i[i]==ni.c[i]) 
          zi.obs.mean.C1[i] <- zi.obs.mean.C[i] <- rnorm(1, zi[i], sd.zij/sqrt(ni.c[i]))
        
        if(x.i[i]==0) 
          zi.obs.mean.C0[i] <- zi.obs.mean.C[i] <- rnorm(1, zi[i], sd.zij/sqrt(ni.c[i]))
        
        x <- NA
        if(x.i[i]>0 & x.i[i]<ni.c[i]){
          x <- rmvnorm(1, c(zi[i],zi[i]), matrix(c(var.zij/x.i[i],rho.c[i]*sqrt((var.zij/x.i[i])*(var.zij/(ni.c[i]-x.i[i]))),rho.c[i]*sqrt((var.zij/x.i[i])*(var.zij/(ni.c[i]-x.i[i]))),var.zij/(ni.c[i]-x.i[i])), ncol=2))
          zi.obs.mean.C1[i] <- x[1]
          zi.obs.mean.C0[i] <- x[2]
          zi.obs.mean.C[i] <- (x.i[i]*zi.obs.mean.C1[i]+(ni.c[i]-x.i[i])*zi.obs.mean.C0[i])/ni.c[i]
        }
      }
      
      ## estimate the mean for each study, the observed Z
      zi.mean <- (zi.obs.mean.T*ni.t + zi.obs.mean.C*ni.c)/ni
      ##  within variance 	
      zi.obs.var.T <- var.zij/ni.t
      zi.obs.var.C <- var.zij/ni.c
      zi.var <- var.zij/ni
      #zi.var <- (zi.obs.var.T*ni.t + zi.obs.var.C*ni.c)/(ni.t+ni.c)
      
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
      var.eta <- 1/y.i+1/(ni.t-y.i) 
      var.xi <- 1/x.i+1/(ni.c-x.i) 
      id.eta <- which(var.eta==Inf)
      if(length(id.eta)>0) ## check for infinite eta.hat variance and correct them
        var.eta[id.eta] <- 1/(y.i[id.eta]+0.5)+1/(ni.t[id.eta]-y.i[id.eta]+0.5)
      id.xi <- which(var.xi==Inf)
      if(length(id.xi)>0) ## check for infinite xi.hat variance and correct them
        var.xi[id.xi] <- 1/(x.i[id.xi]+0.5)+1/(ni.c[id.xi]-x.i[id.xi]+0.5)
      
      ## now compute the covariances at the within-study level
      ## subgroup information available
      cov.obs.etazi <- (zi.obs.mean.T1 - zi.obs.mean.T0)/ni
      cov.obs.xizi <- (zi.obs.mean.C1 - zi.obs.mean.C0)/ni
      cov.obs.etazi[is.na(cov.obs.etazi) | cov.obs.etazi>sqrt(var.eta*zi.var)] <- 0.0 ##when the presence of one single obs does not allow to compute the variance
      cov.obs.xizi[is.na(cov.obs.xizi) | cov.obs.xizi>sqrt(var.xi*zi.var)] <- 0.0 ##when the presence of one single obs does not allow to compute the variance
      
      dati <- data.frame(eta.hat=eta.hat, xi.hat=xi.hat, zi.obs=zi.mean, var.eta=var.eta,
                         cov.etaxi=0, cov.etaz=cov.obs.etazi, cov.xieta=0, var.xi=var.xi, cov.xiz=cov.obs.xizi,
                         cov.zeta=cov.obs.etazi, cov.zxi=cov.obs.xizi, var.z=zi.var, y.i=y.i, x.i=x.i, zi=zi,
                         ni.t=ni.t, ni.c=ni.c, mean.zi.T=zi.obs.mean.T, mean.zi.C=zi.obs.mean.C, var.zi.T=zi.obs.var.T, var.zi.C=zi.obs.var.C)
      # dati <- data.frame(eta.hat=eta.hat, xi.hat=xi.hat, zi.obs=zi.mean, var.eta=var.eta,
      #                    cov.etaxi=0, cov.etaz=0, cov.xieta=0, var.xi=var.xi, cov.xiz=0,
      #                    cov.zeta=0, cov.zxi=0, var.z=zi.var, y.i=y.i, x.i=x.i, zi=zi)
      # dati <- data.frame(eta.hat=eta.hat, xi.hat=xi.hat, zi.obs=zi.mean, var.eta=var.eta,
      #                    cov.etaxi=0, cov.etaz=0, cov.xieta=0, var.xi=var.xi, cov.xiz=0,
      #                    cov.zeta=0, cov.zxi=0, var.z=0, y.i=y.i, x.i=x.i, zi=zi)
      
      colnames(dati) <- c('y.obs', 'x.obs', 'zi.obs', 'var.y', 'cov.yx', 'cov.yz',
                          'cov.yx', 'var.x', 'cov.xz', 'cov.yz', 'cov.xz', 'var.z', 'y', 'x', 'zi', 'ni.t', 'ni.c', 'mean.zi.T', 'mean.zi.C', 'var.zi.T', 'var.zi.C')
      
      dati.all[[b]] <- dati
      
      #save.image(paste('mean_dep_xinorm_n',n[j],'_tau',tau2[k],'.RData', sep='')) ##without subgroup summary
      save.image(paste('mean_dep_xiskewnorm_n',n[j],'_tau',tau2[k],'.RData', sep='')) ##without subgroup summary
      print(b) ## keep track of the loop
    } 
  }
}