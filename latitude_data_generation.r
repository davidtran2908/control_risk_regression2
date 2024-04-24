## generation of data with Z measured without error - latitude

rm(list=ls())
library(mvtnorm,EnvStats,sn)
set.seed(1)
setwd("C:/Users/trant/Desktop/University of Padova/PhD Thesis/software")
#setwd("C:/Users/trant/Dropbox/Thien_Phuc_Tran/software")

source("C:/Users/trant/Dropbox/Thien_Phuc_Tran/software/functions.r")

n <- c(10, 20, 50) #10 is too small 
beta0 <- 0
beta1 <- 1
beta2 <- 0.8 ## associated to Z
tau2 <- c(0.1, 0.5, 1) #0.1 is too small ## tau2 is the variance of epsilon, sigma2 in our notes
tau <- sqrt(tau2)

mu.xi <- 0
sigma2.xi <- 1

# loc.xi <- 0 ## norm
# scale.xi <- 1
# skew.xi <- -5
# 
# mu.xi <- loc.xi+scale.xi*(skew.xi/sqrt(1+skew.xi^2))*sqrt(2/pi) ## centered skew-norm
# sigma2.xi <- scale.xi^2*(1-(2/pi)*(skew.xi/sqrt(1+skew.xi^2))^2)

#theta.true <- c(beta0, beta1, beta2, tau2, mu.xi, sigma2.xi)
#names(theta.true) <- c("beta0", "beta1", "beta2", "tau2", "mu.xi", "sigma2.xi")
n.rep <- 1000 ## B

for(j in 1:length(n)){
  ni.t <- round(runif(n[j], 15, 200)) ## number of treated 
  ni.c <- round(runif(n[j], 15, 200)) ## number of controls
  
  for(k in 1:length(tau2)){
    dati.all <- list(NULL)
    zi <- rnorm(n[j], 0, 1) ## more reasonable!!!
    for(b in 1:n.rep){
      xi <- rnorm(n[j], mu.xi, sqrt(sigma2.xi)) ## normal dist
      #xi <- rsn(n[j], loc.xi, scale.xi, skew.xi) ## skewnorm
      #xi <- xi-1*(5/sqrt(1+5^2))*sqrt(2/pi) ## centered skewnorm
      
      eta <- beta0 + beta1 * xi + beta2*zi  + rnorm(n[j], 0, sqrt(tau2[k])) ## insert Z here 
      p.mort.t <- plogis(eta) ## true.treatmentrisk = exp(eta)/(1+exp(eta))
      p.mort.c <- plogis(xi) ## true.controlrisk = exp(xi)/(1+exp(xi))
      y.i <- apply(matrix(cbind(ni.t, p.mort.t), ncol=2), 1, function(x) rbinom(1, x[1], x[2])) ## number of events in treatment group
      x.i <- apply(matrix(cbind(ni.c, p.mort.c), ncol=2), 1, function(x) rbinom(1, x[1], x[2])) ## number of events in control group
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
      dati <- data.frame(eta.hat=eta.hat, xi.hat=xi.hat, zi=zi,
                         var.eta=var.eta, cov.etaxi=0, cov.xieta=0, var.xi=var.xi, y.i=y.i, x.i=x.i, ni.t=ni.t, ni.c=ni.c)
      colnames(dati) <- c('y.obs', 'x.obs', 'zi', 'var.y', 'cov.yx', 'cov.yx', 'var.x', 'y', 'x', 'ni.t', 'ni.c')
      dati.all[[b]] <- dati
      
      #save.image(paste('latitude_xinorm_B100_n',n[j],'_tau',tau2[k],'.RData', sep='')) ## xi ~ norm 
      save.image(paste('latitude_xinorm_n',n[j],'_tau',tau2[k],'.RData', sep='')) ## xi ~ norm 
      #save.image(paste('latitude_xiskewnorm_n',n[j],'_tau',tau2[k],'.RData', sep='')) ## xi ~ skewnorm
      print(b) ## keep track of the loop
    } 
  }
}