## evaluate the approximate likelihood for error-affected case

rm(list=ls())
#install.packages("ggplot2","xtable")
library(ggplot2)
library(xtable)
set.seed(1)
setwd("C:/Users/trant/Dropbox/Thien_Phuc_Tran/software/simulation")

source('functions.r')

n <- c(10, 20) #c(10, 20, 50) 
tau2 <- c(0.1, 0.5, 1) 
n.rep <- 1000

bias.se.sd.conv <- data.frame(matrix(NA, ncol=4*length(n), nrow=8*length(tau2))) ## table of bias, se, sd and convergence
emp.cov.prob <- data.frame(matrix(NA, ncol=2*length(n), nrow=5*length(tau2))) ## table of emperical coverage probability and convergence

for(j in 1:length(n)){
  for(k in 1:length(tau2)){
    #load(paste('mean_dep_xinorm_n',n[j],'_tau',tau2[k],'_LIK.RData', sep='')) ##subgroup summary are available
    #load(paste('mean_dep_xinorm_n',n[j],'_tau',tau2[k],'_pLIK.RData', sep='')) ##subgroup summary are unavailable
    #load(paste('mean_dep_xinorm_n',n[j],'_tau',tau2[k],'_NAIVE.RData', sep='')) ##subgroup summary are unavailable
    
    #load(paste('mean_dep_xiskewnorm_n',n[j],'_tau',tau2[k],'_LIK.RData', sep='')) ##subgroup summary are available
    #load(paste('mean_dep_xiskewnorm_n',n[j],'_tau',tau2[k],'_pLIK.RData', sep='')) ##subgroup summary are available
    #load(paste('mean_dep_xiskewnorm_n',n[j],'_tau',tau2[k],'_NAIVE.RData', sep='')) ##subgroup summary are available
    
    #load(paste('gender_dep_xinorm_n',n[j],'_tau',tau2[k],'_LIK.RData', sep='')) ##subgroup summary are available
    #load(paste('gender_dep_xinorm_n',n[j],'_tau',tau2[k],'_pLIK.RData', sep='')) ##subgroup summary are available
    #load(paste('gender_dep_xinorm_n',n[j],'_tau',tau2[k],'_NAIVE.RData', sep='')) ##subgroup summary are available
    
    #load(paste('gender_dep_xiskewnorm_n',n[j],'_tau',tau2[k],'_LIK.RData', sep='')) ##subgroup summary are available
    #load(paste('gender_dep_xiskewnorm_n',n[j],'_tau',tau2[k],'_pLIK.RData', sep='')) ##subgroup summary are available
    #load(paste('gender_dep_xiskewnorm_n',n[j],'_tau',tau2[k],'_NAIVE.RData', sep='')) ##subgroup summary are available
    
    # bias.se.sd.conv[((k-1)*8+1):(k*8),((j-1)*4+1):(j*4)] <- evaluate(ris.lik.est, ris.lik.se, value = c(beta0,beta1,beta2,mu.xi,mu.z,tau2[k],sigma2.xi,sigma2.z))[,1:4]
    # emp.cov.prob[((k-1)*5+1):(k*5),((j-1)*2+1):(j*2)] <- emp.coverage(ris.lik.est[,1:5], ris.lik.se[,1:5], value = c(beta0,beta1,beta2,mu.xi,mu.z))
    
    # bias.se.sd.conv[((k-1)*8+1):(k*8),((j-1)*4+1):(j*4)] <- evaluate(ris.lik.est, sand.lik.se, value = c(beta0,beta1,beta2,mu.xi,mu.z,tau2[k],sigma2.xi,sigma2.z))
    # emp.cov.prob[((k-1)*5+1):(k*5),((j-1)*2+1):(j*2)] <- emp.coverage(ris.lik.est[,1:5], sand.lik.se[,1:5], value = c(beta0,beta1,beta2,mu.xi,mu.z))
    
    # bias.se.sd.conv[((k-1)*8+1):(k*8),((j-1)*4+1):(j*4)] <- evaluate(ris.lik.pseudo.approx.est, ris.lik.pseudo.approx.se, value = c(beta0,beta1,beta2,mu.xi,mu.z,tau2[k],sigma2.xi,sigma2.z))[,1:4]
    # emp.cov.prob[((k-1)*5+1):(k*5),((j-1)*2+1):(j*2)] <- emp.coverage(ris.lik.pseudo.approx.est[,1:5], ris.lik.pseudo.approx.se[,1:5], value = c(beta0,beta1,beta2,mu.xi,mu.z))
    
    M<-matrix(NA,nrow=2,ncol=4)
    colnames(M)<-c('bias','se','sd','convergence')
    bias.se.sd.conv[((k-1)*8+1):(k*8),((j-1)*4+1):(j*4)] <- rbind(evaluate(cbind(coef.naive[,-4],NA,NA,coef.naive[,4]), cbind(se.naive[,-4],NA,NA,se.naive[,4]), value = c(beta0,beta1,beta2,mu.xi,mu.z,tau2[k]))[,1:4],M)
    M2<-matrix(NA,nrow=2,ncol=2)
    colnames(M2)<-c('coverage','convergence')
    emp.cov.prob[((k-1)*5+1):(k*5),((j-1)*2+1):(j*2)] <- rbind(emp.coverage(coef.naive[,-4], se.naive[,-4], value = c(beta0,beta1,beta2)),M2)
    
  } 
}

n <- c(10, 20)
parameter <- rep(c("beta0", "beta1", "beta2", "mu.xi", "mu.z", "tau2", "sigma2.xi", "sigma2.z"), length(tau2)) ## extra covariates

bias.se.sd.conv <- cbind(parameter, bias.se.sd.conv)
colnames(bias.se.sd.conv) <- c("parameter", rep(c("bias", "se", "sd", "convergence"), length(n)))
tab <- print(xtable(bias.se.sd.conv,digits = c(0,0,rep(c(3,3,3,0),2))), include.rownames=F) ## latex table

emp.cov.prob <- cbind(parameter[1:5], emp.cov.prob) 
colnames(emp.cov.prob) <- c("parameter", rep(c("coverage", "convergence"), length(n)))

## plots of empirical coverage probability  
emp.cov.prob2 <- list(NULL) 
v <- emp.cov.prob[,-1]
for (i in 1:5) {
  u <- data.frame(matrix(NA,nrow=length(n)*length(tau2),ncol=3))
  for (j in 1:length(n)) {
    for (k in 1:length(tau2)) {
      u[((k-1)*length(n)+j),2] <- as.character(n[j])
      u[((k-1)*length(n)+j),1] <- tau2[k]
      u[((k-1)*length(n)+j),3] <- v[((k-1)*5+i),((j-1)*2+1)]
    }
  }
  emp.cov.prob2[[i]] <- u
  colnames(emp.cov.prob2[[i]]) <- c("tau2", "n", "coverage")
}
names(emp.cov.prob2) <- parameter[1:5]

plot.beta0 <- ggplot(emp.cov.prob2$beta0, aes(x=tau2, y=coverage, shape=n)) + geom_point() + 
  geom_hline(yintercept=0.95, linetype="dashed", color = "black", size=0.5) + 
  labs(title = bquote(beta[0] == ~ .(beta0))) + xlab(bquote(tau^{2}))
plot.beta1 <- ggplot(emp.cov.prob2$beta1, aes(x=tau2, y=coverage, shape=n)) + geom_point() + 
  geom_hline(yintercept=0.95, linetype="dashed", color = "black", size=0.5) + 
  labs(title = bquote(beta[1] == ~ .(beta1))) + xlab(bquote(tau^{2}))
plot.beta2 <- ggplot(emp.cov.prob2$beta2, aes(x=tau2, y=coverage, shape=n)) + geom_point() + 
  geom_hline(yintercept=0.95, linetype="dashed", color = "black", size=0.5) + 
  labs(title = bquote(beta[2] == ~ .(beta2))) + xlab(bquote(tau^{2}))
plot.mu.xi <- ggplot(emp.cov.prob2$mu.xi, aes(x=tau2, y=coverage, shape=n)) + geom_point() + 
  geom_hline(yintercept=0.95, linetype="dashed", color = "black", size=0.5) + 
  labs(title = bquote(mu[xi] == ~ .(mu.xi))) + xlab(bquote(tau^{2}))
plot.mu.z <- ggplot(emp.cov.prob2$mu.z, aes(x=tau2, y=coverage, shape=n)) + geom_point() + 
  geom_hline(yintercept=0.95, linetype="dashed", color = "black", size=0.5) + 
  labs(title = bquote(mu[zeta] == ~ .(mu.z))) + xlab(bquote(tau^{2}))

#save.image('mean_dep_xinorm_lik_RESULT.RData') ##subgroup are available
#save.image('mean_dep_xinorm_plik_RESULT.RData') ##subgroup are available
#save.image('mean_dep_xinorm_naive_RESULT.RData') ##subgroup are available

#save.image('mean_dep_xiskewnorm_lik_RESULT.RData') ##subgroup are available
#save.image('mean_dep_xiskewnorm_plik_RESULT.RData') ##subgroup are available
#save.image('mean_dep_xiskewnorm_naive_RESULT.RData') ##subgroup are available

#save.image('gender_dep_xinorm_lik_RESULT.RData') ##subgroup summary are available
#save.image('gender_dep_xinorm_plik_RESULT.RData') ##subgroup summary are available
#save.image('gender_dep_xinorm_naive_RESULT.RData') ##subgroup summary are available

#save.image('gender_dep_xiskewnorm_lik_RESULT.RData') ##subgroup summary are available
#save.image('gender_dep_xiskewnorm_plik_RESULT.RData') ##subgroup summary are available
#save.image('gender_dep_xiskewnorm_naive_RESULT.RData') ##subgroup summary are available
