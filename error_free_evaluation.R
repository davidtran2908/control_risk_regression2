## evaluate the approximate likelihood for error-free case

rm(list=ls())
#install.packages("ggplot2","xtable")
library(ggplot2)
library(xtable)
set.seed(1)

source('C:/Users/trant/Dropbox/Thien_Phuc_Tran/software/simulation/functions.r')

n <- c(10, 20) #c(10, 20, 50) #10 is too small 
tau2 <- c(0.1, 0.5, 1) #0.1 is too small
#n.rep <- 1000

bias.se.sd.conv <- data.frame(matrix(NA, ncol=4*length(n), nrow=6*length(tau2))) ## table of bias, se, sd and convergence
emp.cov.prob <- data.frame(matrix(NA, ncol=2*length(n), nrow=4*length(tau2))) ## table of emperical coverage probability and convergence

for(j in 1:2){
  for(k in 1:3){
    #load(paste('severity_xinorm_n',n[j],'_tau',tau2[k],'_LIK.RData', sep='')) ## severity 
    #load(paste('severity_xinorm_n',n[j],'_tau',tau2[k],'_NAIVE.RData', sep='')) ## severity 

    #load(paste('severity_xiskewnorm_n',n[j],'_tau',tau2[k],'_LIK.RData', sep='')) ## severity 
    #load(paste('severity_xiskewnorm_n',n[j],'_tau',tau2[k],'_NAIVE.RData', sep='')) ## severity 
    
    #load(paste('latitude_xinorm_n',n[j],'_tau',tau2[k],'_LIK.RData', sep='')) ## latitude
    #load(paste('latitude_xinorm_n',n[j],'_tau',tau2[k],'_NAIVE.RData', sep='')) ## latitude

    #load(paste('latitude_xiskewnorm_n',n[j],'_tau',tau2[k],'_LIK.RData', sep='')) ## LATITUDE
    #load(paste('latitude_xiskewnorm_n',n[j],'_tau',tau2[k],'_NAIVE.RData', sep='')) ## LATITUDE
    
    # bias.se.sd.conv[((k-1)*6+1):(k*6),((j-1)*4+1):(j*4)] <- evaluate(ris.lik.est, ris.lik.se, value = c(beta0,beta1,beta2,mu.xi,tau2[k],sigma2.xi))[,1:4]
    # emp.cov.prob[((k-1)*4+1):(k*4),((j-1)*2+1):(j*2)] <- emp.coverage(ris.lik.est[,1:4], ris.lik.se[,1:4], value = c(beta0,beta1,beta2,mu.xi))

    M<-matrix(NA,nrow=1,ncol=4)
    colnames(M)<-c('bias','se','sd','convergence')
    bias.se.sd.conv[((k-1)*6+1):(k*6),((j-1)*4+1):(j*4)] <- rbind(evaluate(cbind(coef.naive[,-4],NA,coef.naive[,4]), cbind(se.naive[,-4],NA,se.naive[,4]), value = c(beta0,beta1,beta2,mu.xi,tau2[k]))[,1:4],M)
    emp.cov.prob[((k-1)*4+1):(k*4),((j-1)*2+1):(j*2)] <- rbind(emp.coverage(coef.naive[,-4], se.naive[,-4], value = c(beta0,beta1,beta2)),rep(NA,2))
  } 
}

n <- c(10, 20) #10 is too small 

parameter <- rep(c("beta0", "beta1", "beta2", "mu.xi", "tau2", "sigma2.xi"), length(tau2)) ## extra covariates

bias.se.sd.conv <- cbind(parameter, bias.se.sd.conv)
colnames(bias.se.sd.conv) <- c("parameter", rep(c("bias", "se", "sd", "convergence"), length(n)))
tab <- print(xtable(bias.se.sd.conv,digits = c(0,0,rep(c(3,3,3,0),2))), include.rownames=F) ## latex table

emp.cov.prob <- cbind(parameter[1:4], emp.cov.prob) 
colnames(emp.cov.prob) <- c("parameter", rep(c("coverage", "convergence"), length(n)))

## plots of empirical coverage probability  
emp.cov.prob2 <- list(NULL) 
v <- emp.cov.prob[,-1]
for (i in 1:4) {
  u <- data.frame(matrix(NA,nrow=length(n)*length(tau2),ncol=3))
  for (j in 1:length(n)) {
    for (k in 1:length(tau2)) {
      u[((k-1)*length(n)+j),2] <- as.character(n[j])
      u[((k-1)*length(n)+j),1] <- tau2[k]
      u[((k-1)*length(n)+j),3] <- v[((k-1)*4+i),((j-1)*2+1)]
    }
  }
  emp.cov.prob2[[i]] <- u
  colnames(emp.cov.prob2[[i]]) <- c("tau2", "n", "coverage")
}
names(emp.cov.prob2) <- parameter[1:4]

plot.beta0 <- ggplot(emp.cov.prob2$beta0) + geom_point(aes(x=tau2, y=coverage, size=n)) + 
  geom_hline(yintercept=0.95, linetype="dashed", color = "black", size=0.5) + 
  labs(title = bquote(beta[0] == ~ .(beta0))) + xlab(bquote(tau^{2})) + 
  scale_size_manual(values=c(4,7)) + theme(legend.position = "none")
plot.beta1 <- ggplot(emp.cov.prob2$beta1) + geom_point(aes(x=tau2, y=coverage, size=n)) + 
  geom_hline(yintercept=0.95, linetype="dashed", color = "black", size=0.5) + 
  labs(title = bquote(beta[1] == ~ .(beta1))) + xlab(bquote(tau^{2})) + 
  scale_size_manual(values=c(4,7)) + theme(legend.position = "none")
plot.beta2 <- ggplot(emp.cov.prob2$beta2) + geom_point(aes(x=tau2, y=coverage, size=n)) + 
  geom_hline(yintercept=0.95, linetype="dashed", color = "black", size=0.5) + 
  labs(title = bquote(beta[2] == ~ .(beta2))) + xlab(bquote(tau^{2})) + 
  scale_size_manual(values=c(4,7)) + theme(legend.position = "none")
plot.mu.xi <- ggplot(emp.cov.prob2$mu.xi) + geom_point(aes(x=tau2, y=coverage, size=n)) +
  geom_hline(yintercept=0.95, linetype="dashed", color = "black", size=0.5) +
  labs(title = bquote(mu[xi] == ~ .(mu.xi))) + xlab(bquote(tau^{2})) +
  scale_size_manual(values=c(4,7))

#save.image('severity_xinorm_RESULT.RData') ## severity
#save.image('severity_xinorm_naive_RESULT.RData') ## severity

#save.image('severity_xiskewnorm_RESULT.RData') ## severity
#save.image('severity_xiskewnorm_naive_RESULT.RData') ## severity

#save.image('latitude_xinorm_RESULT.RData') ## latitude
#save.image('latitude_xinorm_naive_RESULT.RData') ## latitude

#save.image('latitude_xiskewnorm_RESULT.RData') ## latitude
#save.image('latitude_xiskewnorm_naive_RESULT.RData') ## latitude
