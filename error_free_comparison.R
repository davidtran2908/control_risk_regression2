## compare likelihoods 

rm(list=ls())
#install.packages("ggplot2","xtable")
library(ggplot2)
library(xtable)
set.seed(1)
setwd("C:/Users/trant/Dropbox/Thien_Phuc_Tran/software/simulation")

#load('latitude_xinorm_RESULT.RData') ## latitude
#load('latitude_xiskewnorm_RESULT.RData') ## latitude
#load('severity_xinorm_RESULT.RData') ## severity
#load('severity_xiskewnorm_RESULT.RData') ## severity

eva <- bias.se.sd.conv
cover <- emp.cov.prob2

#load('latitude_xinorm_naive_RESULT.RData') ## latitude
#load('latitude_xiskewnorm_naive_RESULT.RData') ## latitude
#load('severity_xinorm_naive_RESULT.RData') ## severity
#load('severity_xiskewnorm_naive_RESULT.RData') ## severity

eva2 <- bias.se.sd.conv
cover2 <- emp.cov.prob2

eva3 <- data.frame(matrix(NA, ncol=ncol(eva), nrow=2*nrow(eva)))
for (i in 1:nrow(eva)) {
  eva3[2*i-1,] <- eva[i,]
  eva3[2*i,] <- eva2[i,]
}
eva3 <- data.frame(eva3[,1], approach=rep(c("lik", "naive"), nrow(eva)), eva3[,-1])
colnames(eva3) <- c("parameter", "covariance", rep(c("bias", "se", "sd", "convergence"), length(n)))
tab.compare <- print(xtable(eva3, digits = c(0,0,0,rep(c(3,3,3,0),2))), include.rownames=F) ## latex table

cover3 <- list(NULL)
for (s in 1:4) { 
  a <- data.frame(matrix(NA, ncol=ncol(cover[[s]]), nrow=2*nrow(cover[[s]]))) ## table of emperical coverage probability and convergence
  for (i in 1:nrow(cover[[s]])) {
    a[2*i-1,] <- cover[[s]][i,]
    a[2*i,] <- cover2[[s]][i,]
  }
  cover3[[s]] <- data.frame(a[,1:2], approach=rep(c("lik", "naive"), nrow(cover[[s]])), a[,3])
  colnames(cover3[[s]]) <- c("tau2", "n", "approach", "coverage")
}
names(cover3) <- c("beta0", "beta1", "beta2", "mu.xi")

plot.compare.beta0 <- ggplot(data=cover3$beta0) + geom_point(aes(x=tau2, y=coverage, size=n, shape=approach)) +
  geom_hline(yintercept=0.95, linetype="dashed", color = "black", size=0.5) + labs(title = bquote(beta[0] == ~ .(beta0))) + xlab(bquote(tau^{2})) +
  scale_size_manual(values=c(4,7))+ scale_shape_manual(values=c(19, 0)) + theme(legend.position = "none")
plot.compare.beta1 <- ggplot(data=cover3$beta1) + geom_point(aes(x=tau2, y=coverage, size=n, shape=approach)) +
  geom_hline(yintercept=0.95, linetype="dashed", color = "black", size=0.5) + labs(title = bquote(beta[1] == ~ .(beta1))) + xlab(bquote(tau^{2})) +
  scale_size_manual(values=c(4,7))+ scale_shape_manual(values=c(19, 0)) + theme(legend.position = "none")
plot.compare.beta2 <- ggplot(data=cover3$beta2) + geom_point(aes(x=tau2, y=coverage, size=n, shape=approach)) +
  geom_hline(yintercept=0.95, linetype="dashed", color = "black", size=0.5) + labs(title = bquote(beta[2] == ~ .(beta2))) + xlab(bquote(tau^{2})) +
  scale_size_manual(values=c(4,7))+ scale_shape_manual(values=c(19, 0)) + theme(legend.position = "none")
plot.compare.mu.xi <- ggplot(data=cover3$mu.xi) + geom_point(aes(x=tau2, y=coverage, size=n, shape=approach)) +
  geom_hline(yintercept=0.95, linetype="dashed", color = "black", size=0.5) + labs(title = bquote(mu[xi] == ~ .(mu.xi))) + xlab(bquote(tau^{2})) +
  scale_size_manual(values=c(4,7))+ scale_shape_manual(values=c(19, 0))

#save.image('latitude_xinorm_COMPARE.RData') ## latitude
#save.image('latitude_xiskewnorm_COMPARE.RData') ## latitude
#save.image('severity_xinorm_COMPARE.RData') ## severity
#save.image('severity_xiskewnorm_COMPARE.RData') ## severity
