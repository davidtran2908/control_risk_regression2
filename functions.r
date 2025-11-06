
## approximate likelihood function
basic.lik.approx <- function(theta, dati, center=TRUE){
    n <- nrow(dati)
    beta0 <- theta[1]
    beta1 <- theta[2]
    mu.xi <- theta[3]
    sigma2.eta <- theta[4]
    # sigma2.eta <- theta[3]
    # mu.xi <- theta[4]
    sigma2.xi <-theta[5]
    if(theta[4] < 0 | theta[5] < 0){ ## check variances
        return(NA)
    }
    mean.vector <- c(beta0, mu.xi)
    if(center==FALSE){mean.vector <- c(beta0+beta1*mu.xi, mu.xi)}
    #S.matrix <- matrix(dati[,3:6], ncol=2, byrow=TRUE)
    lik <- 0.0
    for(i in 1:n){
        S.matrix <- matrix(c(dati[i,3]+(beta1^2)*sigma2.xi+sigma2.eta, dati[i,4]+beta1*sigma2.xi,
                             dati[i,5]+beta1*sigma2.xi, dati[i,6]+sigma2.xi), ncol=2, byrow=TRUE) ## S.matrix depends on i
        lik <- lik + dmvnorm(dati[i,1:2], mean = mean.vector, sigma = S.matrix, log=TRUE)
    }
    return(lik)
}

basic.lik.approx.repa <- function(theta, dati, center=TRUE){
  n <- nrow(dati)
  beta0 <- theta[1]
  beta1 <- theta[2]
  mu.xi <- theta[3]
  log.sigma2.eta <- theta[4]
  # sigma2.eta <- theta[3]
  # mu.xi <- theta[4]
  log.sigma2.xi <-theta[5]
  mean.vector <- c(beta0, mu.xi)
  if(center==FALSE){mean.vector <- c(beta0+beta1*mu.xi, mu.xi)}
  #S.matrix <- matrix(dati[,3:6], ncol=2, byrow=TRUE)
  lik <- 0.0
  for(i in 1:n){
    S.matrix <- matrix(c(dati[i,3]+(beta1^2)*exp(log.sigma2.xi)+exp(log.sigma2.eta), dati[i,4]+beta1*exp(log.sigma2.xi),
                         dati[i,5]+beta1*exp(log.sigma2.xi), dati[i,6]+exp(log.sigma2.xi)), ncol=2, byrow=TRUE) ## S.matrix depends on i
    lik <- lik + dmvnorm(dati[i,1:2], mean = mean.vector, sigma = S.matrix, log=TRUE)
  }
  return(lik)
}

## approximate likelihood function in case of error-free covariates
error.free.lik.approx <- function(theta, dati){
    n <- nrow(dati)
    beta0 <- theta[1]
    beta1 <- theta[2]
    beta2 <- theta[3]
    mu.xi <- theta[4]
    sigma2.eta <- theta[5]
    # sigma2.eta <- theta[3]
    # mu.xi <- theta[4]
    sigma2.xi <-theta[6]
    if(theta[5] < 0 | theta[6] < 0) ## check variances
        return(NA)
    #S.matrix <- matrix(dati[,3:6], ncol=2, byrow=TRUE)
    lik <- 0.0
    for(i in 1:n){
        mean.vector <- c(beta0+beta1*mu.xi+beta2*dati[i,3], mu.xi)
        S.matrix <- matrix(c(dati[i,4]+(beta1^2)*sigma2.xi+sigma2.eta, dati[i,5]+beta1*sigma2.xi,
                             dati[i,6]+beta1*sigma2.xi, dati[i,7]+sigma2.xi), ncol=2, byrow=TRUE) ## S.matrix depends on i
        lik <- lik + dmvnorm(dati[i,1:2], mean = mean.vector, sigma = S.matrix, log=TRUE)
    }
    return(lik)
}

error.free.lik.approx.repa <- function(theta, dati){
  n <- nrow(dati)
  beta0 <- theta[1]
  beta1 <- theta[2]
  beta2 <- theta[3]
  mu.xi <- theta[4]
  log.sigma2.eta <- theta[5]
  # sigma2.eta <- theta[3]
  # mu.xi <- theta[4]
  log.sigma2.xi <-theta[6]
  #S.matrix <- matrix(dati[,3:6], ncol=2, byrow=TRUE)
  lik <- 0.0
  for(i in 1:n){
    mean.vector <- c(beta0+beta1*mu.xi+beta2*dati[i,3], mu.xi)
    S.matrix <- matrix(c(dati[i,4]+(beta1^2)*exp(log.sigma2.xi)+exp(log.sigma2.eta), dati[i,5]+beta1*exp(log.sigma2.xi),
                         dati[i,6]+beta1*exp(log.sigma2.xi), dati[i,7]+exp(log.sigma2.xi)), ncol=2, byrow=TRUE) ## S.matrix depends on i
    lik <- lik + dmvnorm(dati[i,1:2], mean = mean.vector, sigma = S.matrix, log=TRUE)
  }
  return(lik)
}

error.free.lik.exact <- function(theta, dati, n.node=10){
  n <- nrow(dati)
  beta0 <- theta[1]
  beta1 <- theta[2]
  beta2 <- theta[3]
  mu.xi <- theta[4]
  sigma2.eta <- theta[5]
  sigma2.xi <-theta[6]
  if(theta[5] < 0 | theta[6] < 0){ ## check variances
    return(NA)
  }
  lik <- 0.0
  ## nodes and weights for Gauss-Hermite quadrature
  objGH <- gauss.quad(n.node, "hermite")
  w <- objGH$weights
  w <- cbind(rep(w,rep(n.node,n.node)),rep(w,n.node))
  w <- w[,1]*w[,2]
  node <- objGH$nodes
  node <- cbind(rep(node,rep(n.node,n.node)),rep(node,n.node))
  for(i in 1:n){
    g<-function(x){
      p.eta <- exp(x[1])/(1+exp(x[1]))
      p.xi <- exp(x[2])/(1+exp(x[2]))
      return(dbinom(dati[i,8],dati[i,10],p.eta)*dbinom(dati[i,9],dati[i,11],p.xi)*dnorm(x[1],beta0+beta1*x[2]+beta2*dati[i,3],sqrt(sigma2.eta))*dnorm(x[2],mu.xi,sqrt(sigma2.xi))*exp(x[1]^2+x[2]^2))
    }    
    h <- apply(node, MARGIN=1, FUN=g)
    lik <- lik + log(sum(w*h))
  }
  return(lik)
}

## approximate likelihood function in case of aggregate covariates
error.affect.lik.approx <- function(theta, dati){
  n <- nrow(dati)
  beta0 <- theta[1]
  beta1 <- theta[2]
  beta2 <- theta[3]
  mu.xi <- theta[4]
  sigma2.eta <- theta[6]
  sigma2.xi <-theta[7]
  
  mu.z <- theta[5]
  sigma2.z <- theta[8]
  
  if(theta[6] < 0 | theta[7] < 0 | theta[8] < 0) ## check variances
    return(NA)
  lik <- 0.0
  mean.vector <- c(beta0+beta1*mu.xi+beta2*mu.z, mu.xi, mu.z)
  #mean.vector <- c(beta0, mu.xi, mu.z)
  
  for(i in 1:n){
    S.matrix <- matrix(c(dati[i,4]+(beta1^2)*sigma2.xi+(beta2^2)*sigma2.z+sigma2.eta,
                         dati[i,5]+beta1*sigma2.xi, dati[i,6]+beta2*sigma2.z,
                         dati[i,7]+beta1*sigma2.xi, dati[i,8]+sigma2.xi, dati[i,9],
                         dati[i,10]+beta2*sigma2.z, dati[i,11], sigma2.z+dati[i,12]), ncol=3, byrow=TRUE) ## S.matrix depends on i
    
    lik <- lik + dmvnorm(dati[i,1:3], mean = mean.vector, sigma = S.matrix, log=TRUE)
  }
  return(lik)
}

error.affect.lik.approx.repa <- function(theta, dati, pseudo=TRUE){
  n <- nrow(dati)
  beta0 <- theta[1]
  beta1 <- theta[2]
  beta2 <- theta[3]
  mu.xi <- theta[4]
  log.sigma2.eta <- theta[6]
  log.sigma2.xi <-theta[7]
  
  mu.z <- theta[5]
  log.sigma2.z <- theta[8]
  
  lik <- 0.0
  mean.vector <- c(beta0+beta1*mu.xi+beta2*mu.z, mu.xi, mu.z)
  #mean.vector <- c(beta0, mu.xi, mu.z)
  
  for(i in 1:n){
    if(pseudo==TRUE)
      S.matrix <- matrix(c(dati[i,4]+(beta1^2)*exp(log.sigma2.xi)+(beta2^2)*exp(log.sigma2.z)+exp(log.sigma2.eta),
                           beta1*exp(log.sigma2.xi), beta2*exp(log.sigma2.z),
                           beta1*exp(log.sigma2.xi), dati[i,8]+exp(log.sigma2.xi), 0,
                           beta2*exp(log.sigma2.z), 0, exp(log.sigma2.z)+dati[i,12]), ncol=3, byrow=TRUE) ## S.matrix depends on i
    else
      S.matrix <- matrix(c(dati[i,4]+(beta1^2)*exp(log.sigma2.xi)+(beta2^2)*exp(log.sigma2.z)+exp(log.sigma2.eta),
                           dati[i,5]+beta1*exp(log.sigma2.xi), dati[i,6]+beta2*exp(log.sigma2.z),
                           dati[i,7]+beta1*exp(log.sigma2.xi), dati[i,8]+exp(log.sigma2.xi), dati[i,9],
                           dati[i,10]+beta2*exp(log.sigma2.z), dati[i,11], exp(log.sigma2.z)+dati[i,12]), ncol=3, byrow=TRUE) ## S.matrix depends on i
        
    lik <- lik + dmvnorm(dati[i,1:3], mean = mean.vector, sigma = S.matrix, log=TRUE)
  }
  return(lik)
}

error.affect.lik.exact.repa <- function(theta, dati, n.node, cov.type){
  n <- nrow(dati)
  beta0 <- theta[1]
  beta1 <- theta[2]
  beta2 <- theta[3]
  mu.xi <- theta[4]
  log.sigma2.eta <- theta[6]
  log.sigma2.xi <-theta[7]

  mu.z <- theta[5]
  log.sigma2.z <- theta[8]

  lik <- 0.0

  ## nodes and weights for Gauss-Hermite quadrature
  objGH <- gauss.quad(n.node, "hermite")
  w <- objGH$weights
  w <- cbind(rep(w,rep(n.node^2,n.node)),rep(rep(w,rep(n.node,n.node)),n.node),rep(w,n.node^2))
  w <- w[,1]*w[,2]*w[,3]
  node <- objGH$nodes
  node <- cbind(rep(node,rep(n.node^2,n.node)),rep(rep(node,rep(n.node,n.node)),n.node),rep(node,n.node^2))
  if(cov.type=='log.odds'){
    for(i in 1:n){
      g<-function(x){
        xi <- sqrt(2*dati[i,8])*x[2]+dati[i,2]
        eta <- sqrt(2*dati[i,4])*x[1]+dati[i,1]
        zeta <- sqrt(2*dati[i,12])*x[3]+dati[i,3]
        p.eta <- exp(eta)/(1+exp(eta))
        p.xi <- exp(xi)/(1+exp(xi))
        p.zeta <- exp(zeta)/(1+exp(zeta))
        return(dbinom(dati[i,13],dati[i,16],p.eta)*dbinom(dati[i,14],dati[i,17],p.xi)*dbinom(dati[i,18]+dati[i,19],dati[i,16]+dati[i,17],p.zeta)*dnorm(eta,beta0+beta1*xi+beta2*zeta,sqrt(exp(log.sigma2.eta)))*dnorm(xi,mu.xi,sqrt(exp(log.sigma2.xi)))*dnorm(zeta,mu.z,sqrt(exp(log.sigma2.z)))*exp(x[1]^2+x[2]^2+x[3]^2))
        
        # eta<-x[1]
        # xi<-x[2]
        # zeta<-x[3]
        # p.eta <- exp(eta)/(1+exp(eta))
        # p.xi <- exp(xi)/(1+exp(xi))
        # p.zeta <- exp(zeta)/(1+exp(zeta))
        # return(dbinom(dati[i,13],dati[i,16],p.eta)*dbinom(dati[i,14],dati[i,17],p.xi)*dbinom(dati[i,18]+dati[i,19],dati[i,16]+dati[i,17],p.zeta)*dnorm(eta,beta0+beta1*xi+beta2*zeta,sqrt(exp(log.sigma2.eta)))*dnorm(xi,mu.xi,sqrt(exp(log.sigma2.xi)))*dnorm(zeta,mu.z,sqrt(exp(log.sigma2.z)))*exp(x[1]^2+x[2]^2+x[3]^2))
      }    
      h <- apply(node, MARGIN=1, FUN=g)
      lik <- lik + log(sum(w*h))
    }  
  }
  if(cov.type=='mean'){
    for(i in 1:n){
      g<-function(x){
        xi <- sqrt(2*dati[i,8])*x[2]+dati[i,2]
        eta <- sqrt(2*dati[i,4])*x[1]+dati[i,1]
        zeta <- sqrt(2*dati[i,12])*x[3]+dati[i,3]
        p.eta <- exp(eta)/(1+exp(eta))
        p.xi <- exp(xi)/(1+exp(xi))
        return(dbinom(dati[i,13],dati[i,16],p.eta)*dbinom(dati[i,14],dati[i,17],p.xi)*dnorm(eta,beta0+beta1*xi+beta2*zeta,sqrt(exp(log.sigma2.eta)))*dnorm(xi,mu.xi,sqrt(exp(log.sigma2.xi)))*dnorm(zeta,mu.z,sqrt(exp(log.sigma2.z)))*exp(x[1]^2+x[2]^2))

        # eta<-x[1]
        # xi<-x[2]
        # zeta<-x[3]
        # p.eta <- exp(eta)/(1+exp(eta))
        # p.xi <- exp(xi)/(1+exp(xi))
        # p.zeta <- exp(zeta)/(1+exp(zeta))
        # return(dbinom(dati[i,13],dati[i,16],p.eta)*dbinom(dati[i,14],dati[i,17],p.xi)*dbinom(dati[i,18]+dati[i,19],dati[i,16]+dati[i,17],p.zeta)*dnorm(eta,beta0+beta1*xi+beta2*zeta,sqrt(exp(log.sigma2.eta)))*dnorm(xi,mu.xi,sqrt(exp(log.sigma2.xi)))*dnorm(zeta,mu.z,sqrt(exp(log.sigma2.z)))*exp(x[1]^2+x[2]^2+x[3]^2))
      }    
      h <- apply(node, MARGIN=1, FUN=g)
      lik <- lik + log(sum(w*h))
    }  
  }
  return(lik)
}

error.affect.lik.exact.repa2 <- function(theta, dati, n.node, cov.type){
  n <- nrow(dati)
  beta0 <- theta[1]
  beta1 <- theta[2]
  beta2 <- theta[3]
  mu.xi <- theta[4]
  log.sigma2.eta <- theta[6]
  log.sigma2.xi <-theta[7]
  
  mu.z <- theta[5]
  log.sigma2.z <- theta[8]
  
  lik <- 0.0
  
  ## nodes and weights for Gauss-Hermite quadrature
  objGH <- gauss.quad(n.node, "hermite")
  w <- objGH$weights
  w <- cbind(rep(w,rep(n.node^2,n.node)),rep(rep(w,rep(n.node,n.node)),n.node),rep(w,n.node^2))
  w <- w[,1]*w[,2]*w[,3]
  node <- objGH$nodes
  node <- cbind(rep(node,rep(n.node^2,n.node)),rep(rep(node,rep(n.node,n.node)),n.node),rep(node,n.node^2))
  if(cov.type=='log.odds'){
    for(i in 1:n){
      g<-function(x){
        xi <- x[2]
        eta <- x[1]
        zeta <- x[3]
        p.eta <- exp(eta)/(1+exp(eta))
        p.xi <- exp(xi)/(1+exp(xi))
        p.zeta <- exp(zeta)/(1+exp(zeta))
        return(dbinom(dati[i,13],dati[i,16],p.eta)*dbinom(dati[i,14],dati[i,17],p.xi)*dbinom(dati[i,18]+dati[i,19],dati[i,16]+dati[i,17],p.zeta)*dnorm(eta,beta0+beta1*xi+beta2*zeta,sqrt(exp(log.sigma2.eta)))*dnorm(xi,mu.xi,sqrt(exp(log.sigma2.xi)))*dnorm(zeta,mu.z,sqrt(exp(log.sigma2.z)))*exp(x[1]^2+x[2]^2+x[3]^2))
      }    
      h <- apply(node, MARGIN=1, FUN=g)
      lik <- lik + log(sum(w*h))
    }  
  }
  if(cov.type=='mean'){
    for(i in 1:n){
      g<-function(x){
        xi <- x[2]
        eta <- x[1]
        zeta <- x[3]
        # xi <- sqrt(2*(exp(log.sigma2.xi)+dati[i,8]))*x[2]+mu.xi
        # eta <- sqrt(2*(beta1^2*exp(log.sigma2.xi)+beta2^2*exp(log.sigma2.z)+exp(log.sigma2.eta)+dati[i,4]))*x[1]+beta0+beta1*mu.xi+beta2*mu.z
        # zeta <- sqrt(2*(exp(log.sigma2.z)+dati[i,12]))*x[3]+mu.z
        p.eta <- exp(eta)/(1+exp(eta))
        p.xi <- exp(xi)/(1+exp(xi))
        return(dbinom(dati[i,13],dati[i,16],p.eta)*dbinom(dati[i,14],dati[i,17],p.xi)*dnorm(dati[i,3],zeta,sqrt(dati[i,12]))*dnorm(eta,beta0+beta1*xi+beta2*zeta,sqrt(exp(log.sigma2.eta)))*dnorm(xi,mu.xi,sqrt(exp(log.sigma2.xi)))*dnorm(zeta,mu.z,sqrt(exp(log.sigma2.z)))*exp(x[1]^2+x[2]^2+x[3]^2))
        #return(dbinom(dati[i,13],dati[i,16],p.eta)*dbinom(dati[i,14],dati[i,17],p.xi)*dnorm(zeta,dati[i,3],sqrt(dati[i,12]))*dnorm(eta,beta0+beta1*xi+beta2*zeta,sqrt(exp(log.sigma2.eta)))*dnorm(xi,mu.xi,sqrt(exp(log.sigma2.xi)))*dnorm(zeta,mu.z,sqrt(exp(log.sigma2.z)))*exp(x[1]^2+x[2]^2+x[3]^2)*sqrt(2*(exp(log.sigma2.xi)+dati[i,8]))*sqrt(2*(beta1^2*exp(log.sigma2.xi)+beta2^2*exp(log.sigma2.z)+exp(log.sigma2.eta)+dati[i,4]))*sqrt(2*(exp(log.sigma2.z)+dati[i,12])))
      }    
      h <- apply(node, MARGIN=1, FUN=g)
      lik <- lik + log(sum(w*h))
    }  
  }
  return(lik)
}

## approximate likelihood function + GH approximation
basic.lik.GH <- function(theta, dati, n.node=10, model="approx"){
  n <- nrow(dati)
  beta0 <- theta[1]
  beta1 <- theta[2]
  mu.xi <- theta[3]
  sigma2.eta <- theta[4]
  sigma2.xi <-theta[5]
  if(theta[4] < 0 | theta[5] < 0){ ## check variances
    return(NA)
  }
  lik <- 0.0
  ## nodes and weights for Gauss-Hermite quadrature
  objGH <- gauss.quad(n.node, "hermite")
  w <- objGH$weights
  w <- cbind(rep(w,rep(n.node,n.node)),rep(w,n.node))
  w <- w[,1]*w[,2]
  node <- objGH$nodes
  node <- cbind(rep(node,rep(n.node,n.node)),rep(node,n.node))
  if (model == "approx"){
    for(i in 1:n){
      g <- function(x){dnorm(dati[i,1],x[1],sqrt(dati[i,3]))*dnorm(dati[i,2],x[2],sqrt(dati[i,6]))*dnorm(x[1],beta0+beta1*x[2],sqrt(sigma2.eta))*dnorm(x[2],mu.xi,sqrt(sigma2.xi))*exp(x[1]^2+x[2]^2)}    
      h <- apply(node, MARGIN=1, FUN=g)
      lik <- lik + log(sum(w*h))
    }
  } 
  if (model == "exact"){
    for(i in 1:n){
      g<-function(x){dbinom(dati[i,7],dati[i,9],exp(x[1])/(1+exp(x[1])))*dbinom(dati[i,8],dati[i,10],exp(x[2])/(1+exp(x[2])))*dnorm(x[1],beta0+beta1*x[2],sqrt(sigma2.eta))*dnorm(x[2],mu.xi,sqrt(sigma2.xi))*exp(x[1]^2+x[2]^2)}    
      h <- apply(node, MARGIN=1, FUN=g)
      lik <- lik + log(sum(w*h))
    }
  }
  return(lik)
}

basic.lik.GH3 <- function(theta, dati, n.node=10, model="approx", center=TRUE){
  n <- nrow(dati)
  beta0 <- theta[1]
  beta1 <- theta[2]
  mu.xi <- theta[3]
  sigma2.eta <- theta[4]
  sigma2.xi <-theta[5]
  if(theta[4] < 0 | theta[5] < 0){ ## check variances
    return(NA)
  }
  lik <- 0.0
  ## nodes and weights for Gauss-Hermite quadrature
  objGH <- gauss.quad(n.node, "hermite")
  w <- objGH$weights
  w <- cbind(rep(w,rep(n.node,n.node)),rep(w,n.node))
  w <- w[,1]*w[,2]
  node <- objGH$nodes
  node <- cbind(rep(node,rep(n.node,n.node)),rep(node,n.node))
  if (model == "approx"){
    for(i in 1:n){
      g <- function(x){
        xi <- sqrt(2*dati[i,6])*x[2]+dati[i,2]
        eta <- sqrt(2*dati[i,3])*x[1]+dati[i,1]
        return(dnorm(eta,beta0+beta1*(xi-mu.xi),sqrt(sigma2.eta))*dnorm(xi,mu.xi,sqrt(sigma2.xi)))
      }
      h <- apply(node, MARGIN=1, FUN=g)
      lik <- lik + log(sum(w*h))
    }
  } 
  if (model == "exact"){
    for(i in 1:n){
      g<-function(x){
        xi <- sqrt(2*dati[i,6])*x[2]+dati[i,2]
        eta <- sqrt(2*dati[i,3])*x[1]+dati[i,1]
        p.eta <- exp(eta)/(1+exp(eta))
        p.xi <- exp(xi)/(1+exp(xi))
        return(dbinom(dati[i,7],dati[i,9],p.eta)*dbinom(dati[i,8],dati[i,10],p.xi)*dnorm(eta,beta0+beta1*(xi-mu.xi),sqrt(sigma2.eta))*dnorm(xi,mu.xi,sqrt(sigma2.xi))*exp(x[1]^2+x[2]^2))
      }    
      h <- apply(node, MARGIN=1, FUN=g)
      lik <- lik + log(sum(w*h))
    }
  }
  if(center==FALSE){
    if (model == "approx"){
      for(i in 1:n){
        g <- function(x){
          xi <- sqrt(2*dati[i,6])*x[2]+dati[i,2]
          eta <- sqrt(2*dati[i,3])*x[1]+dati[i,1]
          return(dnorm(eta,beta0+beta1*(xi),sqrt(sigma2.eta))*dnorm(xi,mu.xi,sqrt(sigma2.xi)))
        }
        h <- apply(node, MARGIN=1, FUN=g)
        lik <- lik + log(sum(w*h))
      }
    } 
    if (model == "exact"){
      for(i in 1:n){
        g<-function(x){
          xi <- sqrt(2*dati[i,6])*x[2]+dati[i,2]
          eta <- sqrt(2*dati[i,3])*x[1]+dati[i,1]
          p.eta <- exp(eta)/(1+exp(eta))
          p.xi <- exp(xi)/(1+exp(xi))
          return(dbinom(dati[i,7],dati[i,9],p.eta)*dbinom(dati[i,8],dati[i,10],p.xi)*dnorm(eta,beta0+beta1*(xi),sqrt(sigma2.eta))*dnorm(xi,mu.xi,sqrt(sigma2.xi))*exp(x[1]^2+x[2]^2))
        }    
        h <- apply(node, MARGIN=1, FUN=g)
        lik <- lik + log(sum(w*h))
      }
    }
  }
  return(lik)
}

basic.lik.GH3.repa <- function(theta, dati, n.node=10, center=TRUE){
  n <- nrow(dati)
  beta0 <- theta[1]
  beta1 <- theta[2]
  mu.xi <- theta[3]
  log.sigma2.eta <- theta[4]
  log.sigma2.xi <-theta[5]
  
  lik <- 0.0
  ## nodes and weights for Gauss-Hermite quadrature
  objGH <- gauss.quad(n.node, "hermite")
  w <- objGH$weights
  w <- cbind(rep(w,rep(n.node,n.node)),rep(w,n.node))
  w <- w[,1]*w[,2]
  node <- objGH$nodes
  node <- cbind(rep(node,rep(n.node,n.node)),rep(node,n.node))
  for(i in 1:n){
    g<-function(x){
      xi <- sqrt(2*dati[i,6])*x[2]+dati[i,2]
      eta <- sqrt(2*dati[i,3])*x[1]+dati[i,1]
      p.eta <- exp(eta)/(1+exp(eta))
      p.xi <- exp(xi)/(1+exp(xi))
      if(center==TRUE)
        return(dbinom(dati[i,7],dati[i,9],p.eta)*dbinom(dati[i,8],dati[i,10],p.xi)*dnorm(eta,beta0+beta1*(xi-mu.xi),sqrt(exp(log.sigma2.eta)))*dnorm(xi,mu.xi,sqrt(exp(log.sigma2.xi)))*exp(x[1]^2+x[2]^2)*sqrt(2*dati[i,3]*2*dati[i,6]))
      else
        return(dbinom(dati[i,7],dati[i,9],p.eta)*dbinom(dati[i,8],dati[i,10],p.xi)*dnorm(eta,beta0+beta1*(xi),sqrt(exp(log.sigma2.eta)))*dnorm(xi,mu.xi,sqrt(exp(log.sigma2.xi)))*exp(x[1]^2+x[2]^2)*sqrt(2*dati[i,3]*2*dati[i,6]))
    }    
    h <- apply(node, MARGIN=1, FUN=g)
    lik <- lik + log(sum(w*h))
  }
  return(lik)
}

basic.lik.GH7.repa <- function(theta, dati, n.node=10, center=TRUE){
  n <- nrow(dati)
  beta0 <- theta[1]
  beta1 <- theta[2]
  mu.xi <- theta[3]
  log.sigma2.eta <- theta[4]
  log.sigma2.xi <-theta[5]
  
  lik <- 0.0
  ## nodes and weights for Gauss-Hermite quadrature
  objGH <- gauss.quad(n.node, "hermite")
  w <- objGH$weights
  w <- cbind(rep(w,rep(n.node,n.node)),rep(w,n.node))
  w <- w[,1]*w[,2]
  node <- objGH$nodes
  node <- cbind(rep(node,rep(n.node,n.node)),rep(node,n.node))
  for(i in 1:n){
    g<-function(x){
      xi <- sqrt(2*(dati[i,6]))*x[2]+dati[i,2]
      eta <- sqrt(2*exp(log.sigma2.eta))*x[1]+dati[i,1]
      p.eta <- exp(eta)/(1+exp(eta))
      p.xi <- exp(xi)/(1+exp(xi))
      if(center==TRUE)
        return(dbinom(dati[i,7],dati[i,9],p.eta)*dbinom(dati[i,8],dati[i,10],p.xi)*dnorm(eta,beta0+beta1*(xi-mu.xi),sqrt(exp(log.sigma2.eta)))*dnorm(xi,mu.xi,sqrt(exp(log.sigma2.xi)))*exp(x[1]^2+x[2]^2)*sqrt(2*exp(log.sigma2.eta)*2*(dati[i,6])))
      else
        return(dbinom(dati[i,7],dati[i,9],p.eta)*dbinom(dati[i,8],dati[i,10],p.xi)*dnorm(eta,beta0+beta1*(xi),sqrt(exp(log.sigma2.eta)))*dnorm(xi,mu.xi,sqrt(exp(log.sigma2.xi)))*exp(x[1]^2+x[2]^2)*sqrt(2*exp(log.sigma2.eta)*2*(dati[i,6])))
    }    
    h <- apply(node, MARGIN=1, FUN=g)
    lik <- lik + log(sum(w*h))
  }
  return(lik)
}

## approximate likelihood function + GH approximation for quadratic relationship

quad.lik.GH2 <- function(theta, dati, n.node=10, model="approx"){
  n <- nrow(dati)
  beta0 <- theta[1]
  beta1 <- theta[2]
  beta2 <- theta[3]
  mu.xi <- theta[4]
  sigma2.eta <- theta[5]
  sigma2.xi <-theta[6]
  if(theta[5] < 0 | theta[6] < 0){ ## check variances
    return(NA)
  }
  lik <- 0.0
  ## nodes and weights for Gauss-Hermite quadrature
  objGH <- gauss.quad(n.node, "hermite")
  w <- objGH$weights
  w <- cbind(rep(w,rep(n.node,n.node)),rep(w,n.node))
  w <- w[,1]*w[,2]
  node <- objGH$nodes
  node <- cbind(rep(node,rep(n.node,n.node)),rep(node,n.node))
  if (model == "approx"){
    for(i in 1:n){
      g <- function(x){dnorm(dati[i,1],x[1],sqrt(dati[i,4]))*dnorm(dati[i,2],x[2],sqrt(dati[i,7]))*dnorm(x[1],beta0+beta1*(x[2]-mu.xi)+beta2*(x[2]-mu.xi)^2,sqrt(sigma2.eta))*dnorm(x[2],mu.xi,sqrt(sigma2.xi))*exp(x[1]^2+x[2]^2)}
      h <- apply(node, MARGIN=1, FUN=g)
      lik <- lik + log(sum(w*h))
    }
  } 
  if (model == "exact"){
    for(i in 1:n){
      g<-function(x){
        p.eta <- exp(x[1])/(1+exp(x[1]))
        p.xi <- exp(x[2])/(1+exp(x[2]))
        return(dbinom(dati[i,8],dati[i,10],p.eta)*dbinom(dati[i,9],dati[i,11],p.xi)*dnorm(x[1],beta0+beta1*(x[2]-mu.xi)+beta2*(x[2]-mu.xi)^2,sqrt(sigma2.eta))*dnorm(x[2],mu.xi,sqrt(sigma2.xi))*exp(x[1]^2+x[2]^2))
      }    
      h <- apply(node, MARGIN=1, FUN=g)
      lik <- lik + log(sum(w*h))
    }
  }
  return(lik)
}

quad.lik.GH2.repa <- function(theta, dati, n.node=10, model="approx"){
  n <- nrow(dati)
  beta0 <- theta[1]
  beta1 <- theta[2]
  beta2 <- theta[3]
  mu.xi <- theta[4]
  log.sigma2.eta <- theta[5]
  log.sigma2.xi <-theta[6]
  lik <- 0.0
  ## nodes and weights for Gauss-Hermite quadrature
  objGH <- gauss.quad(n.node, "hermite")
  w <- objGH$weights
  w <- cbind(rep(w,rep(n.node,n.node)),rep(w,n.node))
  w <- w[,1]*w[,2]
  node <- objGH$nodes
  node <- cbind(rep(node,rep(n.node,n.node)),rep(node,n.node))
  if (model == "approx"){
    for(i in 1:n){
      g <- function(x){dnorm(dati[i,1],x[1],sqrt(dati[i,4]))*dnorm(dati[i,2],x[2],sqrt(dati[i,7]))*dnorm(x[1],beta0+beta1*(x[2]-mu.xi)+beta2*(x[2]-mu.xi)^2,sqrt(exp(log.sigma2.eta)))*dnorm(x[2],mu.xi,sqrt(exp(log.sigma2.xi)))*exp(x[1]^2+x[2]^2)}
      h <- apply(node, MARGIN=1, FUN=g)
      lik <- lik + log(sum(w*h))
    }
  } 
  if (model == "exact"){
    for(i in 1:n){
      g<-function(x){
        p.eta <- exp(x[1])/(1+exp(x[1]))
        p.xi <- exp(x[2])/(1+exp(x[2]))
        return(dbinom(dati[i,8],dati[i,10],p.eta)*dbinom(dati[i,9],dati[i,11],p.xi)*dnorm(x[1],beta0+beta1*(x[2]-mu.xi)+beta2*(x[2]-mu.xi)^2,sqrt(exp(log.sigma2.eta)))*dnorm(x[2],mu.xi,sqrt(exp(log.sigma2.xi)))*exp(x[1]^2+x[2]^2))
      }    
      h <- apply(node, MARGIN=1, FUN=g)
      lik <- lik + log(sum(w*h))
    }
  }
  return(lik)
}

quad.lik.GH3 <- function(theta, dati, n.node=10, model="approx"){
  n <- nrow(dati)
  beta0 <- theta[1]
  beta1 <- theta[2]
  beta2 <- theta[3]
  mu.xi <- theta[4]
  sigma2.eta <- theta[5]
  sigma2.xi <-theta[6]
  if(theta[5] < 0 | theta[6] < 0){ ## check variances
    return(NA)
  }
  lik <- 0.0
  ## nodes and weights for Gauss-Hermite quadrature
  objGH <- gauss.quad(n.node, "hermite")
  w <- objGH$weights
  w <- cbind(rep(w,rep(n.node,n.node)),rep(w,n.node))
  w <- w[,1]*w[,2]
  node <- objGH$nodes
  node <- cbind(rep(node,rep(n.node,n.node)),rep(node,n.node))
  if (model == "approx"){
    for(i in 1:n){
      g <- function(x){
        xi <- sqrt(2*dati[i,7])*x[2]+dati[i,2]
        eta <- sqrt(2*dati[i,4])*x[1]+dati[i,1]
        return(dnorm(eta,beta0+beta1*(xi-mu.xi)+beta2*(xi-mu.xi)^2,sqrt(sigma2.eta))*dnorm(xi,mu.xi,sqrt(sigma2.xi))/sqrt(pi*pi))
      }
      h <- apply(node, MARGIN=1, FUN=g)
      lik <- lik + log(sum(w*h))
    }
  } 
  if (model == "exact"){
    for(i in 1:n){
      g<-function(x){
        xi <- sqrt(2*dati[i,7])*x[2]+dati[i,2]
        eta <- sqrt(2*dati[i,4])*x[1]+dati[i,1]
        p.eta <- exp(eta)/(1+exp(eta))
        p.xi <- exp(xi)/(1+exp(xi))
        return(dbinom(dati[i,8],dati[i,10],p.eta)*dbinom(dati[i,9],dati[i,11],p.xi)*dnorm(eta,beta0+beta1*(xi-mu.xi)+beta2*(xi-mu.xi)^2,sqrt(sigma2.eta))*dnorm(xi,mu.xi,sqrt(sigma2.xi))*exp(x[1]^2+x[2]^2)*sqrt(2*dati[i,4]*2*dati[i,7]))
      }    
      h <- apply(node, MARGIN=1, FUN=g)
      lik <- lik + log(sum(w*h))
    }
  }
  return(lik)
}

quad.lik.GH3.repa <- function(theta, dati, n.node=10, model="approx"){
  n <- nrow(dati)
  beta0 <- theta[1]
  beta1 <- theta[2]
  beta2 <- theta[3]
  mu.xi <- theta[4]
  log.sigma2.eta <- theta[5]
  log.sigma2.xi <-theta[6]
  lik <- 0.0
  ## nodes and weights for Gauss-Hermite quadrature
  objGH <- gauss.quad(n.node, "hermite")
  w <- objGH$weights
  w <- cbind(rep(w,rep(n.node,n.node)),rep(w,n.node))
  w <- w[,1]*w[,2]
  node <- objGH$nodes
  node <- cbind(rep(node,rep(n.node,n.node)),rep(node,n.node))
  if (model == "approx"){
    for(i in 1:n){
      g <- function(x){
        xi <- sqrt(2*dati[i,7])*x[2]+dati[i,2]
        eta <- sqrt(2*dati[i,4])*x[1]+dati[i,1]
        return(dnorm(eta,beta0+beta1*(xi-mu.xi)+beta2*(xi-mu.xi)^2,sqrt(exp(log.sigma2.eta)))*dnorm(xi,mu.xi,sqrt(exp(log.sigma2.xi)))/sqrt(pi*pi))
      }
      h <- apply(node, MARGIN=1, FUN=g)
      lik <- lik + log(sum(w*h))
    }
  } 
  if (model == "exact"){
    for(i in 1:n){
      g<-function(x){
        xi <- sqrt(2*dati[i,7])*x[2]+dati[i,2]
        eta <- sqrt(2*dati[i,4])*x[1]+dati[i,1]
        p.eta <- exp(eta)/(1+exp(eta))
        p.xi <- exp(xi)/(1+exp(xi))
        return(dbinom(dati[i,8],dati[i,10],p.eta)*dbinom(dati[i,9],dati[i,11],p.xi)*dnorm(eta,beta0+beta1*(xi-mu.xi)+beta2*(xi-mu.xi)^2,sqrt(exp(log.sigma2.eta)))*dnorm(xi,mu.xi,sqrt(exp(log.sigma2.xi)))*exp(x[1]^2+x[2]^2)*sqrt(2*dati[i,4]*2*dati[i,7]))
      }    
      h <- apply(node, MARGIN=1, FUN=g)
      lik <- lik + log(sum(w*h))
    }
  }
  return(lik)
}

quad.lik.GH4.repa <- function(theta, dati, n.node=10, model="approx"){
  n <- nrow(dati)
  beta0 <- theta[1]
  beta1 <- theta[2]
  beta2 <- theta[3]
  mu.xi <- theta[4]
  log.sigma2.eta <- theta[5]
  log.sigma2.xi <-theta[6]
  lik <- 0.0
  ## nodes and weights for Gauss-Hermite quadrature
  objGH <- gauss.quad(n.node, "hermite")
  w <- objGH$weights
  w <- cbind(rep(w,rep(n.node,n.node)),rep(w,n.node))
  w <- w[,1]*w[,2]
  node <- objGH$nodes
  node <- cbind(rep(node,rep(n.node,n.node)),rep(node,n.node))
  if (model == "approx"){
    for(i in 1:n){
      g <- function(x){
        xi <- sqrt(2*exp(log.sigma2.xi))*x[2]+mu.xi
        eta <- sqrt(2*exp(log.sigma2.eta))*x[1]+beta0+beta1*(xi-mu.xi)+beta2*(xi-mu.xi)^2
        return(dnorm(dati[i,1],eta,dati[i,4])*dnorm(dati[i,2],xi,dati[i,7]))
      }
      h <- apply(node, MARGIN=1, FUN=g)
      lik <- lik + log(sum(w*h))
    }
  } 
  if (model == "exact"){
    for(i in 1:n){
      g<-function(x){
        xi <- sqrt(2*exp(log.sigma2.xi))*x[2]+mu.xi
        eta <- sqrt(2*exp(log.sigma2.eta))*x[1]+beta0+beta1*(xi-mu.xi)+beta2*(xi-mu.xi)^2
        p.eta <- exp(eta)/(1+exp(eta))
        p.xi <- exp(xi)/(1+exp(xi))
        return(dbinom(dati[i,8],dati[i,10],p.eta)*dbinom(dati[i,9],dati[i,11],p.xi))
      }    
      h <- apply(node, MARGIN=1, FUN=g)
      lik <- lik + log(sum(w*h))
    }
  }
  return(lik)
}

quad.lik.GH5.repa <- function(theta, dati, n.node=10, model="approx"){
  n <- nrow(dati)
  beta0 <- theta[1]
  beta1 <- theta[2]
  beta2 <- theta[3]
  mu.xi <- theta[4]
  log.sigma2.eta <- theta[5]
  log.sigma2.xi <-theta[6]
  lik <- 0.0
  ## nodes and weights for Gauss-Hermite quadrature
  objGH <- gauss.quad(n.node, "hermite")
  w <- objGH$weights
  w <- cbind(rep(w,rep(n.node,n.node)),rep(w,n.node))
  w <- w[,1]*w[,2]
  node <- objGH$nodes
  node <- cbind(rep(node,rep(n.node,n.node)),rep(node,n.node))
  if (model == "approx"){
    for(i in 1:n){
      g <- function(x){
        xi <- sqrt(2*exp(log.sigma2.xi))*x[2]+dati[i,2]
        eta <- sqrt(2*exp(log.sigma2.eta))*x[1]+dati[i,1]
        return(dnorm(dati[i,1],eta,sqrt(dati[i,4]))*dnorm(dati[i,2],xi,sqrt(dati[i,7]))*dnorm(eta,beta0+beta1*(xi-mu.xi)+beta2*(xi-mu.xi)^2,sqrt(exp(log.sigma2.eta)))*dnorm(xi,mu.xi,sqrt(exp(log.sigma2.xi)))*exp(x[1]^2+x[2]^2)*sqrt(2*exp(log.sigma2.eta)*2*exp(log.sigma2.xi)))
      }
      h <- apply(node, MARGIN=1, FUN=g)
      lik <- lik + log(sum(w*h))
    }
  } 
  if (model == "exact"){
    for(i in 1:n){
      g<-function(x){
        xi <- sqrt(2*exp(log.sigma2.xi))*x[2]+dati[i,2]
        eta <- sqrt(2*exp(log.sigma2.eta))*x[1]+dati[i,1]
        p.eta <- exp(eta)/(1+exp(eta))
        p.xi <- exp(xi)/(1+exp(xi))
        return(dbinom(dati[i,8],dati[i,10],p.eta)*dbinom(dati[i,9],dati[i,11],p.xi)*dnorm(eta,beta0+beta1*(xi-mu.xi)+beta2*(xi-mu.xi)^2,sqrt(exp(log.sigma2.eta)))*dnorm(xi,mu.xi,sqrt(exp(log.sigma2.xi)))*exp(x[1]^2+x[2]^2)*sqrt(2*exp(log.sigma2.eta)*2*exp(log.sigma2.xi)))
      }    
      h <- apply(node, MARGIN=1, FUN=g)
      lik <- lik + log(sum(w*h))
    }
  }
  return(lik)
}

quad.lik.GH6.repa <- function(theta, dati, n.node=10, model="approx"){
  n <- nrow(dati)
  beta0 <- theta[1]
  beta1 <- theta[2]
  beta2 <- theta[3]
  mu.xi <- theta[4]
  log.sigma2.eta <- theta[5]
  log.sigma2.xi <-theta[6]
  lik <- 0.0
  ## nodes and weights for Gauss-Hermite quadrature
  objGH <- gauss.quad(n.node, "hermite")
  w <- objGH$weights
  w <- cbind(rep(w,rep(n.node,n.node)),rep(w,n.node))
  w <- w[,1]*w[,2]
  node <- objGH$nodes
  node <- cbind(rep(node,rep(n.node,n.node)),rep(node,n.node))
  if (model == "approx"){
    for(i in 1:n){
      g <- function(x){
        xi <- sqrt(2*exp(log.sigma2.xi))*x[2]+mu.xi
        eta <- sqrt(2*(beta1^2*exp(log.sigma2.xi)+2*beta2^2*exp(log.sigma2.xi)^2+exp(log.sigma2.eta)))*x[1]+beta0+beta2*exp(log.sigma2.xi)
        return(dnorm(dati[i,1],eta,sqrt(dati[i,4]))*dnorm(dati[i,2],xi,sqrt(dati[i,7]))*dnorm(eta,beta0+beta1*(xi-mu.xi)+beta2*(xi-mu.xi)^2,sqrt(exp(log.sigma2.eta)))*exp(x[1]^2)*sqrt(beta1^2*exp(log.sigma2.xi)+2*beta2^2*exp(log.sigma2.xi)^2+exp(log.sigma2.eta)))
      }
      h <- apply(node, MARGIN=1, FUN=g)
      lik <- lik + log(sum(w*h))
    }
  } 
  if (model == "exact"){
    for(i in 1:n){
      g<-function(x){
        xi <- sqrt(2*exp(log.sigma2.xi))*x[2]+mu.xi
        eta <- sqrt(2*(beta1^2*exp(log.sigma2.xi)+2*beta2^2*exp(log.sigma2.xi)^2+exp(log.sigma2.eta)))*x[1]+beta0+beta2*exp(log.sigma2.xi)
        p.eta <- exp(eta)/(1+exp(eta))
        p.xi <- exp(xi)/(1+exp(xi))
        return(dbinom(dati[i,8],dati[i,10],p.eta)*dbinom(dati[i,9],dati[i,11],p.xi)*dnorm(eta,beta0+beta1*(xi-mu.xi)+beta2*(xi-mu.xi)^2,sqrt(exp(log.sigma2.eta)))*exp(x[1]^2)*sqrt(2*(beta1^2*exp(log.sigma2.xi)+2*beta2^2*exp(log.sigma2.xi)^2+exp(log.sigma2.eta))))
      }    
      h <- apply(node, MARGIN=1, FUN=g)
      lik <- lik + log(sum(w*h))
    }
  }
  return(lik)
}

quad.lik.GH7.repa <- function(theta, dati, n.node=10, model="approx"){
  n <- nrow(dati)
  beta0 <- theta[1]
  beta1 <- theta[2]
  beta2 <- theta[3]
  mu.xi <- theta[4]
  log.sigma2.eta <- theta[5]
  log.sigma2.xi <-theta[6]
  lik <- 0.0
  ## nodes and weights for Gauss-Hermite quadrature
  objGH <- gauss.quad(n.node, "hermite")
  w <- objGH$weights
  w <- cbind(rep(w,rep(n.node,n.node)),rep(w,n.node))
  w <- w[,1]*w[,2]
  node <- objGH$nodes
  node <- cbind(rep(node,rep(n.node,n.node)),rep(node,n.node))
  if (model == "approx"){
    for(i in 1:n){
      g <- function(x){
        xi <- sqrt(2*(dati[i,7]))*x[2]+dati[i,2]
        eta <- sqrt(2*exp(log.sigma2.eta))*x[1]+dati[i,1]
        return(dnorm(dati[i,1],eta,sqrt(dati[i,4]))*dnorm(dati[i,2],xi,sqrt(dati[i,7]))*dnorm(eta,beta0+beta1*(xi-mu.xi)+beta2*(xi-mu.xi)^2,sqrt(exp(log.sigma2.eta)))*dnorm(xi,mu.xi,sqrt(exp(log.sigma2.xi)))*exp(x[1]^2+x[2]^2)*sqrt(2*exp(log.sigma2.eta)*2*(dati[i,7])))
      }
      h <- apply(node, MARGIN=1, FUN=g)
      lik <- lik + log(sum(w*h))
    }
  } 
  if (model == "exact"){
    for(i in 1:n){
      g<-function(x){
        xi <- sqrt(2*(dati[i,7]))*x[2]+dati[i,2]
        eta <- sqrt(2*exp(log.sigma2.eta))*x[1]+dati[i,1]
        p.eta <- exp(eta)/(1+exp(eta))
        p.xi <- exp(xi)/(1+exp(xi))
        return(dbinom(dati[i,8],dati[i,10],p.eta)*dbinom(dati[i,9],dati[i,11],p.xi)*dnorm(eta,beta0+beta1*(xi-mu.xi)+beta2*(xi-mu.xi)^2,sqrt(exp(log.sigma2.eta)))*dnorm(xi,mu.xi,sqrt(exp(log.sigma2.xi)))*exp(x[1]^2+x[2]^2)*sqrt(2*exp(log.sigma2.eta)*2*(dati[i,7])))
      }    
      h <- apply(node, MARGIN=1, FUN=g)
      lik <- lik + log(sum(w*h))
    }
  }
  return(lik)
}

quad.lik.GH8.repa <- function(theta, dati, n.node=10, model="approx"){
  n <- nrow(dati)
  beta0 <- theta[1]
  beta1 <- theta[2]
  beta2 <- theta[3]
  mu.xi <- theta[4]
  log.sigma2.eta <- theta[5]
  log.sigma2.xi <-theta[6]
  lik <- 0.0
  ## nodes and weights for Gauss-Hermite quadrature
  objGH <- gauss.quad(n.node, "hermite")
  w <- objGH$weights
  w <- cbind(rep(w,rep(n.node,n.node)),rep(w,n.node))
  w <- w[,1]*w[,2]
  node <- objGH$nodes
  node <- cbind(rep(node,rep(n.node,n.node)),rep(node,n.node))
  if (model == "approx"){
    for(i in 1:n){
      g <- function(x){
        xi <- sqrt(2*exp(log.sigma2.xi))*x[2]+dati[i,2]
        eta <- sqrt(2*exp(log.sigma2.eta))*x[1]+dati[i,1]
        return(dnorm(dati[i,1],eta,sqrt(dati[i,4]))*dnorm(dati[i,2],xi,sqrt(dati[i,7]))*dnorm(eta,beta0+beta1*(xi-mu.xi)+beta2*(xi-mu.xi)^2,sqrt(exp(log.sigma2.eta)))*dnorm(xi,mu.xi,sqrt(exp(log.sigma2.xi)))*exp(x[1]^2+x[2]^2)*sqrt(2*exp(log.sigma2.eta)*2*exp(log.sigma2.xi)))
      }
      h <- apply(node, MARGIN=1, FUN=g)
      lik <- lik + log(sum(w*h))
    }
  } 
  if (model == "exact"){
    for(i in 1:n){
      g<-function(x){
        xi <- sqrt(2*exp(log.sigma2.xi))*x[2]+dati[i,2]
        eta <- sqrt(2*exp(log.sigma2.eta))*x[1]+dati[i,1]
        p.eta <- exp(eta)/(1+exp(eta))
        p.xi <- exp(xi)/(1+exp(xi))
        return(dbinom(dati[i,8],dati[i,10],p.eta)*dbinom(dati[i,9],dati[i,11],p.xi)*dnorm(eta,beta0+beta1*(xi-mu.xi)+beta2*(xi-mu.xi)^2,sqrt(exp(log.sigma2.eta)))*dnorm(xi,mu.xi,sqrt(exp(log.sigma2.xi)))*exp(x[1]^2+x[2]^2)*sqrt(2*exp(log.sigma2.eta)*2*exp(log.sigma2.xi)))
      }    
      h <- apply(node, MARGIN=1, FUN=g)
      lik <- lik + log(sum(w*h))
    }
  }
  return(lik)
}

ee.hat <- function(theta, dati, correct=TRUE){
  n <- nrow(dati)
  beta0 <- theta[1]
  beta1 <- theta[2]
  beta2 <- theta[3]
  tau2 <- theta[4]
  ee.hat <- matrix(NA, nrow=4, ncol=1)
  if(correct==FALSE){
    ee.hat[1,1] <- sum(dati[,1]-beta0-beta1*dati[,2]-beta2*dati[,3])
    ee.hat[2,1] <- sum((dati[,1]-beta0-beta1*dati[,2]-beta2*dati[,3])*dati[,2]+beta1*dati[,8])
    ee.hat[3,1] <- sum((dati[,1]-beta0-beta1*dati[,2]-beta2*dati[,3])*dati[,3]+beta2*dati[,12])
    ee.hat[4,1] <- sum(tau2-(dati[,1]-beta0-beta1*dati[,2]-beta2*dati[,3])^2+dati[,4]+beta1^2*dati[,8]+beta2^2*dati[,12])
  }  
  return(ee.hat)
}

## bias, se, sd and convergence 
# basic <- function(x, sex, value=1){
#     conv <- which(!is.na(sex) & !is.na(x))
#     x <- x[conv]
#     sex <- sex[conv]
#     bias <- mean(x)-value
#     se <- mean(sex)
#     SD <- sd(x)
#     return(c(round(bias,3), round(se,3), round(SD,3), length(x)))
# }
    
evaluate <- function(x, sex, value){ ## more than one parameters
    result<-matrix(NA, ncol=7, nrow=ncol(x))
    if ((ncol(x)==ncol(sex)) & (ncol(x)==length(value))){
        for (i in 1:ncol(x)) {
            conv <- which(!is.na(sex[,i]) & !is.na(x[,i])) ## convergent elements
            result[i,1] <- mean(x[conv,i])-value[i] ## bias
            result[i,2] <- mean(sex[conv,i])
            result[i,3] <- sd(x[conv,i])
            #result[i,4] <- sqrt(mean(sex[conv,i]^2))
            result[i,4] <- length(conv)
            if (value[i]!=0) {
              result[i,5] <- result[i,1]/value[i]
            }
            else {
              result[i,5] <- 1000 
            }
            result[i,6] <- sd(x[conv,i])/mean(x[conv,i])
            result[i,7] <- mean((x[conv,i]-value[i])^2)
        }
        result <- data.frame(result)
        colnames(result) <- c("bias", "se", "sd", "convergence", "rebias", "cv", "mse")
        return(result)
    }else{return(NA)}
    # conv <- which(!is.na(sex) & !is.na(x)) ## index of convergent estimates
    # x <- x[conv] ## extract convergent estimates
    # sex <- sex[conv] ## extract corresponding sd
    # bias <- mean(x)-value
    # se <- mean(sex)
    # SD <- sd(x)
    # return(c(round(bias,3), round(se,3), round(SD,3), length(x)))
}

## empirical coverage probability
# coverage <- function(x, sex, value=1, level=0.95){
#     val <- na.omit((x-value)/sex)
#     ris <- mean(abs(val) <= qnorm((1+level)/2))
#     return(c(ris, length(val)))
# }

emp.coverage <- function(x, sex, value, level=0.95){ ## more than one parameters
    result<-matrix(NA, ncol=2, nrow=ncol(x))
    if ((ncol(x)==ncol(sex)) & (ncol(x)==length(value))){
        for (i in 1:ncol(x)) {
            val <- na.omit((x[,i]-value[i])/sex[,i])
            result[i,1] <- mean(abs(val) <= qnorm((1+level)/2))
            result[i,2] <- length(val)
        }
        result <- data.frame(result)
        colnames(result) <- c("coverage", "convergence")
        return(result)
    }else{return(NA)}
    # val <- na.omit((x-value)/sex) ## compute standardized difference of estimate and true value
    # ris <- mean(abs(val) <= qnorm((1+level)/2)) ## compare the standardized difference with quantiles
    # return(c(ris, length(val)))
}

## sandwich estimate of the asymptotic covariance matrix
#library(nlme)

sand.lik <- function(theta.mle, lik, dati){
  
  #a <- nlme::fdHess(theta.mle, lik, dati=dati[1,])
  #values.gradient <- a$gradient 
  G <- matrix(0,length(theta.mle),length(theta.mle))
  #G <- values.gradient%*%t(values.gradient)
  H <- matrix(0,length(theta.mle),length(theta.mle))
  #H <- a$Hessian         
  for(i in 1:nrow(dati)){
    a <- fdHess(theta.mle, lik, dati=dati[i,]) ## compute Hessian matrix
    values.gradient <- a$gradient 
    G <- G + values.gradient%*%t(values.gradient) ## compute J matrix
    H <- H + a$Hessian   
  }
  return(solve(H)%*%G%*%solve(H))
}

sand.lik2 <- function(theta.mle, lik, dati){
  
  #a <- nlme::fdHess(theta.mle, lik, dati=dati[1,])
  #values.gradient <- a$gradient 
  G <- matrix(0,length(theta.mle),length(theta.mle))
  #G <- values.gradient%*%t(values.gradient)
  H <- matrix(0,length(theta.mle),length(theta.mle))
  #H <- a$Hessian         
  for(i in 1:nrow(dati)){
    a <- fdHess(theta.mle, lik, dati=dati[i,]) ## compute Hessian matrix
    values.gradient <- a$gradient 
    G <- G + values.gradient%*%t(values.gradient) ## compute J matrix
    H <- H + a$Hessian   
  }
  return(solve(H,tol=1e-20)%*%G%*%solve(H,tol=1e-20))
}