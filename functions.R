#Functions

colorout <- function(x){
  x <- abs(x)
  out <- "0-2"
  if(x > 2){out <- "2-3"}
  if(x > 3){out <- "3+"}
  return(out)
}

ks <- c(100,500*(1:40))

k <- ks


conditionalmeansig <- function(i, Sigma, mu,Beta){
  Sigma11 <- Sigma[i,i]
  Sigma12 <- Sigma[i,]
  Sigma21 <- Sigma[-i,i]
  Sigma22 <- Sigma[-i,-i]
  Sigma12 <- Sigma[i,-i]
  
  condmu <- mu[i]+Sigma12%*%solve(Sigma22)%*%(Beta[-i,]-as.matrix(mu[-i]))
  
  condsig <- Sigma11 - Sigma12%*%solve(Sigma22)%*%Sigma21
  return(list(condmu, condsig))
}




generateinfo <- function(Beta){
  
  y <- NA
  
  X <- as.matrix(cbind(rep(1,100),rnorm(100), rnorm(100)))
  
  rnk <- rankMatrix(X)[1]
  
  if(rnk < dim(Beta)[1]){
    cat("Rank Matrix Problem")
  }
  
  ny <- 100
  
  z <- NA
  
  for(j in 1:ny){
    z[j] <- rt(1, 7, (X%*%Beta)[j])
    y[j] <- 1
    if(z[j] < 0){
      y[j] <- 0
    }
  }
  
  return(list(y,X,Beta))
  
}

trace <- function(x,seed, k = 20000){
  
  x <- data.frame(x[[seed]],rep(0:k,3))
  
  colnames(x) <- c(colnames(x)[-length(colnames(x))],"Iteration")
  
  x <- gather(x, key = Iteration, value = Beta)
  
  colnames(x) <- c("Iteration", "Beta", "Values")
  
  ggplot(x, aes(x = Iteration, y = Values)) + geom_line() + facet_grid(Beta ~ .)
  
}



Gewekedataframe <- function(x, frac1, frac2, beta, Method, chains){
  z <- 0
  Geweke <- data.frame(0,0,0,0,0,0)
  for(i in 1:3){
    for(j in 1:3){
      for(k in 1:3){
        for(l in 1:4){
          for(m in 1:5){
            z <- z + 1
            Geweke[z,] <- c(x[i,j,k,l,m], frac1[i], frac2[j], beta[k], Method[l], chains[m])
          }
        }
      }
    }
    
  }
  colnames(Geweke) <- c("Values", "Init.Fraction", "End.Fraction", "Beta", "Method", "Chain")
  
  return(Geweke)
}


geweke <- function(x, frac1 = frac1, frac2 = frac2){
  
  out <- array(-999, dim = c(3,3,3))
  
  for(i in 1:3){
    for(j in 1:3){
      out[i,j,] <- geweke.diag(x, frac1[i], frac2[j])$z
    }
  }
  return(out)  
}

gewekewrapper <- function(x, frac1, frac2){
  out <- array(-999, dim = c(3,3,3,4,5))
  for(i in 1:4){
    for(j in 1:5){
      out[,,,i,j] <- geweke(x[[i]][[j]],frac1, frac2)
    }
  }
  return(out)
}



BetaAggregate <- function(BetaChains, name, p = 3){
  
  X <- sapply(BetaChains, ACFWrapper)
  
  X.2 <- data.frame(X, rep(paste("Beta", 1:p), each = length(X[,1])/p))
  
  colnames(X.2) <- c(paste("Chain", 1:length(BetaChains), sep = "."), "Beta")
  
  X.2 <- X.2 %>% gather(key = Beta, value = values)
  
  colnames(X.2)[2] <- "Chain"
  
  X.2 <- data.frame(X.2, paste(name))
  
  colnames(X.2) <- c(colnames(X.2)[-length(colnames(X.2))], "Method")
  
  X.2 <- data.frame(X.2, 1:(dim(X.2)/(length(BetaChains)*p))[1])
  
  colnames(X.2) <- c(colnames(X.2)[-length(colnames(X.2))], "Lag")
  
  return(X.2)
  
}





ACF.Gather <- function(ChainList){
  Chains <- sapply(ChainList, ACFWrapper)
  
  Chains.2 <- data.frame(Chains, rep(paste("Beta", 1:p), each = length(Suff.ACF[,1])/p))
  
  colnames(Chains.2) <- c(paste("Chain", 1:length(ChainList), sep = "."), "Beta")
  
  return(Chains.2)
}



reduce <- function(x,k){
  x <- x[1:k,]
  return(x)
}


mean.func <- function(x){
  out <- apply(x, 2, mean)
}

sds <- function(x){
  out <- apply(x, 2, sd)
  return(out)
}



Gelman <- function(data, M, k){
  N <- k
  
  data <- lapply(data, reduce, k = k)
  
  Thetas <- sapply(data, mean.func)
  
  Thetahat <- apply(Thetas,1,mean)
  
  Thetas-Thetahat
  
  Thetas[1,]-Thetahat[1]
  
  B <-  N/(M-1)*apply((Thetas-Thetahat)^2, 1, sum)
  
  sdsvals <- sapply(data,sds)
  
  W <- apply(sdsvals,1,mean)
  
  V <- (N-1)/N*W + (M+1)/(M*N)*B
  
  scale <- (V/W)^.5
  
  
  
  return(scale)
}


list.to.dataframe <- function(x){
  c <- length(x)
  cbind(list[[]])
}

ACF.Inner <- function(x, plot){
  out <- acf(x, plot = plot)$acf
  return(out)
}

ACFWrapper <- function(x){
  out <- apply(x,2,ACF.Inner, plot = FALSE)
  return(out)
}


GelmanWrapper <- function(data, M, ks){
  out <- array(-999, dim = c(3,length(ks)))
  for(i in 1:length(ks)){
    out[,i] <- Gelman(data, M, ks[i])
  }
  
  colnames(out) <- paste(ks)
  return(out)
}





















#####
#x <- data.frame(Sufficient.BetaChains[[1]])

#x <- data.frame(x, 0:20000)

#colnames(x) <- c(colnames(x)[-length(colnames(x))],"Iteration")

#ggplot(x, aes(x = Iteration, y = Beta.3)) + geom_line()



