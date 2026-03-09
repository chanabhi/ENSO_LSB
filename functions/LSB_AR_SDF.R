source("C:/Users/chana/OneDrive/Desktop/LSB_links.R")
source("C:/Users/chana/OneDrive/Desktop/LSB_ARp.R")

dB <- function(x) 10 * log10(x)

ar.sdf <- function(freqs, phi, sigma2) {
  p <- length(phi)
  sapply(freqs, function(f) {
    z <- exp(-2i * pi * f * (1:p))
    denom <- 1 - sum(phi * z)
    sigma2 / (2 * pi * Mod(denom)^2)
  })
}

my.acf <- function (x, lag, type = "correlation") {
  
  acf(x, plot=FALSE, lag.max=lag, type = type)$acf[lag+1]
}


wstats <- function (x, FUN, B, ...) {
  
  N <- length(x)
  
  bs <- 1:B
  ks <- 0:(N-B-1)
  
  med.B <- median(bs)
  
  structure(list(indexes=ks+med.B,
                 stats=lapply(ks, function (k) FUN(x[k+bs], ...)),
                 x=x, N=N),
            class="wstats")
}


plot.wstats <- function (obj, ts, ...) {
  
  if (missing(ts)) {
    
    plot(obj$indexes, unlist(obj$stats), xlim=range(1:obj$N), type="l", ...)            
    
  } else {
    
    plot(ts[obj$indexes], unlist(obj$stats), xlim=range(ts), type="l", ...)
  }
}


LSB.ARp.sdf <- function(par, X, fs){
  
  require(spectral)
  
  T. <- nrow(X[[1]])
  p <- length(X) - 1
  
  the.beta <- theta.to.beta.list(X, par, p)
  
  the.phi <- beta.list.to.phi.list(X, the.beta, p)
  
  the.sigma.sq <- eta.to.sigma2(calc.eta(X[[p+1]], the.beta[[p+1]]))
  
  ar.sdfs <- 
    lapply(1:(T.-p), 
           function(k) dB(ar.sdf(freqs = fs, 
                                 phi = as.vector(the.phi[[p]][1:p,k]), 
                                 sigma2 = the.sigma.sq[k])))
  
  t(sapply(ar.sdfs, function (x) x))
  
  
}



LSB.AR2.log.sdf <- function(phi.21, phi.22, sigma.sq, fs){
  ## ======================================================================
  ## Calculate log(SDF) from AR2 parameter curves phi21 and phi22
  ## for fixed sigma.sq and frequencies fs
  ## ======================================================================
  
  log.ar.sdfs <- 
    lapply(1:length(phi.21), 
           function(k) log(ar.sdf(freqs = fs, 
                                  phi = c(phi.21[k], phi.22[k]), 
                                  sigma2 = sigma.sq)))
  
  t(sapply(log.ar.sdfs, function (x) x))
  
}



LSB.AR2.phi.log.sdf <- function(X, theta, fs){
  ## ======================================================================
  ## Calculate log(SDF) and AR2 parameter curves phi21 and phi22
  ## from LSB basis parameters 
  ## ======================================================================
  
  n = nrow(X)
  
  beta.tt1 <- theta[1:ncol(X)]
  beta.tt2 <- theta[(ncol(X) + 1):(2 * (ncol(X)))]
  beta.tt3 <- theta[length(theta)]
  
  
  phi.11 <- eta.to.phi(calc.eta(X[-1, ], beta.tt1))
  phi.22 <- eta.to.phi(calc.eta(X[-c(1, 2), ], beta.tt2))
  sigma.sq <- eta.to.sigma2(beta.tt3)
  
  ts <- 3:n
  
  phi.21 <- array(NA, length(phi.22))
  
  phi.21[ts - 2] <- phi.11[ts - 2] - phi.22[ts - 2] * phi.11[ts - 2]
  
  log.ar.sdfs <- 
    lapply(1:(n-2), 
           function(k) log(ar.sdf(freqs = fs, 
                                 phi = c(phi.21[k], phi.22[k]), 
                                 sigma2 = sigma.sq)))
  
  the.log.ar.sdfs <- t(sapply(log.ar.sdfs, function (x) x))
  
  list(phi.21 = phi.21,
       phi.22 = phi.22,
       log.lsb2.sdf = the.log.ar.sdfs)
  
}
