logit <- function (x) {
  
  log(x/(1-x))
}


inv.logit <- function (x) {
  
  exp(x) / (1 + exp(x))
}


delta.to.eta <- function (delta) {
  
  x <- delta + 0.5
  
  log(x/(1-x))
}


eta.to.delta <- function (eta) {
  
  y <- exp(eta) / (1 + exp(eta))
  
  y-0.5
}


phi.to.eta <- function (phi) {
  
  x <- 0.5 * (phi + 1.0)
  
  log(x/(1-x))
}


eta.to.phi <- function (eta) {
  
  y <- exp(eta) / (1.0 + exp(eta))
  
  2.0 * y - 1.0
}



eta.to.sigma2 <- function (eta) {
  
  exp(eta)
}



sigma2.to.eta <- function (sigma2) {
  
  log(sigma2)
}



calc.eta <- function (X, beta) {
  
  drop(X %*% beta)
}


calc.phi <- function (X, beta) {
  
  eta.to.phi(drop(X %*% beta))
}


calc.delta <-function (X, beta) {
  
  eta.to.delta(drop(X %*% beta))
}

poly.basis <- function (N, deg=1, orthogonal=FALSE, normalized=FALSE) {
  
  if (deg==0) {
    
    cbind(rep(1,N))
    
  } else {
    
    ts <- 0:(N-1)
    us <- ts/N    
    if(normalized) {
      
      X <- cbind(1, poly(us, deg, raw = !orthogonal))
      
      t(t(X) / sqrt(apply(X, 2, function (x) sum(x ^ 2))))
      
    } else {
      
      cbind(1, poly(us, deg, raw = !orthogonal))
      
    }
    
  }
}

poly.basis.for.sim <- function(N) {
  
  ts <- 0:(N - 1)
  us <- ts / N
  cbind(1, us, us ^ 2 - (1 / 3))
  
}

fourier.basis <- function(N, deg=1) {
  
  if (deg==0) {
    
    cbind(rep(1,N))
    
  } else {
    
    ts <- 0:(N-1)
    us <- ts/N
    nbasis <- 2 * deg + 1
    d <- seq(1, deg, 1)
    args <- outer(us, d)
    basismat <- matrix(0, N, nbasis)
    basismat[, 1] <- 1
    basismat[, 2 * d] <- sin(args)
    basismat[, (2 * d) + 1] <- cos(args)
    return(basismat)
    
  }
}

spline.basis <- function (N, deg=1, normalized=TRUE) {
  
  if (deg==0) {
    
    cbind(rep(1,N))
    
  } else {
    
    ts <- 0:(N-1)
    us <- ts/N
    
    if (normalized) {
      
      X <- cbind(1, matrix(splines::ns(us, deg), N))
      
      t(t(X) / sqrt(apply(X, 2, function (x) sum(x^2))))
      
    } else {
      
      cbind(1, matrix(splines::ns(us, deg), N))
      
    }
  }
}

build.X <- function(n, X.method, X.deg, orthogonal=TRUE, normalized=TRUE){
  
  if (X.method == "poly") {
    
    poly.basis(N = n, deg = X.deg, orthogonal, normalized)
    
  } else if (X.method == "fourier") {
    
    fourier.basis(N = n, deg = X.deg)
    
  } else if (X.method == "spline") {
    
    spline.basis(N = n, deg = X.deg, normalized)
    
  } else {
    
    stop("X.method has to be one of poly, fourier or spline")
    
  }
}

theta.to.beta.list <- function(X, theta, ord) {
  
  p <- ord
  
  beta <- vector("list", length = p+1)
  beta[[1]] <- theta[1:ncol(as.matrix(X[[1]]))]
  
  beta <- c(list(beta[[1]]), lapply(2:(p + 1), function(i) {
    theta[(sum(sapply(1:(i - 1), function(j)
      ncol(as.matrix(
        X[[j]]
      )))) + 1):(sum(sapply(1:i, function(k)
        ncol(as.matrix(
          X[[k]]
        )))))]
    
  }))
  
  beta
  
}

beta.list.to.phi.list <- function(X, beta, ord){
  
  N <- nrow(X[[1]])
  p <- ord
  
  the.phi <- vector("list", length = p)
  for(i in 1:p){
    the.phi[[i]] <- matrix(NA, i, (N-i))
    the.phi[[i]][i,] <- calc.phi(as.matrix(X[[i]][-seq(1, i, 1),]), beta[[i]])
  }
  
  for (i in 2:p) {
    for (j in 1:(i - 1)) {
      the.phi[[i]][j, ] <- the.phi[[i - 1]][j, -1] -
        the.phi[[i]][i, ] * the.phi[[i - 1]][i - j, -1]
    }
  }
  
  the.phi
  
}


LSB.ARp.sim <- function(X, theta, ord = length(X)-1){
  ## ======================================================================
  ## Simulate an LSB-AR process of order p given a list of 
  ## basis function matrices X
  ## ======================================================================
  
  T. <- nrow(X[[1]])
  p <- ord 
  
  beta <- theta.to.beta.list(X, theta, ord)
  
  phi <- beta.list.to.phi.list(X, beta, ord)
  
  sigma.sq <- eta.to.sigma2(calc.eta(X[[p + 1]], beta[[p + 1]]))
  
  phi_denom = V = array(NA, p)
  
  phi_denom <- sapply(1:p, function(j) (1 - phi[[j]][j, 1] ^ 2))
  V <- sapply(1:p, function(j) sigma.sq[j] / prod(phi_denom[j:p]))
  
  ts1 = 2:p
  ts = (p + 1):T.
  
  y <- array(NA, T.)
  
  y[1] = rnorm(1, 0, sqrt(V[1]))
  
  for (t in ts1) {
    
    y[t] <-
      rnorm(1, sum(phi[[t - 1]][1:(t - 1), 1] * y[(t - 1):1]), sqrt(V[t]))
    
  }
  
  for (t in ts) {
    
    y[t] <-
      rnorm(1, sum(phi[[p]][1:p, (t - p)] * y[(t - 1):(t - p)]), 
            sqrt(sigma.sq[t]))
    
  }
  
  as.vector(y)
  
}


LSB.ARp.lik <- function(y, X, theta, ord = length(X)-1){
  ## ======================================================================
  ## Calculate the exact conditional likelihood for LSB-AR(p) process 
  ## ======================================================================
  
  T. <- length(y)
  k <- length(theta)
  p <- ord 
  
  beta <- theta.to.beta.list(X, theta, ord)
  
  phis <- beta.list.to.phi.list(X, beta, ord)
  
  sigma.sq <- eta.to.sigma2(calc.eta(X[[p+1]], beta[[p+1]]))
  phi_denom = V = array(NA, p)
  
  phi_denom <- sapply(1:p, function(j) (1 - phis[[j]][j,1]^2))
  V <- sapply(1:p, function(j) sigma.sq[j]/prod(phi_denom[j:p]))
  
  
  ts1 = 2:p
  ts = (p+1):T.
  
  - dnorm(y[1], 0, sqrt(V[1]), log = T)  - 
    sum(sapply(ts1, function(ts1){
      sum(dnorm(y[ts1], 
                sum(phis[[ts1-1]][1:(ts1-1), 1]*y[(ts1-1):1]), 
                sqrt(V[ts1]), log = T))
    })) - 
    sum(sapply(ts, function(ts){
      sum(dnorm(y[ts], sum(phis[[p]][1:p, (ts-p)]*y[(ts-1):(ts-p)]), 
                sqrt(sigma.sq[ts]), log = T))
    }))
  
}

LSB.ARp.model.sel <- function(y, p, b, s = NULL, X.method = "spline"){
  ## ======================================================================
  ## Fit basis parameters to time series for given AR order p
  ## and basis order b for time varying parameter curves and 
  ## the variance curve and get NIC and BIC 
  ## ======================================================================
  
  if (is.null(s)) {
    s = b
  }
  beta_ord <- c(rep(b, p), s)
  
  T. <- length(y)
  
  X <- lapply(1:length(beta_ord), function(i)
    build.X(
      n = T.,
      X.deg = beta_ord[i],
      X.method = X.method
    ))  
  
  fit = optim(rep(0.01, sum(beta_ord+1)), LSB.ARp.lik,
              y= y, X = X, method = "L-BFGS-B", hessian = FALSE)
  
  list(
    ar.order = p,
    ar.basis.order = b,
    var.basis.order = s,
    NIC = fit$value / T. + length(fit$par) / T.,
    BIC = fit$value / T. + (length(fit$par) * log(T.)) / (2 * T.),
    par = fit$par
  )
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

dB <- function(x) 10 * log10(x)

ar.sdf <- function(freqs, phi, sigma2) {
  p <- length(phi)
  sapply(freqs, function(f) {
    z <- exp(-2i * pi * f * (1:p))
    denom <- 1 - sum(phi * z)
    sigma2 / (2 * pi * Mod(denom)^2)
  })
}

