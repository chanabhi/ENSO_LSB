source("C:/Users/chana/OneDrive/Desktop/LSB_links.R")

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
  
  # This loop was missing entirely
  for(i in 1:p){
    the.phi[[i]] <- matrix(NA, i, (N-i))
    the.phi[[i]][i,] <- calc.phi(as.matrix(X[[i]][-seq(1, i, 1),]), beta[[i]])
  }
  
  if (p >= 2) {
    for (i in 2:p) {
      for (j in 1:(i - 1)) {
        the.phi[[i]][j, ] <- the.phi[[i - 1]][j, -1] -
          the.phi[[i]][i, ] * the.phi[[i - 1]][i - j, -1]
      }
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
  
  phi_denom <- as.numeric(sapply(1:p, function(j) (1 - phi[[j]][j,1]^2)))
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
  
  phi_denom <- as.numeric(sapply(1:p, function(j) (1 - phis[[j]][j,1]^2)))
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
              y=enso, X = X, method = "L-BFGS-B", hessian = FALSE)
  
  list(
    ar.order = p,
    ar.basis.order = b,
    var.basis.order = s,
    NIC = fit$value / T. + length(fit$par) / T.,
    BIC = fit$value / T. + (length(fit$par) * log(T.)) / (2 * T.),
    par = fit$par
  )
}



