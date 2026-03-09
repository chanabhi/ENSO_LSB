
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
