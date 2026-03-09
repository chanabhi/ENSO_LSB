

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
