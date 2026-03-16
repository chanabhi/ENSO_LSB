# Libraries & sources
library(astsa)
library(spectral)
library(fields)
library(readxl)
library(mvnfast)

source("functions/basis_functions.R")
source("functions/LSB_ARp.R")
source("functions/LSB_AR_SDF.R")

# Load & prepare data
raw <- read_excel("data/ENSO.xlsx")

enso <- as.vector(t(as.matrix(raw[, -1])))
enso <- enso[!is.na(enso)]

T.       <- length(enso)
ts_years <- seq(1951, by = 1/12, length.out = T.)
fs       <- seq(0, 0.5, length = 128)

# Plot raw series
plot(ts_years, enso, type = "l",
     xlab = "Year", ylab = "ENSO index",
     main = "ENSO monthly index")

# Model selection (AR order p, basis order b)
param_grid <- expand.grid(ar.order = 2:6, ar.basis.order = 2:6)

models <- mapply(
  function(p, b) LSB.ARp.model.sel(enso, p, b, X.method = "spline"),
  param_grid$ar.order,
  param_grid$ar.basis.order,
  SIMPLIFY = FALSE
)

best_nic <- models[[which.min(sapply(models, `[[`, "NIC"))]]
best_bic <- models[[which.min(sapply(models, `[[`, "BIC"))]]

cat("NIC-selected model: AR order =", best_nic$ar.order,
    " | basis order =", best_nic$ar.basis.order, "\n")
cat("BIC-selected model: AR order =", best_bic$ar.order,
    " | basis order =", best_bic$ar.basis.order, "\n")

# Select variance basis order (s) via NIC 
b_values   <- 2:6
var_models <- lapply(b_values, function(s)
  LSB.ARp.model.sel(enso, best_nic$ar.order, best_nic$ar.basis.order,
                    s = s, X.method = "spline"))

best_s <- b_values[which.min(sapply(var_models, `[[`, "NIC"))]
cat("Variance basis order:", best_s, "\n")

# Fit time-varying LSB-AR model
p <- best_nic$ar.order
b <- best_nic$ar.basis.order
s <- best_s

build_basis <- function(orders) {
  lapply(orders, function(d)
    build.X(n = T., X.deg = d, X.method = "spline"))
}

fit_lsb <- function(basis_orders, ...) {
  X <- build_basis(basis_orders)
  fit <- optim(
    rep(0.01, sum(basis_orders + 1)),
    LSB.ARp.lik,
    y = enso, X = X,
    method = "L-BFGS-B", ...
  )
  list(fit = fit, X = X)
}

tv  <- fit_lsb(c(rep(b, p), s), hessian = TRUE)   # time-varying
sta <- fit_lsb(rep(0, p + 1))                      # stationary (null)

# Extract parameters & compute time-varying SDF
the.beta   <- theta.to.beta.list(tv$X, tv$fit$par, p)
the.phi    <- beta.list.to.phi.list(tv$X, the.beta, p)
the.sigma2 <- eta.to.sigma2(calc.eta(tv$X[[p + 1]], the.beta[[p + 1]]))

ar.sdfs <- lapply(seq_len(T. - p), function(k)
  dB(ar.sdf(freqs  = fs,
            phi    = as.vector(the.phi[[p]][1:p, k]),
            sigma2 = the.sigma2[k])))

sp.matrix <- t(sapply(ar.sdfs, identity))

# Plot time-varying SDF


image.plot(
  ts_years[(p + 1):T.],
  fs * 12,
  sp.matrix,
  xlab = "Year",
  ylab = "Frequency (cycles/year)",
  main = paste0("Time-varying SDF: LSB-AR(", p, ")"),
  ylim = c(0, 3),
  
)
abline(h   = c(1/7, 1/2, 1), lwd = 2, lty = 2)
text(x      = 1960,
     y      = c(0.25, 0.75, 1.1),
     labels = c("ENSO low", "ENSO high", "Annual"),
     cex    = 0.75)

# SD of time-varying SDF 
cov.mat <- solve(tv$fit$hessian)

par.sim <- mvnfast::rmvn(
  1000,
  mu    = tv$fit$par,
  sigma = cov.mat
)

sdfs <- lapply(1:nrow(par.sim), function(i)
  LSB.ARp.sdf(par = par.sim[i, ], X = tv$X, fs = fs))

spec.sd <- array(
  sapply(1:(nrow(sdfs[[1]]) * ncol(sdfs[[1]])), function(x)
    sd(unlist(lapply(sdfs, "[[", x)))),
  dim(sdfs[[1]])
)

# Plot SD of time-varying SDF
image.plot(
  ts_years[(p + 1):T.],
  fs * 12,
  spec.sd,
  xlab = "Year",
  ylab = "Frequency (cycles/year)",
  main = paste0("SD of Time-varying SDF: LSB-AR(", p, ")"),
  ylim = c(0, 3)
)
abline(h   = c(1/7, 1/2, 1), lwd = 2, lty = 2)
text(x      = 1960,
     y      = c(0.25, 0.75, 1.1),
     labels = c("ENSO low", "ENSO high", "Annual"),
     cex    = 0.75)

# Stationarity likelihood-ratio test 
test_stat <- -2 * (tv$fit$value - sta$fit$value)
df        <- sum(c(rep(b, p), s)) - (p + 1)
p_value   <- pchisq(test_stat, df = df, lower.tail = FALSE)

cat(sprintf("Stationarity test: statistic = %.2f, df = %d, p-value = %.4f\n",
            test_stat, df, p_value))