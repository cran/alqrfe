

#' @title quantile regression
#'
#' @description Estimate quantile regression with fixed effects for one tau
#'
#' @param x Numeric matrix, covariates
#' @param y Numeric vector, outcome.
#' @param subj Numeric vector, identifies the unit to which the observation belongs.
#' @param tau Numeric, identifies the percentile.
#' @param method Factor, "qr" quantile regression, "qrfe" quantile regression with fixed effects, "lqrfe" Lasso quantile regression with fixed effects, "alqrfe" adaptive Lasso quantile regression with fixed effects.
#' @param ngrid Numeric scalar greater than one, number of BIC to test.
#' @param inf Numeric scalar, internal value, small value.
#' @param digt Numeric scalar, internal value greater than one, define "zero" coefficient.
#'
#' @return alpha       Numeric vector, intercepts' coefficients.
#' @return beta        Numeric vector, exploratory variables' coefficients.
#' @return lambda      Numeric, estimated lambda.
#' @return res         Numeric vector, percentile residuals.
#' @return tau         Numeric scalar, the percentile.
#' @return penalty     Numeric scalar, indicate the chosen effect.
#' @return sig2_alpha  Numeric vector, intercepts' standard errors.
#' @return sig2_beta   Numeric vector, exploratory variables' standard errors.
#' @return Tab_alpha   Data.frame, intercepts' summary.
#' @return Tab_beta    Data.frame, exploratory variables' summary.
#' @return Mat_alpha   Numeric matrix, intercepts' summary.
#' @return Mat_beta    Numeric matrix, exploratory variables' summary.
#' @return method      Factor, method applied.
#'
#' @import Rcpp RcppArmadillo MASS stats
#'
#' @examples
#' # Example 1
#' n = 10
#' m = 5
#' d = 4
#' N = n*m
#' L = N*d
#' x = matrix(rnorm(L), ncol=d, nrow=N)
#' subj = rep(1:n, each=m)
#' alpha = rnorm(n)
#' beta = rnorm(d)
#' eps = rnorm(N)
#' y = x %*% beta  + matrix(rep(alpha, each=m) + eps)
#' y = as.vector(y)
#' m1 = qr(x,y,subj,tau=0.75, method="qrfe")
#' m1
#' m2 = qr(x,y,subj,tau=0.3, method="lqrfe", ngrid = 10)
#' m2
#'
#' # Example 2, from MASS package
#' Rabbit = MASS::Rabbit
#' Rabbit$Treatment = ifelse(Rabbit$Treatment=="Control",0,1)
#' Rabbit$Animal = ifelse(Rabbit$Animal == "R1",1,ifelse(Rabbit$Animal == "R2",2,
#' ifelse(Rabbit$Animal == "R3",3,ifelse(Rabbit$Animal == "R4",4,5))))
#' X = matrix(cbind(Rabbit$Dose,Rabbit$Treatment), ncol=2)
#' m3 = qr(x=X, y=Rabbit$BPchange, subj=Rabbit$Animal,tau=0.5, method="alqrfe", ngrid = 10)
#' m3
#'
#' @references
#' Koenker, R. (2004) "Quantile regression for longitudinal data", J. Multivar. Anal., 91(1): 74-89, <doi:10.1016/j.jmva.2004.05.006>
#'
#' @export
qr = function(x, y, subj, tau = 0.5, method = "qr", ngrid = 20, inf = 1e-8, digt = 4) {
  methods = c("qr", "qrfe", "lqrfe", "alqrfe")
  if (!method %in% methods) stop("Invalid method! Choose from: 'qr', 'qrfe', 'lqrfe', 'alqrfe'")

  inf2 = 1 / 10^digt
  d  = ifelse(is.null(dim(x)[2]), 1, dim(x)[2])
  if (d == 1) x = as.matrix(x)

  dclean = clean_data(y, x, id = subj)
  subj  = as.vector(dclean$id)
  y     = as.vector(dclean$y)
  N     = length(y)
  x     = as.matrix(dclean$x)

  beta  = as.vector(MASS::ginv(t(x) %*% x) %*% t(x) %*% y)
  n     = max(subj)
  alpha = rnorm(n)
  z     = make_z(n, N, id = subj)

  opt = switch(method,
               "qr"     = optim_qr(beta, x, y, tau, N, d),
               "qrfe"   = optim_qrfe(beta, alpha, x, y, z, tau, N, d, n),
               "lqrfe"  = optim_lqr(beta, alpha, x, y, z, tau, N, d, n, ngrid, inf2),
               "alqrfe" = {
                 opt1   = optim_qrfe(beta, alpha, x, y, z, tau, N, d, n)
                 walpha = ifelse(abs(opt1$alpha) <= inf2, inf2, abs(opt1$alpha))
                 wbeta  = ifelse(abs(opt1$beta)  <= inf2, inf2, abs(opt1$beta))
                 optim_alqr(opt1$beta, opt1$alpha, wbeta, walpha, x, y, z, tau, N, d, n, ngrid, inf2)
               }
  )

  alpha  = ifelse(abs(opt$alpha) <= inf2, 0, opt$alpha)
  beta   = ifelse(abs(opt$beta)  <= inf2, 0, opt$beta)
  res    = opt$res
  lambda = opt$lambda

  Sigma = q_cov(alpha, beta, d, inf, n, N, res, method, tau, X = x, Z = z)
  sig2_alpha = Sigma$sig2_alpha
  sig2_beta  = Sigma$sig2_beta

  tab_alpha = f_tab(N, n, d, alpha, sig2_alpha, 1, inf, digt)
  tab_beta  = f_tab(N, n, d, beta,  sig2_beta,  2, inf, digt)

  obj = list(
    alpha = alpha, beta = beta, lambda = lambda, res = res, tau = tau,
    sig2_alpha = sig2_alpha, sig2_beta = sig2_beta,
    Tab_alpha = tab_alpha$Core, Tab_beta = tab_beta$Core,
    Mat_alpha = tab_alpha$Matx, Mat_beta = tab_beta$Matx,
    method = method
  )
  class(obj) = "ALQRFE"
  return(obj)
}


#' @title Print an ALQRFE
#'
#' @description Define the visible part of the object class ALQRFE
#'
#' @param x An object of class "ALQRFE"
#' @param ... further arguments passed to or from other methods.
#'
#' @keywords internal
#' @noRd
print.ALQRFE = function(x,...){
    base::print(x$Tab_beta)
    invisible(x)
}

#' @title optim quantile regression
#'
#' @description This function solves a quantile regression
#'
#' @param beta Numeric vector, initials values.
#' @param x Numeric matrix, covariates.
#' @param y Numeric vector, output.
#' @param tau Numeric scalar, the percentile.
#' @param N Numeric integer, sample size.
#' @param d Numeric integer, X number of columns.
#'
#' @import stats
#' @keywords internal
#' @noRd
optim_qr = function(beta,x,y,tau,N,d){
  Opt    = optim(par=beta, fn = loss_qr, method = "L-BFGS-B",
                 N=N, y=y, x=x, tau=tau, d=d)
  if (Opt$convergence != 0) warning("Optimization did not converge.")
  beta   = Opt$par
  if(d==1) res = y -  (x * beta)
  else     res = y -  (x %*% beta)
  return(list(alpha=0, beta=beta, lambda=0, res=res))
}

#' @title optim quantile regression with fixed effects
#'
#' @description This function solves a quantile regression with fixed effects
#'
#' @param beta Numeric vector, initials values
#' @param alpha Numeric vector, initials values
#' @param x Numeric matrix, covariates.
#' @param y Numeric vector, output.
#' @param z Numeric matrix, incidents.
#' @param tau Numeric scalar, the percentile.
#' @param N Numeric integer, sample size.
#' @param d Numeric integer, X number of columns.
#' @param n Numeric integer, length of alpha.
#'
#' @import stats
#' @keywords internal
#' @noRd
optim_qrfe = function(beta, alpha, x, y, z, tau, N, d, n) {
  theta = c(beta, alpha[-n])  # exclude last alpha for sum-to-zero

  Opt = optim(par = theta, fn = loss_qrfe, method = "L-BFGS-B",
                     x = x, y = y, z = z, tau = tau, n = N, d = d, mm = n)

  if (Opt$convergence != 0) warning("Optimization did not converge.")

  beta = Opt$par[1:d]
  alpha_head = Opt$par[(d + 1):(d + n - 1)]
  alpha = c(alpha_head, -sum(alpha_head))

  x <- as.matrix(x)
  res = y - (z %*% alpha) - (x %*% beta)
  res = as.vector(res)

  return(list(alpha = alpha, beta = beta, lambda = 0, res = res))
}



#' @title optim lasso quantile regression with fixed effects
#'
#' @description This function solves a lasso quantile regression with fixed effects
#'
#' @param beta Numeric vector, initials values
#' @param alpha Numeric vector, initials values
#' @param x Numeric matrix, covariates.
#' @param y Numeric vector, output.
#' @param z Numeric matrix, incidents.
#' @param tau Numeric scalar, the percentile.
#' @param N Numeric integer, sample size.
#' @param d Numeric integer, X number of columns.
#' @param n Numeric integer, length of alpha.
#' @param ngrid Numeric integer, number of iteractions of BIC.
#' @param inf Numeric, internal small quantity.
#'
#' @import stats
#' @keywords internal
#' @noRd
optim_lqr = function(beta, alpha, x, y, z, tau, N, d, n, ngrid, inf) {
  theta = c(beta, alpha[-n])
  p = length(theta)

  lambda_min = sqrt(1 / N)
  lambda_max = max(c(tau, 1 - tau)) * (N / n)

  lambda_tex = exp(seq(log(lambda_min), log(lambda_max), length.out = ngrid))

  bic_tex = numeric(ngrid)
  theta_tex = matrix(0, nrow = p, ncol = ngrid)

  x <- as.matrix(x)

  for (t in seq_len(ngrid)) {
    Opt = optim(par = theta, fn = loss_lqr, method = "L-BFGS-B",
                       n = N, y = y, x = x, z = z, tau = tau, d = d, mm = n, lambda = lambda_tex[t])
    theta_tex[, t] = Opt$par

    beta = theta_tex[1:d, t]
    alpha_head = theta_tex[(d + 1):(d + n - 1), t]
    alpha = c(alpha_head, -sum(alpha_head))

    res = y - (z %*% alpha) - (x %*% beta)

    bic_tex[t] = bic_hat(res, theta = theta_tex[, t], tau, N, p, inf)
  }

  maxw = which.max(bic_tex)
  lambda = lambda_tex[maxw]
  theta = theta_tex[, maxw]

  phi = ifelse(abs(theta) > inf, 1, 0)  # active set indicator

  beta = theta[1:d]
  alpha_head = theta[(d + 1):(d + n - 1)]
  alpha = c(alpha_head, -sum(alpha_head))

  res = y - (z %*% alpha) - (x %*% beta)

  return(list(alpha = alpha, beta = beta, lambda = lambda, res = res, phi = phi))
}


#' @title optim adaptive lasso quantile regression with fixed effects
#'
#' @description This function solves an adaptive lasso quantile regression with fixed effects
#'
#' @param beta Numeric vector, initials values
#' @param alpha Numeric vector, initials values
#' @param wbeta Numeric vector, beta weigths
#' @param walpha Numeric vector, alpha weigths
#' @param x Numeric matrix, covariates.
#' @param y Numeric vector, output.
#' @param z Numeric matrix, incidents.
#' @param tau Numeric scalar, the percentile.
#' @param N Numeric integer, sample size.
#' @param d Numeric integer, X number of columns.
#' @param n Numeric integer, length of alpha.
#' @param ngrid Numeric integer, number of iteractions of BIC.
#' @param inf Numeric, internal small quantity.
#'
#' @import stats
#' @keywords internal
#' @noRd
optim_alqr = function(beta, alpha, wbeta, walpha, x, y, z, tau, N, d, n, ngrid, inf) {
  theta = c(beta, alpha[-n])
  w = c(wbeta, walpha)
  w = ifelse(abs(w) < 1e-10, 1e-10, w)  # Avoid division by zero in penalty

  p = length(theta)
  lambda_min = sqrt(1 / N)
  lambda_max1 = max(c(tau, 1 - tau)) * (N / n)
  lambda_max2 = max(abs(2 * t(x) %*% y / N))
  lambda_max = min(c(lambda_max1, lambda_max2))

  # Log-spaced lambda grid
  lambda_tex = exp(seq(log(lambda_min), log(lambda_max), length.out = ngrid))

  bic_tex = numeric(ngrid)
  theta_tex = matrix(0, nrow = p, ncol = ngrid)

  x <- as.matrix(x)

  for (t in seq_len(ngrid)) {
    Opt = optim(par = theta, fn = loss_alqr, method = "L-BFGS-B",
                       n = N, y = y, x = x, z = z, tau = tau, d = d, mm = n,
                       lambda = lambda_tex[t], w = w)
    theta_tex[, t] = Opt$par

    beta_t = theta_tex[1:d, t]
    alpha_head = theta_tex[(d + 1):(d + n - 1), t]
    alpha_t = c(alpha_head, -sum(alpha_head))

    res = y - (z %*% alpha_t) - (x %*% beta_t)

    bic_tex[t] = bic_hat(res, theta = theta_tex[, t], tau, N, p, inf)
  }

  maxw = which.max(bic_tex)
  lambda = lambda_tex[maxw]
  theta = theta_tex[, maxw]

  # Active set: 1 if parameter magnitude > inf threshold, else 0
  phi = ifelse(abs(theta) > inf, 1, 0)

  beta = theta[1:d]
  alpha_head = theta[(d + 1):(d + n - 1)]
  alpha = c(alpha_head, -sum(alpha_head))

  res = y - (z %*% alpha) - (x %*% beta)

  return(list(alpha = alpha, beta = beta, lambda = lambda, res = res, phi = phi))
}


#' @title Degrees of freedom
#'
#' @description This function estimates the degrees of freedom
#'
#' @param theta Numeric vector
#' @param N Numeric scalar, sample size
#' @param p Numeric scalar, parameter dimension
#' @param inf Numeric scalar, threshold parameter
#'
#' @keywords internal
#' @noRd
df_hat = function(theta,N,p,inf){
  s = sum(ifelse(abs(theta)<inf,0,1))
  m = min(c(s,N,p))
  return(m)
}

#' @title Bayesian Information Criteria
#'
#' @param res Numeric vector, residuals.
#' @param theta Numeric vector, parameters.
#' @param tau Numeric scalar, the percentile.
#' @param N Numeric integer, sample size.
#' @param p Numeric integer, parameter length.
#' @param inf Numeric, internal small quantity.
#'
#' @keywords internal
#' @noRd
bic_hat = function(res,theta,tau,N,p,inf){
  df  = df_hat(theta,N,p,inf)
  rho = rho_koenker(res, tau)
  ss  = sum(rho)
  bi  = ss +log(N)*df
  return(bi)
}

#' @title Incident matrix Z
#'
#' @description Create an Incident matrix Z
#'
#' @param n Numeric integer, number of incidents (subjects, units or individuals).
#' @param N Numeric integer, sample size.
#' @param id Numeric vector of integer, incident identification.
#'
#' @keywords internal
#' @noRd
make_z = function(n,N,id){
  z = matrix(0, nrow = N, ncol = n)
  j = 0
  for(i in 1:N){
    j = id[i]
    z[i,j] = 1
  }
  return(z)
}

#' @title Kernel density
#'
#' @param x Numeric vector.
#' @param inf Numeric, internal small quantity.
#'
#' @import stats
#'
#' @keywords internal
#' @noRd
f_den = function(x, inf) {
  dh = density(x, kernel = "epanechnikov")
  dx = dh$x
  dy = dh$y
  # Interpolate density values at points x
  y = approx(dx, dy, xout = x, rule = 2)$y
  # Replace any values below 'inf' with 'inf'
  y[y < inf] = inf
  return(y)
}


#' @title Covariance
#'
#' @description Estimate Covariance matrix
#'
#' @param alpha Numeric vector.
#' @param beta Numeric vector.
#' @param d length of beta.
#' @param inf Numeric scalar, internal value, small value.
#' @param n length of alpha.
#' @param N sample size.
#' @param res Numeric vector, residuals.
#' @param method Factor, "qr" quantile regression, "qrfe" quantile regression with fixed effects, "lqrfe" Lasso quantile regression with fixed effects, "alqr" adaptive Lasso quantile regression with fixed effects.
#' @param tau Numeric, identifies the percentile.
#' @param X Numeric matrix, covariates.
#' @param Z Numeric matrix, incident matrix.
#'
#' @keywords internal
#' @noRd
q_cov = function(alpha, beta, d, inf, n, N, res, method, tau, X, Z){
  omega = tau * (1 - tau)

  if(method %in% c("lqrfe", "alqrfe")){
    d_old = d
    beta_ind = ifelse(beta == 0, 0, 1)
    d = sum(beta_ind)
    X = X[, which(beta_ind == 1), drop=FALSE]
  }

  zz  = n * t(Z) %*% Z
  zx  = sqrt(n) * t(Z) %*% X
  xz  = t(zx)
  xx  = t(X) %*% X

  D0  = (omega / N) * rbind(cbind(zz, zx), cbind(xz, xx))

  fij = f_den(res, inf) + inf  # consider removing jitter for stability
  Phi = diag(fij)

  zpz = n * t(Z) %*% Phi %*% Z
  zpx = sqrt(n) * t(Z) %*% Phi %*% X
  xpz = t(zpx)
  xpx = t(X) %*% Phi %*% X

  D1  = (1 / N) * rbind(cbind(zpz, zpx), cbind(xpz, xpx))

  D1inv = MASS::ginv(D1)

  if(method == "qr"){
    sig2_alpha = 0
    if(d == 1){
      sig2_beta = (omega / N) * (N^2) * solve(xpx) %*% xx %*% solve(xpx)
      sig2_beta = as.vector(sig2_beta)
    } else {
      sig2_beta = diag((omega / N) * (N^2) * solve(xpx) %*% xx %*% solve(xpx))
    }
  } else {
    Sigma = D1inv %*% D0 %*% D1inv
    sig2_alpha = diag(Sigma[1:n, 1:n])
    if(d == 1){
      sig2_beta = Sigma[n + 1, n + 1]
    } else {
      sig2_beta = diag(Sigma[(n + 1):(n + d), (n + 1):(n + d)])
    }

    if(method %in% c("lqrfe", "alqrfe")){
      sig2_beta_small = sig2_beta
      sig2_beta = rep(Inf, d_old)
      count = 1
      for(i in seq_len(d_old)){
        if(beta_ind[i] == 1){
          sig2_beta[i] = sig2_beta_small[count]
          count = count + 1
        }
      }
    }
  }

  return(list(sig2_alpha = sig2_alpha, sig2_beta = sig2_beta))
}


#' @title Tabular function
#'
#' @param N sample size.
#' @param n length of alpha.
#' @param d length of beta.
#' @param theta Numeric vector.
#' @param sig2 Numeric vector.
#' @param kind Numeric, 1 means alpha, 2 means beta
#' @param inf Numeric scalar, internal value, small value.
#' @param digt Numeric integer, round.
#'
#' @import stats
#'
#' @keywords internal
#' @noRd
f_tab = function(N, n, d, theta, sig2, kind, inf, digt){
  m = N / n            # Observations per group
  len = qnorm(0.975)  # 1.96 for 95% CI
  p = length(theta)

  # Replace invalid or NA variances with Inf to avoid sqrt errors
  sig2[is.na(sig2) | sig2 <= 0] <- inf

  # Compute standard errors depending on 'kind'
  if(kind == 1) SE = sqrt(sig2 / m)   # kind==1: for alphas
  if(kind == 2) SE = sqrt(sig2 / N)   # kind==2: for betas

  infb = theta - len * SE   # Lower 95% CI bound
  supb = theta + len * SE   # Upper 95% CI bound
  zval = theta / SE
  pval = 2 * pnorm(abs(zval), lower.tail = FALSE)

  # Compute significance stars (sgf presumably user-defined)
  sig0 = sgf(as.vector(pval))

  # Create matrix and then data frame with results, rounding as needed
  Matx = cbind(theta, SE, infb, supb, zval, pval)
  Core = data.frame(round(Matx, digt), Sig = sig0)
  colnames(Core) = c("Coef", "Std. Error", "Inf CI95%", "Sup CI95%", "z value", "Pr(>|z|)", "Sig")

  # Assign row names depending on kind and length
  if(kind == 1){
    if(p == 1) rownames(Core) <- "alpha 1"
    else       rownames(Core) <- paste("alpha", 1:p)
  } else if(kind == 2){
    if(d == 1) rownames(Core) <- "beta 1"
    else       rownames(Core) <- paste("beta", 1:p)
  }

  return(list(Core = Core, Matx = Matx))
}


#' @title Identify significance
#'
#' @param x Numeric vector.
#'
#' @keywords internal
#' @noRd
sgf = function(x){
  m = length(x)
  y = rep(" ", m)
  for(i in 1:m){
    if(is.na(x[i]))x[i]=1
    if(x[i]<0.001) y[i] = "***"
    if(x[i]>=0.001 && x[i]<0.01) y[i] = "**"
    if(x[i]>=0.01 && x[i]<0.05) y[i] = "*"
    if(x[i]>=0.05 && x[i]<0.1) y[i] = "."
  }
  return(y)
}

#' @title Clean missings
#'
#' @param y Numeric vector, outcome.
#' @param x Numeric matrix, covariates
#' @param id Numeric vector, identifies the unit to which the observation belongs.
#'
#' @returns list with the same objects y, x, id, but without missings.
#'
#' @examples
#' n = 10
#' m = 4
#' d = 3
#' N = n*m
#' L = N*d
#' x = matrix(rnorm(L), ncol=d, nrow=N)
#' subj = rep(1:n, each=m)
#' alpha = rnorm(n)
#' beta = rnorm(d)
#' eps = rnorm(N)
#' y = x %*% beta  + matrix(rep(alpha, each=m) + eps)
#' y = as.vector(y)
#' x[1,3] = NA
#' clean_data(y=y, x=x, id=subj)
#'
#' @import stats
#'
#' @export
clean_data = function(y, x, id) {
  # Combine into data frame
  mega = data.frame(y = y, id = id, x)
  # Remove rows with any NA
  mega = na.omit(mega)
  # Sort by id
  mega = mega[order(mega$id), ]
  # Renumber id to consecutive integers
  mega$id = as.integer(factor(mega$id, levels = unique(mega$id)))

  y_new = mega$y
  x_new = as.matrix(mega[, -(1:2)])  # drop y and id columns
  id_new = mega$id

  return(list(y = y_new, x = x_new, id = id_new))
}


#' @title multiple penalized quantile regression
#'
#' @description Estimate QR for several taus
#'
#' @param y Numeric vector, outcome.
#' @param x Numeric matrix, covariates
#' @param subj Numeric vector, identifies the unit to which the observation belongs.
#' @param tau Numeric vector, identifies the percentiles.
#' @param method Factor, "qr" quantile regression, "qrfe" quantile regression with fixed effects, "lqrfe" Lasso quantile regression with fixed effects, "alqr" adaptive Lasso quantile regression with fixed effects.
#' @param ngrid Numeric scalar greater than one, number of BIC to test.
#' @param inf Numeric scalar, internal value, small value.
#' @param digt Numeric scalar, internal value greater than one, define "zero" coefficient.
#'
#' @return Beta Numeric array, with three dimmensions: 1) tau, 2) coef., lower bound, upper bound, 3) exploratory variables.
#'
#' @examples
#' n = 10
#' m = 5
#' d = 4
#' N = n*m
#' L = N*d
#' x = matrix(rnorm(L), ncol=d, nrow=N)
#' subj = rep(1:n, each=m)
#' alpha = rnorm(n)
#' beta = rnorm(d)
#' eps = rnorm(N)
#' y = x %*% beta  + matrix(rep(alpha, each=m) + eps)
#' y = as.vector(y)
#'
#' Beta = mqr(x,y,subj,tau=1:9/10, method="qr", ngrid = 10)
#' Beta
#'
#' @export
mqr = function(x,y,subj,tau=1:9/10, method="qr", ngrid = 20, inf = 1e-08, digt=4){
  ntau   = length(tau)
  d      = ifelse(is.null(dim(x)[2]), 1, dim(x)[2])
  Beta   = array(dim = c(ntau, 3, d))
  for(i in 1:ntau){
    Est = qr(x,y,subj,tau[i], method, ngrid, inf, digt)
    Beta[i,1,] = Est$Mat_beta[,1]
    Beta[i,2,] = Est$Mat_beta[,3]
    Beta[i,3,] = Est$Mat_beta[,4]
  }
  return(Beta)
}

#' @title multiple penalized quantile regression - alpha
#'
#' @description  Estimate QR intercepts for several taus
#'
#' @param y Numeric vector, outcome.
#' @param x Numeric matrix, covariates
#' @param subj Numeric vector, identifies the unit to which the observation belongs.
#' @param tau Numeric vector, identifies the percentiles.
#' @param method Factor, "qr" quantile regression, "qrfe" quantile regression with fixed effects, "lqrfe" Lasso quantile regression with fixed effects, "alqr" adaptive Lasso quantile regression with fixed effects.
#' @param ngrid Numeric scalar greater than one, number of BIC to test.
#' @param inf Numeric scalar, internal value, small value.
#' @param digt Numeric scalar, internal value greater than one, define "zero" coefficient.
#'
#' @return Alpha Numeric array, with three dimmensions: 1) tau, 2) coef., lower bound, upper bound, 3) exploratory variables.
#'
#' @examples
#' n = 10
#' m = 5
#' d = 4
#' N = n*m
#' L = N*d
#' x = matrix(rnorm(L), ncol=d, nrow=N)
#' subj = rep(1:n, each=m)
#' alpha = rnorm(n)
#' beta = rnorm(d)
#' eps = rnorm(N)
#' y = x %*% beta  + matrix(rep(alpha, each=m) + eps)
#' y = as.vector(y)
#'
#' Alpha = mqr(x,y,subj,tau=1:9/10, method="qr", ngrid = 10)
#' Alpha
#'
#' @export
mqr_alpha = function(x,y,subj,tau=1:9/10,  method="qr", ngrid = 20, inf = 1e-08, digt=4){
  ntau   = length(tau)
  Est0   = qr(x,y,subj,0.5, method, ngrid, inf, digt)
  n      = length(Est0$Mat_alpha[,1])
  Alpha  = array(dim = c(ntau, 3, n))
  for(i in 1:ntau){
    Est = qr(x,y,subj,tau[i], method, ngrid, inf, digt)
    Alpha[i,1,] = Est$Mat_alpha[,1]
    Alpha[i,2,] = Est$Mat_alpha[,3]
    Alpha[i,3,] = Est$Mat_alpha[,4]
  }
  return(Alpha)
}

#' @title plot multiple penalized quantile regression
#'
#' @description  plot QR for several taus
#'
#' @param Beta Numeric array, with three dimmensions: 1) tau, 2) coef., lower bound, upper bound, 3) exploratory variables.
#' @param tau Numeric vector, identifies the percentiles.
#' @param D covariate's number.
#' @param col color.
#' @param lwd line width.
#' @param lty line type.
#' @param pch point character.
#' @param cex.axis cex axis length.
#' @param cex.lab cex axis length.
#' @param main title.
#'
#' @examples
#' n = 10
#' m = 5
#' d = 4
#' N = n*m
#' L = N*d
#' x = matrix(rnorm(L), ncol=d, nrow=N)
#' subj = rep(1:n, each=m)
#' alpha = rnorm(n)
#' beta = rnorm(d)
#' eps = rnorm(N)
#' y = x %*% beta  + matrix(rep(alpha, each=m) + eps)
#' y = as.vector(y)
#'
#' Beta = mqr(x,y,subj,tau=1:9/10, method="qr", ngrid = 10)
#' plot_taus(Beta,tau=1:9/10,D=1)
#'
#' @export
plot_taus = function(Beta, tau=1:9/10, D, col=2, lwd=1, lty=2, pch=1, cex.axis=1, cex.lab=1, main=""){
  ntau  = dim(Beta)[1]
  d     = dim(Beta)[3]
  Beta  = matrix(Beta[,,D], ncol = 3, nrow = ntau)
  Mbeta = max(Beta)
  mbeta = min(Beta)
  Mtau  = max(tau)
  mtau  = min(tau)
  graphics::plot(c(mtau,Mtau),c(mbeta, Mbeta), xlab=expression(tau), ylab=expression(paste(beta,"(",tau,")")), main=main, type="n", cex.axis=cex.axis, cex.lab=cex.lab)
  graphics::polygon(c(tau,tau[ntau:1]), c(Beta[,2],Beta[ntau:1,3]), col="gray90", border = NA)
  graphics::lines(tau, Beta[,1], col=col, lty=lty, lwd=lwd)
  graphics::lines(tau, Beta[,1], col=col, type = "p", pch=pch)
  graphics::lines(tau, rep(0,ntau), lty=3, lwd=lwd)
}

#' @title plot multiple penalized quantile regression - alpha
#'
#' @description plot QR intercepts for several taus
#'
#' @param Beta Numeric array, with three dimmensions: 1) tau, 2) coef., lower bound, upper bound, 3) exploratory variables.
#' @param tau Numeric vector, identifies the percentiles.
#' @param D intercept's number.
#' @param ylab y legend
#' @param col color.
#' @param lwd line width.
#' @param lty line type.
#' @param pch point character.
#' @param cex.axis cex axis length.
#' @param cex.lab cex axis length.
#' @param main title.
#'
#' @examples
#' n = 10
#' m = 5
#' d = 4
#' N = n*m
#' L = N*d
#' x = matrix(rnorm(L), ncol=d, nrow=N)
#' subj = rep(1:n, each=m)
#' alpha = rnorm(n)
#' beta = rnorm(d)
#' eps = rnorm(N)
#' y = x %*% beta  + matrix(rep(alpha, each=m) + eps)
#' y = as.vector(y)
#'
#' Beta = mqr_alpha(x,y,subj,tau=1:9/10, method="qr", ngrid = 10)
#' plot_alpha(Beta,tau=1:9/10,D=1)
#'
#' @export
plot_alpha = function(Beta, tau=1:9/10, D,ylab=expression(alpha[1]), col=2, lwd=1, lty=2, pch=1, cex.axis=1, cex.lab=1, main=""){
  ntau  = dim(Beta)[1]
  d     = dim(Beta)[3]
  Beta  = matrix(Beta[,,D], ncol = 3, nrow = ntau)
  Mbeta = max(Beta)
  mbeta = min(Beta)
  Mtau  = max(tau)
  mtau  = min(tau)
  graphics::plot(c(mtau,Mtau),c(mbeta, Mbeta), xlab=expression(tau), ylab=ylab, main=main, type="n", cex.axis=cex.axis, cex.lab=cex.lab)
  graphics::polygon(c(tau,tau[ntau:1]), c(Beta[,2],Beta[ntau:1,3]), col="gray90", border = NA)
  graphics::lines(tau, Beta[,1], col=col, lty=lty, lwd=lwd)
  graphics::lines(tau, Beta[,1], col=col, type = "p", pch=pch)
  graphics::lines(tau, rep(0,ntau), lty=3, lwd=lwd)
}

