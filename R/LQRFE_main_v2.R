#' @title quantile regression
#'
#' @description Estimate quantile regression with fixed effects for one tau
#'
#' @param x Numeric matrix, covariates
#' @param y Numeric vector, outcome.
#' @param subj Numeric vector, identifies the unit to which the observation belongs.
#' @param tau Numeric, identifies the percentile.
#' @param method Factor, "qr" quantile regression, "qrfe" quantile regression with fixed effects, "lqrfe" Lasso quantile regression with fixed effects, "alqr" adaptive Lasso quantile regression with fixed effects.
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
#' @importFrom MASS stats graphics Sitka89 Rabbit
#' @import Rcpp RcppArmadillo
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
qr =  function(x,y,subj,tau=0.5, method="qr", ngrid = 20, inf = 1e-08, digt=4){
  methods = c("qr","qrfe","lqrfe","alqrfe")
  if(any(method == methods)==FALSE){
    print("invalid method!")
  }
  inf2 = 1/10^digt 
  d  = ifelse(is.null(dim(x)[2]), 1, dim(x)[2])
  if(d==1) x = as.matrix(x)
  dclean = clean_data(y, x, id=subj)
  subj  = as.vector(dclean$id)
  y     = as.vector(dclean$y)
  N     = length(y)
  x     = as.matrix(dclean$x)
  beta  = as.vector(MASS::ginv(t(x)%*%x)%*%t(x)%*%y) 
  n     = max(subj)
  alpha = stats::rnorm(n)
  z     = make_z(n,N,id=subj)
  opt = optim_qr(beta,x,y,tau,N,d)
  beta = opt$beta 
  if(method=="qrfe"){
    opt = optim_qrfe(beta,alpha,x,y,z,tau,N,d,n)
  } 
  if(method=="lqrfe"){
    opt    = optim_lqr(beta,alpha,x,y,z,tau,N,d,n,ngrid,inf2)
  }
  if(method=="alqrfe"){
    opt = optim_qrfe(beta,alpha,x,y,z,tau,N,d,n)
    walpha = ifelse(abs(opt$alpha) <= inf2, inf2, abs(opt$alpha)) 
    wbeta  = ifelse(abs(opt$beta) <= inf2, inf2, abs(opt$beta))
    opt    = optim_alqr(beta,alpha,wbeta,walpha,x,y,z,tau,N,d,n,ngrid,inf2)
  }
  alpha  = ifelse(abs(opt$alpha)<=inf2, 0, opt$alpha) 
  beta   = ifelse(abs(opt$beta)<=inf2, 0, opt$beta) 
  res    = opt$res
  lambda = opt$lambda    
  Sigma  = q_cov(alpha, beta, d, inf, n, N, res, method, tau, X=x, Z=z)
  sig2_alpha = Sigma$sig2_alpha
  sig2_beta = Sigma$sig2_beta
  Tab_alpha = f_tab(N, n, d, alpha, sig2_alpha, 1, inf, digt)$Core
  Tab_beta  = f_tab(N, n, d, beta,  sig2_beta,  2, inf, digt)$Core
  Mat_alpha = f_tab(N, n, d, alpha,  sig2_alpha,  1, inf, digt)$Matx
  Mat_beta  = f_tab(N, n, d, beta,  sig2_beta,  2, inf, digt)$Matx
  obj = list(alpha=alpha, beta=beta, lambda=lambda, res=res,tau=tau, 
             c=c, sig2_alpha=sig2_alpha,sig2_beta=sig2_beta, Tab_alpha=Tab_alpha, 
             Tab_beta=Tab_beta, Mat_beta=Mat_beta, Mat_alpha=Mat_alpha, method=method)
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
#' @export
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
#' @return parametric vector and residuals.
#'
optim_qr = function(beta,x,y,tau,N,d){
  Opt    = stats::optim(par=beta, fn = loss_qr, method = "L-BFGS-B",  
                 N=N, y=y, x=x, tau=tau, d=d)
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
#' @return parametric vector and residuals
#' 
optim_qrfe = function(beta,alpha,x,y,z,tau,N,d,n){
  theta  = c(beta,alpha)
  Opt    = stats::optim(par=theta, fn = loss_qrfe, method = "L-BFGS-B",  
            n=N, y=y, x=x, z=z, tau=tau, d=d, mm=n)
  beta   = Opt$par[1:d]
  alpha  = Opt$par[(d+1):(d+n)]
  if(d==1) res = y -  (z %*% alpha) -  (x * beta)
  else     res = y -  (z %*% alpha) -  (x %*% beta) 
  return(list(alpha=alpha, beta=beta, lambda=0, res=res))
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
#' @return parametric vector and residuals
#'
optim_lqr = function(beta,alpha,x,y,z,tau,N,d,n,ngrid,inf){
  theta      = c(beta,alpha)
  p          = length(theta)
  lambda_min = sqrt(1/N)
  lambda_max = max(c(tau,1-tau))*(N/n)
  cgrid      = (log(lambda_max) - log(lambda_min))/(ngrid-1)
  lambda_tex = rep(lambda_min, ngrid)
  bic_tex    = rep(0, ngrid)
  theta_tex  = matrix(0, ncol=ngrid, nrow = p)
  for(t in 2:ngrid){
    lambda_tex[t] = lambda_tex[t-1]/exp(-cgrid)
  }
  for(t in 1:ngrid){
    Opt    = stats::optim(par=theta, fn = loss_lqr, method = "L-BFGS-B",   
                   n=N, y=y, x=x, z=z, tau=tau, d=d, mm=n,lambda=lambda_tex[t])
    theta_tex[,t] = Opt$par
    beta   = theta_tex[1:d,t]
    alpha  = theta_tex[(d+1):(d+n),t]
    if(d==1) res = y -  (z %*% alpha) -  (x * beta)
    else     res = y -  (z %*% alpha) -  (x %*% beta) 
    bic_tex[t] = bic_hat(res,theta=alpha,tau,N,p,inf)
  }
  maxw   = which.max(bic_tex)
  lambda = lambda_tex[maxw]
  theta  = theta_tex[,maxw]
  phi    = ifelse(abs(theta)<inf,0,1)
  beta   = theta[1:d]
  alpha  = theta[(d+1):(d+n)]
  if(d==1) res = y -  (z %*% alpha) -  (x * beta)
  else     res = y -  (z %*% alpha) -  (x %*% beta)   
  return(list(alpha=alpha, beta=beta, lambda=lambda, res=res, phi=phi))
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
#' @return parametric vector and residuals
#' 
optim_alqr = function(beta,alpha,wbeta,walpha,x,y,z,tau,N,d,n,ngrid,inf){
  theta      = c(beta,alpha)
  w          = c(wbeta,walpha)
  p          = length(theta)
  lambda_min = sqrt(1/N)
  lambda_max1 = max(c(tau,1-tau))*(N/n)
  lambda_max2 = max(abs(2*t(x)%*%y/N))
  lambda_max = min(c(lambda_max1,lambda_max2))
  cgrid      = (log(lambda_max) - log(lambda_min))/(ngrid-1)
  lambda_tex = rep(lambda_min, ngrid)
  bic_tex    = rep(0, ngrid)
  theta_tex  = matrix(0, ncol=ngrid, nrow = p)
  for(t in 2:ngrid){
    lambda_tex[t] = lambda_tex[t-1]/exp(-cgrid)
  }
  for(t in 1:ngrid){
    Opt    = stats::optim(par=theta, fn = loss_alqr, method = "L-BFGS-B",   
                   n=N, y=y, x=x, z=z, tau=tau, d=d, mm=n,lambda=lambda_tex[t], w=w)
    theta_tex[,t] = Opt$par
    beta   = theta_tex[1:d,t]
    alpha  = theta_tex[(d+1):(d+n),t]
    if(d==1) res = y -  (z %*% alpha) -  (x * beta)
    else     res = y -  (z %*% alpha) -  (x %*% beta) 
    bic_tex[t] = bic_hat(res,theta=theta_tex[,t],tau,N,p,inf)
  }
  maxw   = which.max(bic_tex)
  lambda = lambda_tex[maxw]
  theta  = theta_tex[,maxw]
  phi    = ifelse(abs(theta)<inf,0,1)
  beta   = theta[1:d]
  alpha  = theta[(d+1):(d+n)]
  if(d==1) res = y -  (z %*% alpha) -  (x * beta)
  else     res = y -  (z %*% alpha) -  (x %*% beta)   
  return(list(alpha=alpha, beta=beta, lambda=lambda, res=res, phi=phi))
}

#' @title degrees of fredom
#'
#' @description This function estimates the degrees of fredom
#' 
#' @param theta Numeric vector, parameters to be test
#' @param N Numeric integer, sample size.
#' @param p Numeric integer, length of theta.
#' @param inf Numeric, internal small quantity.
#' 
#' @return degrees of fredom
#' 
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
#' @return BIC value
#' 
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
#' @return Z matrix.
#' 
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
#' @return y vector, kernel density estimation.
#' 
#' @examples 
#' x = rnorm(10)
#' f_den(x, 0.0001)
#' 
f_den = function(x, inf){
  n  = length(x)
  dh = stats::density(x, kernel = "epanechnikov")
  dx = dh$x
  dy = dh$y
  M  = length(dx)
  y  = rep(max(min(dy),inf),n)
  for(i in 1:n){
    for(j in 2:(M-2)){
      if(x[i]>=dx[j-1] && x[i]<=dx[j+1]) y[i] = dy[j] 
    }
  }
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
q_cov = function(alpha, beta, d, inf, n, N, res, method, tau, X, Z){
  omega = tau*(1-tau)
  if(method=="lqrfe" || method=="alqrfe"){
    X_old = X
    d_old = d
    beta_ind = ifelse(beta==0, 0, 1) 
    d = sum(beta_ind)
    X = X[,c((1:d_old)*beta_ind)]
  }
  zz  = matrix(n * t(Z) %*% Z, nrow = n, ncol = n)
  zx  = matrix(sqrt(n) * t(Z) %*% X, nrow = n, ncol = d)
  xz  = matrix(sqrt(n) * t(X) %*% Z, nrow = d, ncol = n)
  xx  = matrix(t(X) %*% X, nrow = d, ncol = d) 
  D0  = drop(omega/N) * matrix(rbind(cbind(zz,zx), cbind(xz,xx)), nrow = (n+d), ncol = (n+d))
  fij = f_den(res, inf) + inf + stats::runif(N,inf,2*inf)
  Phi = diag(fij)
  zpz = matrix(n * t(Z) %*% Phi %*% Z, nrow = n, ncol = n)
  zpx = matrix(sqrt(n) * t(Z) %*% Phi %*% X, nrow = n, ncol = d)
  xpz = matrix(sqrt(n) * t(X) %*% Phi %*% Z, nrow = d, ncol = n)
  xpx = matrix(t(X) %*% Phi %*% X, nrow = d, ncol = d) 
  D1  = (1/N) * matrix(rbind(cbind(zpz,zpx), cbind(xpz,xpx)), nrow = (n+d), ncol = (n+d))
  D1inv = MASS::ginv(D1)
  if(method=="qr"){
    sig2_alpha = 0
    if(d==1) sig2_beta = (omega/N) * (N^2)*solve(xpx) %*% xx %*% solve(xpx)
    if(d> 1) sig2_beta = diag((omega/N) * (N^2)*solve(xpx) %*% xx %*% solve(xpx)) 
  }else{
    Sigma = D1inv %*% D0 %*% D1inv
    sig2_alpha = diag(Sigma[(1:n),(1:n)])
    if(d==1) sig2_beta  = Sigma[(n+1),(n+1)]
    if(d> 1) sig2_beta  = diag(Sigma[(n+1):(n+d),(n+1):(n+d)])
    if(method=="lqrfe" ||method=="alqrfe"){
      sig2_beta_small = sig2_beta
      sig2_beta = rep(Inf,d_old)
      c = 1
      for(i in 1:d_old){
        if(beta_ind[i]==1){
          sig2_beta[i] = sig2_beta_small[c]
          c = c+1
        }
      }
    }  
  }
  return(list(sig2_alpha=sig2_alpha, sig2_beta=sig2_beta))
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
f_tab = function(N, n, d, theta, sig2, kind, inf, digt){
  m = N/n
  len   = stats::qnorm(0.975)
  p     = length(theta)
  for(i in 1:p){
    if(is.na(sig2[i])) sig2[i] = inf
    if(sig2[i]<=0)     sig2[i] = inf
  }
  if(kind==1) SE = sqrt(sig2/(m))
  if(kind==2) SE = sqrt(sig2/(N))  
  infb  = theta - (len * SE) 
  supb  = theta + (len * SE)
  zval  = theta/SE
  pval  = 2 * stats::pnorm(abs(theta/SE), lower.tail = F)
  sig0  = sgf(as.vector(pval))
  Matx  = matrix(cbind(theta, SE, infb, supb, zval, pval), ncol=6, nrow=p)
  Core  = data.frame(cbind(round(Matx,digt), sig0))
  colnames(Core) = c("Coef", "Std. Error", "Inf CI95%", "Sup CI95%", "z value", "Pr(>|z|)", "Sig")
  if(kind==1){
    if(p==1) rownames(Core)[1] = "alpha 1"
    if(p> 1) rownames(Core) = paste("alpha", 1:p)    
  } 
  if(kind==2){
    if(d==1) rownames(Core)[1] = "beta 1"
    if(d> 1) rownames(Core) = paste("beta", 1:p)
  }   
  return(list(Core=Core, Matx=Matx))
}

#' @title Identify significance
#'
#' @param x Numeric vector.
#' 
#' @return y vector Factor, symbol flag of significant p-values.
#' 
#' @examples 
#' n = 10
#' pvalue = rgamma(10,1,10)
#' sgf(pvalue) 
#'   
sgf = function(x){
  m = length(x)
  y = rep(" ", m)
  for(i in 1:m){
    if(is.na(x[i]))x[i]=0
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
clean_data = function(y, x, id){
  mega = cbind(y,id,x)
  mega = stats::na.omit(mega)
  mega = as.data.frame(mega)
  mega = mega[order(mega$id),]
  n    = max(mega$id)
  N    = length(mega$id)
  d    = dim(x)[2]
  it   = 1
  new  = 1
  iv   = rep(1,N)
  for(i in 1:N){
    if(mega$id[i]==it){
      iv[i] = new
    }else{
      it    = mega$id[i]
      new   = new+1
      iv[i] = new 
    }
  }    
  y = as.vector(mega$y)
  x = as.matrix(mega[,3:(d+2)])
  id = as.vector(iv)
  return(list(y=y,x=x,id=id))
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

