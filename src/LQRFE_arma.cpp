#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
using namespace arma;

//' Rho Koenker (internal function)
 //'
 //' @param x Numeric vector
 //' @param tau Numeric scalar, percentile (between 0 and 1)
 //'
 //' @return A numeric vector with the Koenker rho function applied element-wise.
 //' @keywords internal
 //' @noRd
 // [[Rcpp::export]]
arma::vec rho_koenker(arma::vec x, double tau){
  int n = x.n_elem;
  arma::vec y(n);
  for(int i = 0; i < n; ++i){
    if(x(i)<0){
      y(i) = x(i)*(tau-1);
    } else {
      y(i) = x(i)*tau;
    }
  }
  return(y);
}

//' Loss quantile regression
//'
//' @param beta initial values
//' @param x design matrix
//' @param y vector output
//' @param tau percentile
//' @param N sample size
//' @param d columns of x
//'
//' @return Scalar quantile regression loss value.
//' @keywords internal
//' @noRd
// [[Rcpp::export]]
double loss_qr(arma::vec beta, arma::mat x, arma::vec y, double tau, int N, int d){
  double eta = 0;
  arma::vec res(N);
  arma::vec rho(N);
  res = y - (x * beta);
  rho = rho_koenker(res,tau);
  eta = accu(rho);
  return(eta);
}

//' Loss quantile regression with fixed effects
//'
//' @param theta initial values
//' @param x design matrix
//' @param y vector output
//' @param z incident matrix
//' @param tau percentile
//' @param n N sample size
//' @param d columns of x
//' @param mm n columns of z
//'
//' @return Scalar quantile regression loss value.
//' @keywords internal
//' @noRd
// [[Rcpp::export]]
double loss_qrfe(arma::vec theta,arma::mat x,arma::vec y,arma::mat z,double tau,int n,int d,int mm){
  double eta;
  arma::vec beta(d);
  arma::vec alpha(mm);
  arma::vec res(n);
  arma::vec rho(n);
  beta = theta.subvec(0, d - 1);
  alpha.head(mm - 1) = theta.subvec(d, d + mm - 2);     // alpha: next mm - 1 values
  alpha(mm - 1) = -arma::accu(alpha.head(mm - 1));      // sum-to-zero constraint
  res = y -  (z * alpha) -  (x * beta);
  rho = rho_koenker(res,tau);
  eta = accu(rho);
  return(eta);
}

//' Loss lasso quantile regression with fixed effects
//'
//' @param theta initial values
//' @param x design matrix
//' @param y vector output
//' @param z incident matrix
//' @param tau percentile
//' @param n N sample size
//' @param d columns of x
//' @param mm n columns of z
//' @param lambda constriction parameter
//'
//' @return Scalar loss value including quantile loss and L1 penalty on beta.
//' @keywords internal
//' @noRd
// [[Rcpp::export]]
double loss_lqr(arma::vec theta,arma::mat x,arma::vec y,arma::mat z,double tau,int n,int d,int mm, double lambda){
  double eta;
  arma::vec beta(d);
  arma::vec alpha(mm);
  arma::vec res(n);
  arma::vec rho(n);
  beta = theta.subvec(0, d - 1);
  alpha.head(mm - 1) = theta.subvec(d, d + mm - 2);     // alpha: next mm - 1 values
  alpha(mm - 1) = -arma::accu(alpha.head(mm - 1));      // sum-to-zero constraint
  res = y -  (z * alpha) -  (x * beta);
  rho = rho_koenker(res,tau);
  eta = accu(rho) + lambda * accu(abs(alpha));
  return(eta);
}

//' @title Loss Function for Adaptive Lasso Quantile Regression with Fixed Effects
//'
//' @description Computes the penalized quantile regression loss with adaptive lasso penalty
//' and fixed effects under a sum-to-zero constraint.
//'
//' @param theta Numeric vector. Parameters: beta (length d) and alpha[1:(mm-1)] (length mm - 1).
//' @param x Numeric matrix. Covariate design matrix (n x d).
//' @param y Numeric vector. Response vector (length n).
//' @param z Numeric matrix. Incidence matrix for fixed effects (n x mm).
//' @param tau Numeric scalar. Quantile level.
//' @param n Integer. Sample size.
//' @param d Integer. Number of covariates (columns in x).
//' @param mm Integer. Number of fixed effects (columns in z).
//' @param lambda Numeric. Lasso penalty parameter.
//' @param w Numeric vector. Adaptive lasso weights for all parameters (length d + mm).
//'
//' @return Loss scalar (quantile loss + adaptive lasso penalty)
//' @keywords internal
//' @noRd
// [[Rcpp::export]]
double loss_alqr(arma::vec theta,arma::mat x,arma::vec y,arma::mat z,double tau,int n,int d,int mm,double lambda,arma::vec w){
  double eta;
  double penalty;
  arma::vec beta(d);
  arma::vec alpha(mm);
  arma::vec res(n);
  arma::vec rho(n);
  arma::vec theta2(d + mm);

  beta = theta.subvec(0, d - 1);                  // First d are beta
  alpha.head(mm - 1) = theta.subvec(d, d + mm - 2);         // Fill first mm-1 alpha
  alpha(mm - 1) = -arma::accu(alpha.head(mm - 1));          // Sum-to-zero constraint
  theta2.head(d) = beta;
  theta2.tail(mm) = alpha;

  res = y - z * alpha - x * beta;
  rho = rho_koenker(res, tau);

  penalty = arma::accu(arma::abs(theta2 / w));       // element-wise division
  eta = arma::accu(rho) + lambda * penalty;
  return eta;
}


