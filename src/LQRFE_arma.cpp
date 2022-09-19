#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;

arma::vec fill_v(int a, int b, arma::vec x){
  int n = b - a + 1;
  arma::vec y(n);
  for(int i=0; i<n; i++){
    y(i) = x(i+a-1);
  }
  return(y);
}

//' Rho Koenker
//'
//' @param x generic vector
//' @param tau percentile
// [[Rcpp::export(rho_koenker)]]
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
// [[Rcpp::export(loss_qr)]]
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
// [[Rcpp::export(loss_qrfe)]]
double loss_qrfe(arma::vec theta,arma::mat x,arma::vec y,arma::mat z,double tau,int n,int d,int mm){
  double eta;
  arma::vec beta(d);
  arma::vec alpha(mm);
  arma::vec res(n);
  arma::vec rho(n);
  beta = fill_v(1,d, theta); 
  alpha = fill_v((d+1),d+mm, theta); 
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
// [[Rcpp::export(loss_lqr)]]
double loss_lqr(arma::vec theta,arma::mat x,arma::vec y,arma::mat z,double tau,int n,int d,int mm, double lambda){
  double eta;
  arma::vec beta(d);
  arma::vec alpha(mm);
  arma::vec res(n);
  arma::vec rho(n);
  beta = fill_v(1,d, theta); 
  alpha = fill_v(d+1,d+mm, theta); 
  res = y -  (z * alpha) -  (x * beta);   
  rho = rho_koenker(res,tau);
  eta = accu(rho) + lambda * accu(abs(alpha));
  return(eta);
}

//' Loss adaptive lasso quantile regression with fixed effects
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
//' @param w weights
// [[Rcpp::export(loss_alqr)]]
double loss_alqr(arma::vec theta,arma::mat x,arma::vec y,arma::mat z,double tau,
                int n,int d,int mm, double lambda, arma::vec w){
  double eta;
  arma::vec beta(d);
  arma::vec alpha(mm);
  arma::vec res(n);
  arma::vec rho(n);
  beta = fill_v(1,d, theta); 
  alpha = fill_v(d+1,d+mm, theta); 
  res = y -  (z * alpha) -  (x * beta);   
  rho = rho_koenker(res,tau);
  eta = accu(rho) + lambda * accu(abs(theta/w));
  return(eta);
}



