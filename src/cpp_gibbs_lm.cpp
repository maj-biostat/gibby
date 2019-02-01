// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; indent-tabs-mode: nil; -*-

// [[Rcpp::depends("RcppArmadillo")]]


#include <RcppArmadillo.h>
//#include <boost/random.hpp>
//#include <boost/math/distributions/inverse_gamma.hpp>

// not supposed to retain these if adding funcs to package...
//#include <math.h>
//#include <cmath>

const double log2pi = std::log(2.0 * M_PI);

using namespace Rcpp;


// [[Rcpp::export]]
arma::mat rcpp_gibbs_lm(const arma::colvec y, const arma::mat X, 
                 int iter = 2, 
                 double initphi = 5, 
                 double hyp_phi_alph = 5, double hyp_phi_gam = 8){

  // parameters for the gamma dist 
  double a_shape = 0;
  double b_rate = 0;
  
  
  arma::mat xprimexinv = arma::inv(X.t() * X);
  arma::mat phi = arma::zeros(iter, 1);
  arma::vec mu = arma::zeros(X.n_cols);
  arma::mat b = arma::zeros(iter, X.n_cols); 
  
  arma::mat chain = arma::zeros(iter, X.n_cols + 1); 

  
  // initial value for sigma (I convert into a matrix for generating mv norm)
  phi[0] = initphi; // arma::as_scalar(arma::randu(1));
  mu = arma::inv(X.t() * X) * (X.t() * y); 
  // Rcpp::Rcout << "mu " << mu <<  std::endl;
  
  for(int i = 1; i < iter; i++) {
    
    // see:
    // https://en.wikipedia.org/wiki/Multivariate_normal_distribution#Drawing_values_from_the_distribution
    // http://gallery.rcpp.org/articles/simulate-multivariate-normal/
    // multivariate normal - there is probably a single command to do this but not 
    // sure if will be as quicl so doing it from scratch  :/
    // b.row(i) = mu.t() + (arma::randn(1, X.n_cols)  * phi[i-1] * arma::eye(X.n_cols, X.n_cols));
    
    // hmm, there is: 
    
    // see page 207 (or 227 for multivariable case)
    // beta coefficients have covariance matrix sigma * (X'X)^{-1}
    // e = Y - hat(Y) = Y - HY = (I - H)Y
    // sigma is estimated by s2[e] = MSE(I - H) where H is the hat matrix.
    b.row(i) =  arma::trans(mvnrnd(mu, phi[i-1] * xprimexinv, 1));
    
    // Rcpp::Rcout << "sig "<< std::endl;
    // Rcpp::Rcout << phi[i-1] * xprimexinv << std::endl;
    // Rcpp::Rcout <<  std::endl;
    // Rcpp::Rcout << "b.row(" << i << ") "<< b.row(i)  << std::endl;
 
    a_shape = (y.n_rows/2) + hyp_phi_alph;
    b_rate = arma::as_scalar((0.5 * trans(y - X * b.row(i).t())*(y - X * b.row(i).t())) + hyp_phi_gam);
    
    // Rcpp::Rcout <<  std::endl;
    // Rcpp::Rcout << "a_shape " << a_shape << std::endl;
    // Rcpp::Rcout << "b_rate  " << b_rate << std::endl;
    
    phi[i] = 1/arma::as_scalar(arma::randg(1, arma::distr_param(a_shape, 1/b_rate)));

    // Rcpp::Rcout << "phi[" << i << "] "<< phi[i]  << std::endl;
    
  }
  
  // joins across rows 
  chain = arma::join_rows(b, phi); 

  return chain;
}






// [[Rcpp::export]]
arma::mat mvrnormArma(int n, arma::vec mu, arma::mat sigma) {
  int ncols = sigma.n_cols;
  arma::mat Y = arma::randn(n, ncols);
  return arma::repmat(mu, 1, n).t() + Y * arma::chol(sigma);
}


// [[Rcpp::export]]
arma::mat mygamma(int n = 10, 
                     double a_shape = 1, double b_rate = 1){
  
  
  return arma::randg(n, arma::distr_param(a_shape, 1/b_rate));
  
}


// [[Rcpp::export]]
arma::mat myinvgamma(int n = 10, 
                  double a_shape = 1, double b_rate = 1){
  
  // https://stats.stackexchange.com/questions/224714/sampling-from-an-inverse-gamma-distribution
  // http://www.math.yorku.ca/Who/Faculty/Madras/4281.dir/4281simgen.pdf
  
  return 1/arma::randg(n, arma::distr_param(a_shape, 1/b_rate));
  
}


