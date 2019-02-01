// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; indent-tabs-mode: nil; -*-


#include <RcppDist.h>

// [[Rcpp::depends(RcppArmadillo, RcppDist)]]



// [[Rcpp::export]]
arma::mat rcpp_gibbs_probit(const arma::colvec y, const arma::mat X,
                       int iter = 2){
  
  
  int n1 = sum(y);
  int n0 = y.n_rows - sum(y);
  
  // hyperparameters
  arma::vec theta_0 = arma::zeros(X.n_cols);
  arma::vec q_0 = arma::zeros(X.n_cols);
  q_0.fill(10);
  arma::mat Q_0 = arma::diagmat(q_0);
  
  arma::mat prec_0 = arma::zeros(X.n_cols, X.n_cols); 
  arma::mat V = arma::zeros(X.n_cols, X.n_cols); 
  arma::mat S = arma::zeros(X.n_cols, X.n_rows); 
  arma::vec M = arma::zeros(X.n_cols);
  double m = 0;
  
  
  // initial values for parameters
  arma::vec theta = arma::zeros(X.n_cols);
  arma::vec z = arma::zeros(X.n_rows);
  double z_old = 0;

  arma::vec h = arma::zeros(X.n_rows);
  arma::vec w = arma::zeros(X.n_rows);
  arma::vec u = arma::zeros(X.n_rows);
  
  arma::mat samp = arma::zeros(iter, X.n_cols); 
  
  
  // inverse
  prec_0 = arma::inv(Q_0);
  V = arma::inv(prec_0 + X.t() * X);
  
  // D x N matrix
  S = V * X.t();  
    

  for (int j = 0; j < X.n_rows; j++){
    h[j] = arma::as_scalar(X.row(j) * S.col(j));
    w[j] = h[j] / (1 - h[j]);
    u[j] = w[j] + 1;
  }
  
  
  for(int i = 0; i < X.n_rows; i++){
    
    if(y[i] == 0){
      z[i] = r_truncnorm(0, 1, R_NegInf, 0);
    } else {
      z[i] = r_truncnorm(0, 1, 0, R_PosInf);
    }
    
  }
  
  M = S * z;
  
  for(int t = 1; t < iter; t++){
    
    for(int j = 0; j < X.n_rows; j++){
      
      z_old = z[j];
      
      m = arma::as_scalar(X.row(j) * M);
      m = m - (w[j] * (z[j] - m));
      
      if(y[j] == 0){
        z[j] = r_truncnorm(m, arma::as_scalar(u[j]), R_NegInf, 0);
      } else {
        z[j] = r_truncnorm(m, arma::as_scalar(u[j]), 0, R_PosInf);
      }
      
      M = M + (z[j] - z_old) * S.col(j);
      
    }
    samp.row(t) = arma::trans(arma::mvnrnd(M, V)); 
  }

  
  return samp;
  
}

