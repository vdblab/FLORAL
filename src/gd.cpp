// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; indent-tabs-mode: nil; -*-

#include "RcppArmadillo.h"
// [[Rcpp::depends(RcppArmadillo)]]
using namespace arma;
using namespace std;
// [[Rcpp::export]]

double softthreshold(double x, double lambda){
  
  double y = abs(x) - lambda;
  
  if (y <= 0) y = 0;
  
  y = sign(x)*y;
  
  return y;
  
}

// [[Rcpp::export]]
arma::vec gd_naive(arma::mat x, arma::vec y, double l, arma::vec beta){
  
  int p = beta.n_elem;
  int n = x.n_rows;
  
  for (int j = 0; j < p; ++j){
    
    arma::vec r = y - x*beta;
    
    double z = beta(j) + accu(dot(x.col(j),r))/n;
    
    beta(j) = softthreshold(z, l);
    
  }
  
  return beta;
  
}

// [[Rcpp::export]]
arma::vec gd_cov(arma::mat xx, arma::vec xy, int n, double l, arma::vec beta){
  
  int p = beta.n_elem;
  arma::uvec idxy = arma::find(abs(beta) > 0);
  
  for (int j = 0; j < p; ++j){
    
    arma::uvec idxx = {j};
    
    double z = beta(j) + (xy(j) - dot(xx(idxx,idxy), beta(idxy)))/n;
    
    beta(j) = softthreshold(z, l);
    
  }
  
  return beta;
  
}