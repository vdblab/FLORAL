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
  
  for (int j = 0; j < p; ++j){
    
    arma::uvec idxy = arma::find(abs(beta) > 0);
    arma::uvec idxx = {j};
    
    double z = beta(j) + (xy(j) - dot(xx(idxx,idxy), beta(idxy)))/n;
    
    beta(j) = softthreshold(z, l);
    
  }
  
  return beta;
  
}


// [[Rcpp::export]]
arma::vec gd_cov_al(arma::mat xx, arma::vec xy, int n, double l, arma::vec beta, double mu, double alpha, bool intercept){
  
  int p = beta.n_elem;
  
  if (intercept){
    
    arma::uvec idb = arma::find(abs(beta) > 0);
    arma::uvec idxx = {0};
    
    beta(0) = ((xy(0) - dot(xx(idxx,idb), beta(idb)) + dot(xx(idxx,idxx), beta(idxx)) )/n)/(accu(xx(idxx,idxx))/n);
    
    for (int j = 1; j < p; ++j){
      
      idb = arma::find(abs(beta) > 0);
      idxx = {j};
      
      double z = (xy(j) - dot(xx(idxx,idb), beta(idb)) + dot(xx(idxx,idxx), beta(idxx)) )/n - mu*(accu(beta(idb)) - accu(beta(idxx)) + alpha);
      
      double denom = accu(xx(idxx,idxx))/n + mu;
      
      beta(j) = softthreshold(z, l)/denom;
      
    }
    
  }else{
    
    for (int j = 0; j < p; ++j){
      
      arma::uvec idb = arma::find(abs(beta) > 0);
      arma::uvec idxx = {j};
      
      double z = (xy(j) - dot(xx(idxx,idb), beta(idb)) + dot(xx(idxx,idxx), beta(idxx)) )/n - mu*(accu(beta(idb)) - accu(beta(idxx)) + alpha);
      
      double denom = accu(xx(idxx,idxx))/n + mu;
      
      beta(j) = softthreshold(z, l)/denom;
      
    }
    
  }
  
  return beta;
  
}

// [[Rcpp::export]]
Rcpp::List linear_lasso(arma::mat x, arma::vec y, int len){
  
  int n = y.n_elem;
  int p = x.n_cols;
  double lambda0 = log10(max(x.t() * y/n));
  double inc = (lambda0 + 2)/len;
  
  arma::mat beta;
  beta.zeros(p,len);
  arma::vec lambda = vec(len);
  arma::vec loss = vec(len);
  arma::vec mse = vec(len);
  
  arma::mat xx = x.t()*x;
  arma::vec xy = x.t()*y;
  
  for (int i=0; i<len; ++i){
    
    if (i >= 1) lambda0 = lambda0 - inc;
    double l = pow(10,lambda0);
    lambda(i) = l;
    
    arma::vec betai = beta.col(i);
    if (i >= 1) betai = beta.col(i-1);
    
    arma::vec diff = y - x*betai;
    double loss0 = accu(diff % diff)/(2*n) + l*accu(abs(betai));
    arma::vec betatemp = gd_cov(xx,xy,n,l,betai);
    arma::vec difftemp = y - x*betatemp;
    double lossnew = accu(difftemp % difftemp)/(2*n) + l*accu(abs(betatemp));
    
    while (abs(lossnew - loss0) > 1e-7){
      
      loss0 = lossnew;
      betai = betatemp;
      betatemp = gd_cov(xx,xy,n,l,betai);
      difftemp = y - x*betatemp;
      lossnew = accu(difftemp % difftemp)/(2*n) + l*accu(abs(betatemp));
      
    }
    
    loss(i) = lossnew;
    mse(i) = accu(difftemp % difftemp)/n;
    beta.col(i) = betatemp;
    
  }
  
  Rcpp::List ret;
  ret["beta"] = beta;
  ret["lambda"] = lambda;
  ret["loss"] = loss;
  ret["mse"] = mse;
  
  return ret;
  
}


// [[Rcpp::export]]
Rcpp::List linear_lasso_al(arma::mat x, arma::vec y, int len, double mu, int ub, arma::vec lambda, bool intercept){
  
  int n = y.n_elem;
  int p = x.n_cols;
  //double lambda0 = log10(max(x.t() * y/n));
  //double inc = (lambda0 + 2)/len;
  
  arma::mat beta;
  beta.zeros(p,len);
  //arma::vec lambda = vec(len);
  arma::vec loss = vec(len);
  arma::vec mse = vec(len);
  
  arma::mat xx = x.t()*x;
  arma::vec xy = x.t()*y;
  
  for (int i=0; i<len; ++i){
    
    //if (i >= 1) lambda0 = lambda0 - inc;
    //double l = pow(10,lambda0);
    
    double l = lambda(i);
    //lambda(i) = l;
    
    arma::vec betai = beta.col(i);
    if (i >= 1) betai = beta.col(i-1);
    
    double alpha = 0;
    int k = 0;
    
    arma::vec diff = y - x*betai;
    double loss0 = accu(diff % diff)/(2*n) + l*accu(abs(betai)) + mu*pow(accu(betai) + alpha,2)/2;
    arma::vec betatemp = gd_cov_al(xx,xy,n,l,betai,mu,alpha,intercept);
    arma::vec difftemp = y - x*betatemp;
    double lossnew = accu(difftemp % difftemp)/(2*n) + l*accu(abs(betatemp)) + mu*pow(accu(betatemp) + alpha,2)/2;
    
    while (abs(lossnew - loss0) > 1e-7){
      loss0 = lossnew;
      betatemp = gd_cov_al(xx,xy,n,l,betatemp,mu,alpha,intercept);
      difftemp = y - x*betatemp;
      lossnew = accu(difftemp % difftemp)/(2*n) + l*accu(abs(betatemp)) + mu*pow(accu(betatemp) + alpha,2)/2;
    }
    
    while (mean(abs(betatemp - betai)) > 1e-7 and k < ub){
      
      k = k+1;
      // Rcpp::Rcout << "k: " << k << endl;
      
      betai = betatemp;
      alpha = alpha + accu(betai);
      
      diff = y - x*betai;
      loss0 = accu(diff % diff)/(2*n) + l*accu(abs(betai)) + mu*pow(accu(betai) + alpha,2)/2;
      betatemp = gd_cov_al(xx,xy,n,l,betai,mu,alpha,intercept);
      difftemp = y - x*betatemp;
      lossnew = accu(difftemp % difftemp)/(2*n) + l*accu(abs(betatemp)) + mu*pow(accu(betatemp) + alpha,2)/2;
      
      while (abs(lossnew - loss0) > 1e-7){
        loss0 = lossnew;
        betatemp = gd_cov_al(xx,xy,n,l,betatemp,mu,alpha,intercept);
        difftemp = y - x*betatemp;
        lossnew = accu(difftemp % difftemp)/(2*n) + l*accu(abs(betatemp)) + mu*pow(accu(betatemp) + alpha,2)/2;
      }
    }
    
    loss(i) = lossnew;
    mse(i) = accu(difftemp % difftemp)/n;
    beta.col(i) = betatemp;
    
  }
  
  Rcpp::List ret;
  ret["beta"] = beta;
  ret["lambda"] = lambda;
  ret["loss"] = loss;
  ret["mse"] = mse;
  
  return ret;
  
}


// [[Rcpp::export]]
Rcpp::List logistic_lasso_al(arma::mat x, arma::vec y, int len, double mu, int ub, arma::vec lambda){
  
  int n = y.n_elem;
  int p = x.n_cols;
  
  arma::vec beta0;
  beta0.zeros(len);
  arma::mat beta;
  beta.zeros(p,len);
  arma::vec loss = vec(len);
  arma::vec mse = vec(len);
  
  for (int i=0; i<len; ++i){
    
    double l = lambda(i);
    
    double beta0i = beta0(i);
    arma::vec betai = beta.col(i);
    if (i >= 1) {
      beta0i = beta0(i-1);
      betai = beta.col(i-1);
    }
    
    double alpha = 0;
    int k = 0;
    
    // Rcpp::Rcout << "i: " << i<< endl;
    
    arma::vec prob = exp(beta0i + x*betai)/(1 + exp(beta0i + x*betai));
    arma::vec sfun = y - prob;
    arma::vec hfun = prob % (1-prob) + 1e-8;
    arma::vec z0 = beta0i + x*betai + sfun/hfun;
    arma::vec z = z0 - mean(z0);
    arma::mat xx = x.t()*diagmat(hfun)*x;
    arma::vec xz = x.t()*(z % hfun);
    
    arma::vec diff = z - x*betai;
    double loss0 = accu(hfun % diff % diff)/(2*n) + l*accu(abs(betai)) + mu*pow(accu(betai) + alpha,2)/2;
    arma::vec betatemp = gd_cov_al(xx,xz,n,l,betai,mu,alpha,false);
    arma::vec difftemp = z - x*betatemp;
    double lossnew = accu(hfun % difftemp % difftemp)/(2*n) + l*accu(abs(betatemp)) + mu*pow(accu(betatemp) + alpha,2)/2;
    
    while (abs(lossnew - loss0) > 1e-7){
      loss0 = lossnew;
      betatemp = gd_cov_al(xx,xz,n,l,betatemp,mu,alpha,false);
      difftemp = z - x*betatemp;
      lossnew = accu(hfun % difftemp % difftemp)/(2*n) + l*accu(abs(betatemp)) + mu*pow(accu(betatemp) + alpha,2)/2;
    }
    
    while (mean(abs(betatemp - betai)) > 1e-7 and k < ub){
      
      k = k+1;
      // Rcpp::Rcout << "k: " << k << endl;
      
      betai = betatemp;
      alpha = alpha + accu(betai);
      
      diff = z - x*betai;
      loss0 = accu(hfun % diff % diff)/(2*n) + l*accu(abs(betai)) + mu*pow(accu(betai) + alpha,2)/2;
      betatemp = gd_cov_al(xx,xz,n,l,betai,mu,alpha,false);
      difftemp = z - x*betatemp;
      lossnew = accu(hfun % difftemp % difftemp)/(2*n) + l*accu(abs(betatemp)) + mu*pow(accu(betatemp) + alpha,2)/2;
      
      while (abs(lossnew - loss0) > 1e-7){
        loss0 = lossnew;
        betatemp = gd_cov_al(xx,xz,n,l,betatemp,mu,alpha,false);
        difftemp = z - x*betatemp;
        lossnew = accu(hfun % difftemp % difftemp)/(2*n) + l*accu(abs(betatemp)) + mu*pow(accu(betatemp) + alpha,2)/2;
      }
    }
    
    beta0(i) = mean(z0 - x*betatemp);
    loss(i) = lossnew;
    mse(i) = accu(difftemp % difftemp)/n;
    beta.col(i) = betatemp;
    
  }
  
  Rcpp::List ret;
  ret["beta"] = beta;
  ret["beta0"] = beta0;
  ret["lambda"] = lambda;
  ret["loss"] = loss;
  ret["mse"] = mse;
  
  return ret;
  
}


// [[Rcpp::export]]
Rcpp::List cox_lasso_al(arma::mat x, arma::vec t, arma::vec d, arma::vec tj, int len, double mu, int ub, arma::vec lambda){
  
  int n = t.n_elem;
  int p = x.n_cols;
  int m = tj.n_elem;
  
  arma::mat beta;
  beta.zeros(p,len);
  arma::vec loss = vec(len);
  arma::vec mse = vec(len);
  arma::vec loglik;
  loglik.zeros(len);
  
  for (int i=0; i<len; ++i){
    
    double l = lambda(i);
    
    arma::vec betai = beta.col(i);
    if (i >= 1) {
      betai = beta.col(i-1);
    }
    
    double alpha = 0;
    int k = 0;
    int k1 = 0;
    
    // Rcpp::Rcout << "i: " << i<< endl;
    
    arma::vec link = exp(x*betai);
    arma::vec denomj = vec(m);
    arma::vec d1;
    d1.zeros(n);
    arma::vec d2;
    d2.zeros(n);
    
    arma::uvec idx = find(t >= tj(0));
    denomj(0) = accu(link(idx)) + 1e-8;
    
    arma::uvec tidx;
    arma::uvec widx = find(t >= tj(0));
    arma::uvec uidx = find(t == tj(0));
    
    if (i >= 1) {
      loglik(i-1) = loglik(i-1) + accu(x.rows(uidx)*betai) - log(accu(link(widx))+1e-8);
    }
    d1(widx) = d1(widx) + 1/denomj(0);
    d2(widx) = d2(widx) + 1/pow(denomj(0),2);
    
    for (int j=1; j<m; ++j){
      idx = find(t >= tj(j-1) && t < tj(j));
      denomj(j) = denomj(j-1) - accu(link(idx)) + 1e-8;
      
      widx = find(t >= tj(j));
      if (i >= 1){
        uidx = find(t == tj(j));
        loglik(i-1) = loglik(i-1) + accu(x.rows(uidx)*betai) - log(accu(link(widx))+1e-8);
      }
      
      d1(widx) = d1(widx) + 1/denomj(j);
      d2(widx) = d2(widx) + 1/pow(denomj(j),2);
      
    }
    
    // denomj = denomj + 1e-8;
    
    // for (int k=0; k<n; ++k){
    // 
    //   arma::uvec tidx = find(tj <= t(k));
    //   d1(k) = accu(1/denomj(tidx));
    //   d2(k) = accu(1/pow(denomj(tidx),2));
    // 
    // }
    
    arma::vec sfun = d - link % d1;
    arma::vec hfun = link % d1 - d2 % pow(link,2) + 1e-8;
    arma::vec z =  x*betai + sfun/hfun;
    arma::mat xx = x.t()*diagmat(hfun)*x;
    arma::vec xz = x.t()*(z % hfun);
    
    arma::vec diff = z - x*betai;
    double loss0 = accu(hfun % diff % diff)/(2*n) + l*accu(abs(betai)) + mu*pow(accu(betai) + alpha,2)/2;
    arma::vec betatemp = gd_cov_al(xx,xz,n,l,betai,mu,alpha,false);
    arma::vec difftemp = z - x*betatemp;
    double lossnew = accu(hfun % difftemp % difftemp)/(2*n) + l*accu(abs(betatemp)) + mu*pow(accu(betatemp) + alpha,2)/2;
    
    while (abs(lossnew - loss0) > 1e-7){
      loss0 = lossnew;
      betatemp = gd_cov_al(xx,xz,n,l,betatemp,mu,alpha,false);
      difftemp = z - x*betatemp;
      lossnew = accu(hfun % difftemp % difftemp)/(2*n) + l*accu(abs(betatemp)) + mu*pow(accu(betatemp) + alpha,2)/2;
    }
    
    while (mean(abs(betatemp - betai)) > 1e-7 and k < ub){
      
      k = k+1;
      betai = betatemp;
      alpha = alpha + accu(betai);
      
      diff = z - x*betai;
      loss0 = accu(hfun % diff % diff)/(2*n) + l*accu(abs(betai)) + mu*pow(accu(betai) + alpha,2)/2;
      betatemp = gd_cov_al(xx,xz,n,l,betai,mu,alpha,false);
      difftemp = z - x*betatemp;
      lossnew = accu(hfun % difftemp % difftemp)/(2*n) + l*accu(abs(betatemp)) + mu*pow(accu(betatemp) + alpha,2)/2;
      
      // k1 = 0;
      
      while (abs(lossnew - loss0) > 1e-7){
        // Rcpp::Rcout << "lossnew: " << lossnew << endl;
        // k1 = k1 + 1;
        loss0 = lossnew;
        betatemp = gd_cov_al(xx,xz,n,l,betatemp,mu,alpha,false);
        difftemp = z - x*betatemp;
        lossnew = accu(hfun % difftemp % difftemp)/(2*n) + l*accu(abs(betatemp)) + mu*pow(accu(betatemp) + alpha,2)/2;
      }
    }
    
    loss(i) = lossnew;
    mse(i) = accu(difftemp % difftemp)/n;
    beta.col(i) = betatemp;
    
    if (i == len-1){
      
      link = exp(x*betatemp);
      
      for (int j=0; j<m; ++j){
        widx = find(t >= tj(j));
        uidx = find(t == tj(j));
        loglik(i) = loglik(i) + accu(x.rows(uidx)*betai) - log(accu(link(widx))+1e-8);
      }
      
    }
    
  }
  
  Rcpp::List ret;
  ret["beta"] = beta;
  ret["lambda"] = lambda;
  ret["loss"] = loss;
  ret["mse"] = mse;
  ret["loglik"] = loglik;
  
  return ret;
  
}