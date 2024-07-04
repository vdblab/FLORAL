// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; indent-tabs-mode: nil; -*-

#include "RcppArmadillo.h"
// [[Rcpp::depends(RcppArmadillo)]]
using namespace arma;
using namespace std;
// [[Rcpp::depends(RcppProgress)]]
#include <progress.hpp>
#include <progress_bar.hpp>

// [[Rcpp::export]]
double softthreshold(double x, double lambda){
  
  double y = abs(x) - lambda;
  
  if (y <= 0) y = 0;
  
  y = sign(x)*y;
  
  return y;
  
}

// // [[Rcpp::export]]
// arma::vec gd_naive(arma::mat x, arma::vec y, double l, arma::vec beta){
//   
//   unsigned int p = beta.n_elem;
//   int n = x.n_rows;
//   
//   for (unsigned int j = 0; j < p; ++j){
//     
//     arma::vec r = y - x*beta;
//     
//     double z = beta(j) + accu(dot(x.col(j),r))/n;
//     
//     beta(j) = softthreshold(z, l);
//     
//   }
//   
//   return beta;
//   
// }

// [[Rcpp::export]]
arma::vec gd_cov(arma::mat xx, arma::vec xy, int n, double l, arma::vec beta){
  
  unsigned int p = beta.n_elem;
  
  for (unsigned int j = 0; j < p; ++j){
    
    arma::uvec idxy = arma::find(abs(beta) > 0);
    arma::uvec idxx = {j};
    
    double z = beta(j) + (xy(j) - dot(xx(idxx,idxy), beta(idxy)))/n;
    
    beta(j) = softthreshold(z, l);
    
  }
  
  return beta;
  
}


// [[Rcpp::export]]
arma::vec gd_cov_al(arma::mat xx, arma::vec xy, int n, double l, double a, arma::vec beta, double mu, double alpha, bool adjust=false, unsigned int ncov=0, double wcov=0){
  
  unsigned int p = beta.n_elem;
  
  if (adjust){
    
    for (unsigned int i = 0; i < ncov; ++i){
      
      arma::uvec idb = arma::find(abs(beta) > 0);
      arma::uvec idxx = {i};
      
      double z = (xy(i) - dot(xx(idxx,idb), beta(idb)) + dot(xx(idxx,idxx), beta(idxx)) )/n;
      double denom = accu(xx(idxx,idxx))/n + (1-a)*l*wcov;
      
      beta(i) = softthreshold(z, a*l*wcov)/denom;
      
    }
    
    for (unsigned int j = ncov; j < p; ++j){
      
      arma::uvec idb = arma::find(abs(beta) > 0);
      arma::uvec idxx = {j};
      // cout << idxx << endl;
      // idxx = j;
      // cout << idxx << endl;
      
      double z = (xy(j) - dot(xx(idxx,idb), beta(idb)) + dot(xx(idxx,idxx), beta(idxx)) )/n - mu*(accu(beta(idb)) - accu(beta(idxx)) - accu(beta.subvec(0,ncov-1)) + alpha);
      
      double denom = accu(xx(idxx,idxx))/n + (1-a)*l + mu;
      
      beta(j) = softthreshold(z, a*l)/denom;
      
    }
    
  }else{
    
    for (unsigned int j = 0; j < p; ++j){
      
      arma::uvec idb = arma::find(abs(beta) > 0);
      arma::uvec idxx = {j};
      
      double z = (xy(j) - dot(xx(idxx,idb), beta(idb)) + dot(xx(idxx,idxx), beta(idxx)) )/n - mu*(accu(beta(idb)) - accu(beta(idxx)) + alpha);
      
      double denom = accu(xx(idxx,idxx))/n + (1-a)*l + mu;
      
      beta(j) = softthreshold(z, a*l)/denom;
      
    }
    
  }
  
  return beta;
  
}

// [[Rcpp::export]]
Rcpp::List linear_enet_al(arma::mat x, arma::vec y, int len, double mu, int ub, arma::vec lambda, double wcov, double a, bool adjust, unsigned int ncov, bool display_progress=true){
  
  int n = y.n_elem;
  unsigned int p = x.n_cols;
  //double lambda0 = log10(max(x.t() * y/n));
  //double inc = (lambda0 + 2)/len;
  
  arma::mat beta;
  beta.zeros(p,len);
  //arma::vec lambda = vec(len);
  arma::vec loss = vec(len);
  arma::vec mse = vec(len);
  arma::vec tol = vec(len);
  
  arma::mat xx = x.t()*x;
  arma::vec xy = x.t()*y;
  
  Progress prog(len, display_progress);
  
  for (int i=0; i<len; ++i){
    
    prog.increment();
    
    //if (i >= 1) lambda0 = lambda0 - inc;
    //double l = pow(10,lambda0);
    
    double l = lambda(i);
    //lambda(i) = l;
    
    arma::vec betai = beta.col(i);
    if (i >= 1) betai = beta.col(i-1);
    
    double alpha = 0;
    int k = 0;
    int k1 = 0;
    
    arma::vec diff = y - x*betai;

    double loss0 = 0;
    if (adjust){
      loss0 = accu(diff % diff)/(2*n)+ a*l*accu(abs(betai.subvec(ncov,p-1))) + (1-a)*l*accu(pow(betai.subvec(ncov,p-1),2))/2 + mu*pow(accu(betai.subvec(ncov,p-1)) + alpha,2)/2 + wcov*a*l*accu(abs(betai.subvec(0,ncov-1))) + wcov*(1-a)*l*accu(pow(betai.subvec(0,ncov-1),2))/2;
    }else{
      loss0 = accu(diff % diff)/(2*n)+ a*l*accu(abs(betai.subvec(ncov,p-1))) + (1-a)*l*accu(pow(betai.subvec(ncov,p-1),2))/2 + mu*pow(accu(betai.subvec(ncov,p-1)) + alpha,2)/2;    
    }
    
    arma::vec betatemp = gd_cov_al(xx,xy,n,l,a,betai,mu,alpha,adjust,ncov,wcov);
    arma::vec difftemp = y - x*betatemp;
    
    double lossnew = 0;
    if (adjust){
      lossnew = accu(difftemp % difftemp)/(2*n)+ a*l*accu(abs(betatemp.subvec(ncov,p-1))) + (1-a)*l*accu(pow(betatemp.subvec(ncov,p-1),2))/2 + mu*pow(accu(betatemp.subvec(ncov,p-1)) + alpha,2)/2 + wcov*a*l*accu(abs(betatemp.subvec(0,ncov-1))) + wcov*(1-a)*l*accu(pow(betatemp.subvec(0,ncov-1),2))/2;
    }else{
      lossnew = accu(difftemp % difftemp)/(2*n)+ a*l*accu(abs(betatemp.subvec(ncov,p-1))) + (1-a)*l*accu(pow(betatemp.subvec(ncov,p-1),2))/2 + mu*pow(accu(betatemp.subvec(ncov,p-1)) + alpha,2)/2;
    }
    
    while (abs(lossnew - loss0) > 1e-7 and k1 < ub){
      // Rcpp::Rcout << "lossnew: " << lossnew << endl;
      k1 = k1 + 1;
      loss0 = lossnew;
      betatemp = gd_cov_al(xx,xy,n,l,a,betatemp,mu,alpha,adjust,ncov,wcov);
      difftemp = y - x*betatemp;
      
      if (adjust){
        lossnew = accu(difftemp % difftemp)/(2*n)+ a*l*accu(abs(betatemp.subvec(ncov,p-1))) + (1-a)*l*accu(pow(betatemp.subvec(ncov,p-1),2))/2 + mu*pow(accu(betatemp.subvec(ncov,p-1)) + alpha,2)/2 + wcov*a*l*accu(abs(betatemp.subvec(0,ncov-1))) + wcov*(1-a)*l*accu(pow(betatemp.subvec(0,ncov-1),2))/2;
      }else{
        lossnew = accu(difftemp % difftemp)/(2*n)+ a*l*accu(abs(betatemp.subvec(ncov,p-1))) + (1-a)*l*accu(pow(betatemp.subvec(ncov,p-1),2))/2 + mu*pow(accu(betatemp.subvec(ncov,p-1)) + alpha,2)/2;
      }
    }
    
    while (mean(abs(betatemp - betai)) > 1e-7 and k < ub){
      
      k = k+1;
      // Rcpp::Rcout << "k: " << k << endl;
      
      betai = betatemp;
      alpha = alpha + accu(betai.subvec(ncov,p-1));
      
      diff = y - x*betai;
      loss0 = accu(diff % diff)/(2*n) + l*accu(abs(betai)) + mu*pow(accu(betai) + alpha,2)/2;
      betatemp = gd_cov_al(xx,xy,n,l,a,betai,mu,alpha,adjust,ncov,wcov);
      difftemp = y - x*betatemp;
      
      if (adjust){
        lossnew = accu(difftemp % difftemp)/(2*n)+ a*l*accu(abs(betatemp.subvec(ncov,p-1))) + (1-a)*l*accu(pow(betatemp.subvec(ncov,p-1),2))/2 + mu*pow(accu(betatemp.subvec(ncov,p-1)) + alpha,2)/2 + wcov*a*l*accu(abs(betatemp.subvec(0,ncov-1))) + wcov*(1-a)*l*accu(pow(betatemp.subvec(0,ncov-1),2))/2;
      }else{
        lossnew = accu(difftemp % difftemp)/(2*n)+ a*l*accu(abs(betatemp.subvec(ncov,p-1))) + (1-a)*l*accu(pow(betatemp.subvec(ncov,p-1),2))/2 + mu*pow(accu(betatemp.subvec(ncov,p-1)) + alpha,2)/2;
      }
      
      k1 = 0;
      
      while (abs(lossnew - loss0) > 1e-7 and k1 < ub){
        // Rcpp::Rcout << "lossnew: " << lossnew << endl;
        k1 = k1 + 1;
        loss0 = lossnew;
        betatemp = gd_cov_al(xx,xy,n,l,a,betatemp,mu,alpha,adjust,ncov,wcov);
        difftemp = y - x*betatemp;
        
        if (adjust){
          lossnew = accu(difftemp % difftemp)/(2*n)+ a*l*accu(abs(betatemp.subvec(ncov,p-1))) + (1-a)*l*accu(pow(betatemp.subvec(ncov,p-1),2))/2 + mu*pow(accu(betatemp.subvec(ncov,p-1)) + alpha,2)/2 + wcov*a*l*accu(abs(betatemp.subvec(0,ncov-1))) + wcov*(1-a)*l*accu(pow(betatemp.subvec(0,ncov-1),2))/2;
        }else{
          lossnew = accu(difftemp % difftemp)/(2*n)+ a*l*accu(abs(betatemp.subvec(ncov,p-1))) + (1-a)*l*accu(pow(betatemp.subvec(ncov,p-1),2))/2 + mu*pow(accu(betatemp.subvec(ncov,p-1)) + alpha,2)/2;
        }
      }
    }
    
    loss(i) = lossnew;
    mse(i) = accu(difftemp % difftemp)/n;
    beta.col(i) = betatemp;
    tol(i) = mean(abs(betatemp - betai));
    
  }
  
  Rcpp::List ret;
  ret["beta"] = beta;
  ret["lambda"] = lambda;
  ret["loss"] = loss;
  ret["mse"] = mse;
  ret["tol"] = tol;
  
  return ret;
  
}


// [[Rcpp::export]]
Rcpp::List logistic_enet_al(arma::mat x, arma::vec y, int len, double mu, int ub, arma::vec lambda, double wcov, double a, bool adjust, unsigned int ncov, bool display_progress=true, bool loop1=false, bool loop2=false){
  
  int n = y.n_elem;
  unsigned int p = x.n_cols;
  
  arma::vec beta0;
  beta0.zeros(len);
  arma::mat beta;
  beta.zeros(p,len);
  arma::vec loss = vec(len);
  arma::vec mse = vec(len);
  arma::vec tol = vec(len);
  
  Progress prog(len, display_progress);
  
  for (int i=0; i<len; ++i){
    
    prog.increment();
    
    double l = lambda(i);
    
    double beta0i = beta0(i);
    arma::vec betai = beta.col(i);
    if (i >= 1) {
      beta0i = beta0(i-1);
      betai = beta.col(i-1);
    }
    
    double alpha = 0;
    int k = 0;
    int k1 = 0;
    
    // Rcpp::Rcout << "i: " << i<< endl;
    
    arma::vec prob = exp(beta0i + x*betai)/(1 + exp(beta0i + x*betai));
    arma::vec sfun = y - prob;
    arma::vec hfun = prob % (1-prob) + 1e-8;
    arma::vec z0 = beta0i + x*betai + sfun/hfun;
    arma::vec z = z0 - mean(z0);
    arma::mat xx = x.t()*diagmat(hfun)*x;
    arma::vec xz = x.t()*(z % hfun);
    
    arma::vec diff = z - x*betai;
    
    double loss0 = 0;
    if (adjust){
      loss0 = accu(hfun % diff % diff)/(2*n)+ a*l*accu(abs(betai.subvec(ncov,p-1))) + (1-a)*l*accu(pow(betai.subvec(ncov,p-1),2))/2 + mu*pow(accu(betai.subvec(ncov,p-1)) + alpha,2)/2 + wcov*a*l*accu(abs(betai.subvec(0,ncov-1))) + wcov*(1-a)*l*accu(pow(betai.subvec(0,ncov-1),2))/2;
    }else{
      loss0 = accu(hfun % diff % diff)/(2*n)+ a*l*accu(abs(betai.subvec(ncov,p-1))) + (1-a)*l*accu(pow(betai.subvec(ncov,p-1),2))/2 + mu*pow(accu(betai.subvec(ncov,p-1)) + alpha,2)/2;    
    }
    
    arma::vec betatemp = gd_cov_al(xx,xz,n,l,a,betai,mu,alpha,adjust,ncov,wcov);
    arma::vec difftemp = z - x*betatemp;
    
    double lossnew = 0;
    if (adjust){
      lossnew = accu(hfun % difftemp % difftemp)/(2*n)+ a*l*accu(abs(betatemp.subvec(ncov,p-1))) + (1-a)*l*accu(pow(betatemp.subvec(ncov,p-1),2))/2 + mu*pow(accu(betatemp.subvec(ncov,p-1)) + alpha,2)/2 + wcov*a*l*accu(abs(betatemp.subvec(0,ncov-1))) + wcov*(1-a)*l*accu(pow(betatemp.subvec(0,ncov-1),2))/2;
    }else{
      lossnew = accu(hfun % difftemp % difftemp)/(2*n)+ a*l*accu(abs(betatemp.subvec(ncov,p-1))) + (1-a)*l*accu(pow(betatemp.subvec(ncov,p-1),2))/2 + mu*pow(accu(betatemp.subvec(ncov,p-1)) + alpha,2)/2;    
    }
    
    k1 = 0;
    
    while (abs(lossnew - loss0) > 1e-7 and k1 < ub){
      
      if (loop2){
        beta0i = mean(z0 - x*betatemp);
        prob = exp(beta0i + x*betatemp)/(1 + exp(beta0i + x*betatemp));
        sfun = y - prob;
        hfun = prob % (1-prob) + 1e-8;
        z0 = beta0i + x*betatemp + sfun/hfun;
        z = z0 - mean(z0);
        xx = x.t()*diagmat(hfun)*x;
        xz = x.t()*(z % hfun);
      }
      
      // Rcpp::Rcout << "lossnew: " << lossnew << endl;
      k1 = k1 + 1;
      loss0 = lossnew;
      betatemp = gd_cov_al(xx,xz,n,l,a,betatemp,mu,alpha,adjust,ncov,wcov);
      difftemp = z - x*betatemp;
      if (adjust){
        lossnew = accu(hfun % difftemp % difftemp)/(2*n)+ a*l*accu(abs(betatemp.subvec(ncov,p-1))) + (1-a)*l*accu(pow(betatemp.subvec(ncov,p-1),2))/2 + mu*pow(accu(betatemp.subvec(ncov,p-1)) + alpha,2)/2 + wcov*a*l*accu(abs(betatemp.subvec(0,ncov-1))) + wcov*(1-a)*l*accu(pow(betatemp.subvec(0,ncov-1),2))/2;
      }else{
        lossnew = accu(hfun % difftemp % difftemp)/(2*n)+ a*l*accu(abs(betatemp.subvec(ncov,p-1))) + (1-a)*l*accu(pow(betatemp.subvec(ncov,p-1),2))/2 + mu*pow(accu(betatemp.subvec(ncov,p-1)) + alpha,2)/2;    
      }
    }
    
    while (mean(abs(betatemp - betai)) > 1e-7 and k < ub){
      
      if (loop1){
        beta0i = mean(z0 - x*betatemp);
        prob = exp(beta0i + x*betatemp)/(1 + exp(beta0i + x*betatemp));
        sfun = y - prob;
        hfun = prob % (1-prob) + 1e-8;
        z0 = beta0i + x*betatemp + sfun/hfun;
        z = z0 - mean(z0);
        xx = x.t()*diagmat(hfun)*x;
        xz = x.t()*(z % hfun);
      }
      
      k = k+1;
      // Rcpp::Rcout << "k: " << k << endl;
      
      betai = betatemp;
      alpha = alpha + accu(betai.subvec(ncov,p-1));
      
      diff = z - x*betai;
      loss0 = accu(hfun % diff % diff)/(2*n) + l*accu(abs(betai)) + mu*pow(accu(betai) + alpha,2)/2;
      betatemp = gd_cov_al(xx,xz,n,l,a,betai,mu,alpha,adjust,ncov,wcov);
      difftemp = z - x*betatemp;
      if (adjust){
        lossnew = accu(hfun % difftemp % difftemp)/(2*n)+ a*l*accu(abs(betatemp.subvec(ncov,p-1))) + (1-a)*l*accu(pow(betatemp.subvec(ncov,p-1),2))/2 + mu*pow(accu(betatemp.subvec(ncov,p-1)) + alpha,2)/2 + wcov*a*l*accu(abs(betatemp.subvec(0,ncov-1))) + wcov*(1-a)*l*accu(pow(betatemp.subvec(0,ncov-1),2))/2;
      }else{
        lossnew = accu(hfun % difftemp % difftemp)/(2*n)+ a*l*accu(abs(betatemp.subvec(ncov,p-1))) + (1-a)*l*accu(pow(betatemp.subvec(ncov,p-1),2))/2 + mu*pow(accu(betatemp.subvec(ncov,p-1)) + alpha,2)/2;    
      }
      
      k1 = 0;
      
      while (abs(lossnew - loss0) > 1e-7 and k1 < ub){
        
        if (loop2){
          beta0i = mean(z0 - x*betatemp);
          prob = exp(beta0i + x*betatemp)/(1 + exp(beta0i + x*betatemp));
          sfun = y - prob;
          hfun = prob % (1-prob) + 1e-8;
          z0 = beta0i + x*betatemp + sfun/hfun;
          z = z0 - mean(z0);
          xx = x.t()*diagmat(hfun)*x;
          xz = x.t()*(z % hfun);
        }
        
        // Rcpp::Rcout << "lossnew: " << lossnew << endl;
        k1 = k1 + 1;
        loss0 = lossnew;
        betatemp = gd_cov_al(xx,xz,n,l,a,betatemp,mu,alpha,adjust,ncov,wcov);
        difftemp = z - x*betatemp;
        if (adjust){
          lossnew = accu(hfun % difftemp % difftemp)/(2*n)+ a*l*accu(abs(betatemp.subvec(ncov,p-1))) + (1-a)*l*accu(pow(betatemp.subvec(ncov,p-1),2))/2 + mu*pow(accu(betatemp.subvec(ncov,p-1)) + alpha,2)/2 + wcov*a*l*accu(abs(betatemp.subvec(0,ncov-1))) + wcov*(1-a)*l*accu(pow(betatemp.subvec(0,ncov-1),2))/2;
        }else{
          lossnew = accu(hfun % difftemp % difftemp)/(2*n)+ a*l*accu(abs(betatemp.subvec(ncov,p-1))) + (1-a)*l*accu(pow(betatemp.subvec(ncov,p-1),2))/2 + mu*pow(accu(betatemp.subvec(ncov,p-1)) + alpha,2)/2;    
        }
      }
      
      // Rcpp::Rcout << "k1: " << k1<< endl;
      
    }
    
    beta0(i) = mean(z0 - x*betatemp);
    loss(i) = lossnew;
    mse(i) = accu(difftemp % difftemp)/n;
    beta.col(i) = betatemp;
    tol(i) = mean(abs(betatemp - betai));
    
  }
  
  Rcpp::List ret;
  ret["beta"] = beta;
  ret["beta0"] = beta0;
  ret["lambda"] = lambda;
  ret["loss"] = loss;
  ret["mse"] = mse;
  ret["tol"] = tol;
  
  return ret;
  
}


// [[Rcpp::export]]
Rcpp::List cox_enet_al(arma::mat x, arma::vec t, arma::vec d, arma::vec tj, int len, double mu, int ub, arma::vec lambda, double wcov, double a, bool adjust, unsigned int ncov, double devnull, bool display_progress=true, bool loop1=false, bool loop2=false, bool notcv=true){
  
  int n = t.n_elem;
  unsigned int p = x.n_cols;
  int m = tj.n_elem;
  
  arma::mat beta;
  beta.zeros(p,len);
  arma::vec loss = vec(len);
  arma::vec mse = vec(len);
  arma::vec tol = vec(len);
  arma::vec loglik;
  loglik.zeros(len);
  
  Progress prog(len, display_progress);
  
  for (int i=0; i<len; ++i){
    
    prog.increment();
    
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
    
    if (i >= 1){
      // Rcpp::Rcout << "Delta Dev: " << -2*loglik(i-1) - devnull<< endl;
      // Rcpp::Rcout << "0.99Devnull: " << 0.99*devnull << endl;
      
      if (notcv){
        if (-2*loglik(i-1) - devnull >= 0.99*devnull){
          break;
        }
      }
      
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
    
    double loss0 = 0;
    if (adjust){
      loss0 = accu(hfun % diff % diff)/(2*n)+ a*l*accu(abs(betai.subvec(ncov,p-1))) + (1-a)*l*accu(pow(betai.subvec(ncov,p-1),2))/2 + mu*pow(accu(betai.subvec(ncov,p-1)) + alpha,2)/2 + wcov*a*l*accu(abs(betai.subvec(0,ncov-1))) + wcov*(1-a)*l*accu(pow(betai.subvec(0,ncov-1),2))/2;
    }else{
      loss0 = accu(hfun % diff % diff)/(2*n)+ a*l*accu(abs(betai.subvec(ncov,p-1))) + (1-a)*l*accu(pow(betai.subvec(ncov,p-1),2))/2 + mu*pow(accu(betai.subvec(ncov,p-1)) + alpha,2)/2;    
    }
    
    arma::vec betatemp = gd_cov_al(xx,xz,n,l,a,betai,mu,alpha,adjust,ncov,wcov);
    arma::vec difftemp = z - x*betatemp;
    
    double lossnew = 0;
    if (adjust){
      lossnew = accu(hfun % difftemp % difftemp)/(2*n)+ a*l*accu(abs(betatemp.subvec(ncov,p-1))) + (1-a)*l*accu(pow(betatemp.subvec(ncov,p-1),2))/2 + mu*pow(accu(betatemp.subvec(ncov,p-1)) + alpha,2)/2 + wcov*a*l*accu(abs(betatemp.subvec(0,ncov-1))) + wcov*(1-a)*l*accu(pow(betatemp.subvec(0,ncov-1),2))/2;
    }else{
      lossnew = accu(hfun % difftemp % difftemp)/(2*n)+ a*l*accu(abs(betatemp.subvec(ncov,p-1))) + (1-a)*l*accu(pow(betatemp.subvec(ncov,p-1),2))/2 + mu*pow(accu(betatemp.subvec(ncov,p-1)) + alpha,2)/2;    
    }
    
    while (abs(lossnew - loss0) > 1e-7 and k1 < ub){
      
      if (loop2){
        
        d1.zeros(n);
        d2.zeros(n);
        
        widx = find(t >= tj(0));
        denomj(0) = accu(link(widx)) + 1e-8;
        
        d1(widx) = d1(widx) + 1/denomj(0);
        d2(widx) = d2(widx) + 1/pow(denomj(0),2);
        
        for (int j=1; j<m; ++j){
          idx = find(t >= tj(j-1) && t < tj(j));
          denomj(j) = denomj(j-1) - accu(link(idx)) + 1e-8;
          
          widx = find(t >= tj(j));
          d1(widx) = d1(widx) + 1/denomj(j);
          d2(widx) = d2(widx) + 1/pow(denomj(j),2);
          
        }
        
        link = exp(x*betatemp);
        sfun = d - link % d1;
        hfun = link % d1 - d2 % pow(link,2) + 1e-8;
        z =  x*betatemp + sfun/hfun;
        xx = x.t()*diagmat(hfun)*x;
        xz = x.t()*(z % hfun);
      }
      
      k1 = k1 + 1;
      loss0 = lossnew;
      betatemp = gd_cov_al(xx,xz,n,l,a,betatemp,mu,alpha,adjust,ncov,wcov);
      difftemp = z - x*betatemp;
      if (adjust){
        lossnew = accu(hfun % difftemp % difftemp)/(2*n)+ a*l*accu(abs(betatemp.subvec(ncov,p-1))) + (1-a)*l*accu(pow(betatemp.subvec(ncov,p-1),2))/2 + mu*pow(accu(betatemp.subvec(ncov,p-1)) + alpha,2)/2 + wcov*a*l*accu(abs(betatemp.subvec(0,ncov-1))) + wcov*(1-a)*l*accu(pow(betatemp.subvec(0,ncov-1),2))/2;
      }else{
        lossnew = accu(hfun % difftemp % difftemp)/(2*n)+ a*l*accu(abs(betatemp.subvec(ncov,p-1))) + (1-a)*l*accu(pow(betatemp.subvec(ncov,p-1),2))/2 + mu*pow(accu(betatemp.subvec(ncov,p-1)) + alpha,2)/2;    
      }
    }
    
    while (mean(abs(betatemp - betai)) > 1e-7 and k < ub){
      
      if (loop1){
        
        d1.zeros(n);
        d2.zeros(n);
        
        widx = find(t >= tj(0));
        denomj(0) = accu(link(widx)) + 1e-8;
        
        d1(widx) = d1(widx) + 1/denomj(0);
        d2(widx) = d2(widx) + 1/pow(denomj(0),2);
        
        for (int j=1; j<m; ++j){
          idx = find(t >= tj(j-1) && t < tj(j));
          denomj(j) = denomj(j-1) - accu(link(idx)) + 1e-8;
          
          widx = find(t >= tj(j));
          d1(widx) = d1(widx) + 1/denomj(j);
          d2(widx) = d2(widx) + 1/pow(denomj(j),2);
          
        }
        
        link = exp(x*betatemp);
        sfun = d - link % d1;
        hfun = link % d1 - d2 % pow(link,2) + 1e-8;
        z =  x*betatemp + sfun/hfun;
        xx = x.t()*diagmat(hfun)*x;
        xz = x.t()*(z % hfun);
      }
      
      k = k+1;
      betai = betatemp;
      alpha = alpha + accu(betai.subvec(ncov,p-1));
      
      diff = z - x*betai;
      loss0 = accu(hfun % diff % diff)/(2*n) + l*accu(abs(betai)) + mu*pow(accu(betai) + alpha,2)/2;
      betatemp = gd_cov_al(xx,xz,n,l,a,betai,mu,alpha,adjust,ncov,wcov);
      difftemp = z - x*betatemp;
      if (adjust){
        lossnew = accu(hfun % difftemp % difftemp)/(2*n)+ a*l*accu(abs(betatemp.subvec(ncov,p-1))) + (1-a)*l*accu(pow(betatemp.subvec(ncov,p-1),2))/2 + mu*pow(accu(betatemp.subvec(ncov,p-1)) + alpha,2)/2 + wcov*a*l*accu(abs(betatemp.subvec(0,ncov-1))) + wcov*(1-a)*l*accu(pow(betatemp.subvec(0,ncov-1),2))/2;
      }else{
        lossnew = accu(hfun % difftemp % difftemp)/(2*n)+ a*l*accu(abs(betatemp.subvec(ncov,p-1))) + (1-a)*l*accu(pow(betatemp.subvec(ncov,p-1),2))/2 + mu*pow(accu(betatemp.subvec(ncov,p-1)) + alpha,2)/2;    
      }
      
      k1 = 0;
      
      while (abs(lossnew - loss0) > 1e-7 and k1 < ub){
        
        if (loop2){
          
          d1.zeros(n);
          d2.zeros(n);
          
          widx = find(t >= tj(0));
          denomj(0) = accu(link(widx)) + 1e-8;
          
          d1(widx) = d1(widx) + 1/denomj(0);
          d2(widx) = d2(widx) + 1/pow(denomj(0),2);
          
          for (int j=1; j<m; ++j){
            idx = find(t >= tj(j-1) && t < tj(j));
            denomj(j) = denomj(j-1) - accu(link(idx)) + 1e-8;
            
            widx = find(t >= tj(j));
            d1(widx) = d1(widx) + 1/denomj(j);
            d2(widx) = d2(widx) + 1/pow(denomj(j),2);
            
          }
          
          link = exp(x*betatemp);
          sfun = d - link % d1;
          hfun = link % d1 - d2 % pow(link,2) + 1e-8;
          z =  x*betatemp + sfun/hfun;
          xx = x.t()*diagmat(hfun)*x;
          xz = x.t()*(z % hfun);
        }
        
        // Rcpp::Rcout << "lossnew: " << lossnew << endl;
        k1 = k1 + 1;
        loss0 = lossnew;
        betatemp = gd_cov_al(xx,xz,n,l,a,betatemp,mu,alpha,adjust,ncov,wcov);
        difftemp = z - x*betatemp;
        
        if (adjust){
          lossnew = accu(hfun % difftemp % difftemp)/(2*n)+ a*l*accu(abs(betatemp.subvec(ncov,p-1))) + (1-a)*l*accu(pow(betatemp.subvec(ncov,p-1),2))/2 + mu*pow(accu(betatemp.subvec(ncov,p-1)) + alpha,2)/2 + wcov*a*l*accu(abs(betatemp.subvec(0,ncov-1))) + wcov*(1-a)*l*accu(pow(betatemp.subvec(0,ncov-1),2))/2;
        }else{
          lossnew = accu(hfun % difftemp % difftemp)/(2*n)+ a*l*accu(abs(betatemp.subvec(ncov,p-1))) + (1-a)*l*accu(pow(betatemp.subvec(ncov,p-1),2))/2 + mu*pow(accu(betatemp.subvec(ncov,p-1)) + alpha,2)/2;    
        }
        
      }
    }
    
    loss(i) = lossnew;
    mse(i) = accu(difftemp % difftemp)/n;
    beta.col(i) = betatemp;
    tol(i) = mean(abs(betatemp - betai));
    
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
  ret["tol"] = tol;
  
  return ret;
  
}


// [[Rcpp::export]]
Rcpp::List cox_timedep_enet_al(arma::mat x, arma::vec t0, arma::vec t1, arma::vec d, arma::vec tj, int len, double mu, int ub, arma::vec lambda, double wcov, double a, bool adjust, unsigned int ncov, double devnull, bool display_progress=true){
  
  int n = t0.n_elem;
  unsigned int p = x.n_cols;
  int m = tj.n_elem;
  
  arma::mat beta;
  beta.zeros(p,len);
  arma::vec loss = vec(len);
  arma::vec mse = vec(len);
  arma::vec tol = vec(len);
  arma::vec loglik;
  loglik.zeros(len);
  
  Progress prog(len, display_progress);
  
  for (int i=0; i<len; ++i){
    
    prog.increment();
    
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
    
    arma::uvec idx = find(t1 >= tj(0) and t0 < tj(0));
    denomj(0) = accu(link(idx)) + 1e-8;
    
    arma::uvec tidx;
    arma::uvec widx = find(t1 >= tj(0) and t0 < tj(0));
    arma::uvec uidx = find(t1 == tj(0));
    
    if (i >= 1) {
      loglik(i-1) = loglik(i-1) + accu(x.rows(uidx)*betai) - log(accu(link(widx))+1e-8);
    }
    d1(widx) = d1(widx) + 1/denomj(0);
    d2(widx) = d2(widx) + 1/pow(denomj(0),2);
    
    for (int j=1; j<m; ++j){
      
      widx = find(t1 >= tj(j) and t0 < tj(j));
      denomj(j) = accu(link(widx)) + 1e-8;
      
      if (i >= 1){
        uidx = find(t1 == tj(j));
        loglik(i-1) = loglik(i-1) + accu(x.rows(uidx)*betai) - log(denomj(j));
      }
      
      d1(widx) = d1(widx) + 1/denomj(j);
      d2(widx) = d2(widx) + 1/pow(denomj(j),2);
      
    }
    
    if (i >= 1){
      // Rcpp::Rcout << "Delta Dev: " << -2*loglik(i-1) - devnull<< endl;
      // Rcpp::Rcout << "0.99Devnull: " << 0.99*devnull << endl;
      
      if (-2*loglik(i-1) - devnull >= 0.99*devnull){
        break;
      }
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
    
    double loss0 = 0;
    if (adjust){
      loss0 = accu(hfun % diff % diff)/(2*n)+ a*l*accu(abs(betai.subvec(ncov,p-1))) + (1-a)*l*accu(pow(betai.subvec(ncov,p-1),2))/2 + mu*pow(accu(betai.subvec(ncov,p-1)) + alpha,2)/2 + wcov*a*l*accu(abs(betai.subvec(0,ncov-1))) + wcov*(1-a)*l*accu(pow(betai.subvec(0,ncov-1),2))/2;
    }else{
      loss0 = accu(hfun % diff % diff)/(2*n)+ a*l*accu(abs(betai.subvec(ncov,p-1))) + (1-a)*l*accu(pow(betai.subvec(ncov,p-1),2))/2 + mu*pow(accu(betai.subvec(ncov,p-1)) + alpha,2)/2;    
    }
    
    arma::vec betatemp = gd_cov_al(xx,xz,n,l,a,betai,mu,alpha,adjust,ncov,wcov);
    arma::vec difftemp = z - x*betatemp;
    double lossnew = 0;
    if (adjust){
      lossnew = accu(hfun % difftemp % difftemp)/(2*n)+ a*l*accu(abs(betatemp.subvec(ncov,p-1))) + (1-a)*l*accu(pow(betatemp.subvec(ncov,p-1),2))/2 + mu*pow(accu(betatemp.subvec(ncov,p-1)) + alpha,2)/2 + wcov*a*l*accu(abs(betatemp.subvec(0,ncov-1))) + wcov*(1-a)*l*accu(pow(betatemp.subvec(0,ncov-1),2))/2;
    }else{
      lossnew = accu(hfun % difftemp % difftemp)/(2*n)+ a*l*accu(abs(betatemp.subvec(ncov,p-1))) + (1-a)*l*accu(pow(betatemp.subvec(ncov,p-1),2))/2 + mu*pow(accu(betatemp.subvec(ncov,p-1)) + alpha,2)/2;    
    }
    
    while (abs(lossnew - loss0) > 1e-7 and k1 < ub){
      k1 = k1 + 1;
      loss0 = lossnew;
      betatemp = gd_cov_al(xx,xz,n,l,a,betatemp,mu,alpha,adjust,ncov,wcov);
      difftemp = z - x*betatemp;
      
      if (adjust){
        lossnew = accu(hfun % difftemp % difftemp)/(2*n)+ a*l*accu(abs(betatemp.subvec(ncov,p-1))) + (1-a)*l*accu(pow(betatemp.subvec(ncov,p-1),2))/2 + mu*pow(accu(betatemp.subvec(ncov,p-1)) + alpha,2)/2 + wcov*a*l*accu(abs(betatemp.subvec(0,ncov-1))) + wcov*(1-a)*l*accu(pow(betatemp.subvec(0,ncov-1),2))/2;
      }else{
        lossnew = accu(hfun % difftemp % difftemp)/(2*n)+ a*l*accu(abs(betatemp.subvec(ncov,p-1))) + (1-a)*l*accu(pow(betatemp.subvec(ncov,p-1),2))/2 + mu*pow(accu(betatemp.subvec(ncov,p-1)) + alpha,2)/2;    
      }
      
    }
    
    while (mean(abs(betatemp - betai)) > 1e-7 and k < ub){
      
      k = k+1;
      betai = betatemp;
      alpha = alpha + accu(betai.subvec(ncov,p-1));
      
      diff = z - x*betai;
      loss0 = accu(hfun % diff % diff)/(2*n) + l*accu(abs(betai)) + mu*pow(accu(betai) + alpha,2)/2;
      betatemp = gd_cov_al(xx,xz,n,l,a,betai,mu,alpha,adjust,ncov,wcov);
      difftemp = z - x*betatemp;
      
      if (adjust){
        lossnew = accu(hfun % difftemp % difftemp)/(2*n)+ a*l*accu(abs(betatemp.subvec(ncov,p-1))) + (1-a)*l*accu(pow(betatemp.subvec(ncov,p-1),2))/2 + mu*pow(accu(betatemp.subvec(ncov,p-1)) + alpha,2)/2 + wcov*a*l*accu(abs(betatemp.subvec(0,ncov-1))) + wcov*(1-a)*l*accu(pow(betatemp.subvec(0,ncov-1),2))/2;
      }else{
        lossnew = accu(hfun % difftemp % difftemp)/(2*n)+ a*l*accu(abs(betatemp.subvec(ncov,p-1))) + (1-a)*l*accu(pow(betatemp.subvec(ncov,p-1),2))/2 + mu*pow(accu(betatemp.subvec(ncov,p-1)) + alpha,2)/2;    
      }

      k1 = 0;
      
      while (abs(lossnew - loss0) > 1e-7 and k1 < ub){
        // Rcpp::Rcout << "lossnew: " << lossnew << endl;
        k1 = k1 + 1;
        loss0 = lossnew;
        betatemp = gd_cov_al(xx,xz,n,l,a,betatemp,mu,alpha,adjust,ncov,wcov);
        difftemp = z - x*betatemp;
        
        if (adjust){
          lossnew = accu(hfun % difftemp % difftemp)/(2*n)+ a*l*accu(abs(betatemp.subvec(ncov,p-1))) + (1-a)*l*accu(pow(betatemp.subvec(ncov,p-1),2))/2 + mu*pow(accu(betatemp.subvec(ncov,p-1)) + alpha,2)/2 + wcov*a*l*accu(abs(betatemp.subvec(0,ncov-1))) + wcov*(1-a)*l*accu(pow(betatemp.subvec(0,ncov-1),2))/2;
        }else{
          lossnew = accu(hfun % difftemp % difftemp)/(2*n)+ a*l*accu(abs(betatemp.subvec(ncov,p-1))) + (1-a)*l*accu(pow(betatemp.subvec(ncov,p-1),2))/2 + mu*pow(accu(betatemp.subvec(ncov,p-1)) + alpha,2)/2;    
        }
        
      }
    }
    
    loss(i) = lossnew;
    mse(i) = accu(difftemp % difftemp)/n;
    beta.col(i) = betatemp;
    tol(i) = mean(abs(betatemp - betai));
    
    if (i == len-1){
      
      link = exp(x*betatemp);
      
      for (int j=0; j<m; ++j){
        widx = find(t1 >= tj(j) and t0 < tj(j));
        uidx = find(t1 == tj(j));
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
  ret["tol"] = tol;
  
  return ret;
  
}


// [[Rcpp::export]]
Rcpp::List fg_enet_al(arma::mat x, arma::vec t0, arma::vec t1, arma::vec d, arma::vec tj, arma::vec w, int len, double mu, int ub, arma::vec lambda, double wcov, double a, bool adjust, unsigned int ncov, double devnull, bool display_progress=true){
  
  int n = t0.n_elem;
  unsigned int p = x.n_cols;
  int m = tj.n_elem;
  
  arma::mat beta;
  beta.zeros(p,len);
  arma::vec loss = vec(len);
  arma::vec mse = vec(len);
  arma::vec tol = vec(len);
  arma::vec loglik;
  loglik.zeros(len);
  
  Progress prog(len, display_progress);
  
  for (int i=0; i<len; ++i){
    
    prog.increment();
    
    double l = lambda(i);
    
    arma::vec betai = beta.col(i);
    if (i >= 1) {
      betai = beta.col(i-1);
    }
    
    double alpha = 0;
    int k = 0;
    int k1 = 0;
    
    // Rcpp::Rcout << "i: " << i<< endl;
    
    arma::vec link = w % exp(x*betai);
    arma::vec denomj = vec(m);
    arma::vec d1;
    d1.zeros(n);
    arma::vec d2;
    d2.zeros(n);
    
    arma::uvec idx = find(t1 >= tj(0) and t0 < tj(0));
    denomj(0) = accu(link(idx)) + 1e-8;
    
    arma::uvec tidx;
    arma::uvec widx = find(t1 >= tj(0) and t0 < tj(0));
    arma::uvec uidx = find(t1 == tj(0));
    
    if (i >= 1) {
      loglik(i-1) = loglik(i-1) + accu(w(uidx) * x.rows(uidx)*betai) - log(accu(link(widx))+1e-8);
    }
    d1(widx) = d1(widx) + 1/denomj(0);
    d2(widx) = d2(widx) + 1/pow(denomj(0),2);
    
    for (int j=1; j<m; ++j){
      
      widx = find(t1 >= tj(j) and t0 < tj(j));
      denomj(j) = accu(link(widx)) + 1e-8;
      
      if (i >= 1){
        uidx = find(t1 == tj(j));
        loglik(i-1) = loglik(i-1) + accu(w(uidx) * x.rows(uidx)*betai) - log(denomj(j));
      }
      
      d1(widx) = d1(widx) + 1/denomj(j);
      d2(widx) = d2(widx) + 1/pow(denomj(j),2);
      
    }
    
    if (i >= 1){
      // Rcpp::Rcout << "Delta Dev: " << -2*loglik(i-1) - devnull<< endl;
      // Rcpp::Rcout << "0.99Devnull: " << 0.99*devnull << endl;
      
      if (-2*loglik(i-1) - devnull >= 0.99*devnull){
        break;
      }
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
    double loss0 = 0;
    if (adjust){
      loss0 = accu(hfun % diff % diff)/(2*n)+ a*l*accu(abs(betai.subvec(ncov,p-1))) + (1-a)*l*accu(pow(betai.subvec(ncov,p-1),2))/2 + mu*pow(accu(betai.subvec(ncov,p-1)) + alpha,2)/2 + wcov*a*l*accu(abs(betai.subvec(0,ncov-1))) + wcov*(1-a)*l*accu(pow(betai.subvec(0,ncov-1),2))/2;
    }else{
      loss0 = accu(hfun % diff % diff)/(2*n)+ a*l*accu(abs(betai.subvec(ncov,p-1))) + (1-a)*l*accu(pow(betai.subvec(ncov,p-1),2))/2 + mu*pow(accu(betai.subvec(ncov,p-1)) + alpha,2)/2;    
    }
    
    arma::vec betatemp = gd_cov_al(xx,xz,n,l,a,betai,mu,alpha,adjust,ncov,wcov);
    arma::vec difftemp = z - x*betatemp;
    double lossnew = 0;
    
    if (adjust){
      lossnew = accu(hfun % difftemp % difftemp)/(2*n)+ a*l*accu(abs(betatemp.subvec(ncov,p-1))) + (1-a)*l*accu(pow(betatemp.subvec(ncov,p-1),2))/2 + mu*pow(accu(betatemp.subvec(ncov,p-1)) + alpha,2)/2 + wcov*a*l*accu(abs(betatemp.subvec(0,ncov-1))) + wcov*(1-a)*l*accu(pow(betatemp.subvec(0,ncov-1),2))/2;
    }else{
      lossnew = accu(hfun % difftemp % difftemp)/(2*n)+ a*l*accu(abs(betatemp.subvec(ncov,p-1))) + (1-a)*l*accu(pow(betatemp.subvec(ncov,p-1),2))/2 + mu*pow(accu(betatemp.subvec(ncov,p-1)) + alpha,2)/2;    
    }

    while (abs(lossnew - loss0) > 1e-7 and k1 < ub){
      k1 = k1 + 1;
      loss0 = lossnew;
      betatemp = gd_cov_al(xx,xz,n,l,a,betatemp,mu,alpha,adjust,ncov,wcov);
      difftemp = z - x*betatemp;
      
      if (adjust){
        lossnew = accu(hfun % difftemp % difftemp)/(2*n)+ a*l*accu(abs(betatemp.subvec(ncov,p-1))) + (1-a)*l*accu(pow(betatemp.subvec(ncov,p-1),2))/2 + mu*pow(accu(betatemp.subvec(ncov,p-1)) + alpha,2)/2 + wcov*a*l*accu(abs(betatemp.subvec(0,ncov-1))) + wcov*(1-a)*l*accu(pow(betatemp.subvec(0,ncov-1),2))/2;
      }else{
        lossnew = accu(hfun % difftemp % difftemp)/(2*n)+ a*l*accu(abs(betatemp.subvec(ncov,p-1))) + (1-a)*l*accu(pow(betatemp.subvec(ncov,p-1),2))/2 + mu*pow(accu(betatemp.subvec(ncov,p-1)) + alpha,2)/2;    
      }
      
    }
    
    while (mean(abs(betatemp - betai)) > 1e-7 and k < ub){
      
      k = k+1;
      betai = betatemp;
      alpha = alpha + accu(betai.subvec(ncov,p-1));
      
      diff = z - x*betai;
      loss0 = accu(hfun % diff % diff)/(2*n) + l*accu(abs(betai)) + mu*pow(accu(betai) + alpha,2)/2;
      betatemp = gd_cov_al(xx,xz,n,l,a,betai,mu,alpha,adjust,ncov,wcov);
      difftemp = z - x*betatemp;
      
      if (adjust){
        lossnew = accu(hfun % difftemp % difftemp)/(2*n)+ a*l*accu(abs(betatemp.subvec(ncov,p-1))) + (1-a)*l*accu(pow(betatemp.subvec(ncov,p-1),2))/2 + mu*pow(accu(betatemp.subvec(ncov,p-1)) + alpha,2)/2 + wcov*a*l*accu(abs(betatemp.subvec(0,ncov-1))) + wcov*(1-a)*l*accu(pow(betatemp.subvec(0,ncov-1),2))/2;
      }else{
        lossnew = accu(hfun % difftemp % difftemp)/(2*n)+ a*l*accu(abs(betatemp.subvec(ncov,p-1))) + (1-a)*l*accu(pow(betatemp.subvec(ncov,p-1),2))/2 + mu*pow(accu(betatemp.subvec(ncov,p-1)) + alpha,2)/2;    
      }

      k1 = 0;
      
      while (abs(lossnew - loss0) > 1e-7 and k1 < ub){
        // Rcpp::Rcout << "lossnew: " << lossnew << endl;
        k1 = k1 + 1;
        loss0 = lossnew;
        betatemp = gd_cov_al(xx,xz,n,l,a,betatemp,mu,alpha,adjust,ncov,wcov);
        difftemp = z - x*betatemp;
        
        if (adjust){
          lossnew = accu(hfun % difftemp % difftemp)/(2*n)+ a*l*accu(abs(betatemp.subvec(ncov,p-1))) + (1-a)*l*accu(pow(betatemp.subvec(ncov,p-1),2))/2 + mu*pow(accu(betatemp.subvec(ncov,p-1)) + alpha,2)/2 + wcov*a*l*accu(abs(betatemp.subvec(0,ncov-1))) + wcov*(1-a)*l*accu(pow(betatemp.subvec(0,ncov-1),2))/2;
        }else{
          lossnew = accu(hfun % difftemp % difftemp)/(2*n)+ a*l*accu(abs(betatemp.subvec(ncov,p-1))) + (1-a)*l*accu(pow(betatemp.subvec(ncov,p-1),2))/2 + mu*pow(accu(betatemp.subvec(ncov,p-1)) + alpha,2)/2;    
        }
      }
    }
    
    loss(i) = lossnew;
    mse(i) = accu(difftemp % difftemp)/n;
    beta.col(i) = betatemp;
    tol(i) = mean(abs(betatemp - betai));
    
    if (i == len-1){
      
      link = exp(x*betatemp);
      
      for (int j=0; j<m; ++j){
        widx = find(t1 >= tj(j) and t0 < tj(j));
        uidx = find(t1 == tj(j));
        loglik(i) = loglik(i) + accu(w(uidx) * x.rows(uidx)*betai) - log(accu(link(widx))+1e-8);
      }
      
    }
    
  }
  
  Rcpp::List ret;
  ret["beta"] = beta;
  ret["lambda"] = lambda;
  ret["loss"] = loss;
  ret["mse"] = mse;
  ret["loglik"] = loglik;
  ret["tol"] = tol;
   
  return ret;
  
}