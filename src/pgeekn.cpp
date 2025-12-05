// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; indent-tabs-mode: nil; -*-

// [[Rcpp::plugins(cpp17)]]
// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include <Rcpp.h>
using namespace arma;
using namespace std;
using namespace Rcpp;
// [[Rcpp::depends(RcppProgress)]]
#include <progress.hpp>
#include <progress_bar.hpp>

// Forward declaration of gee_cor from pgee.cpp
Rcpp::List gee_cor(double N, arma::vec nt, arma::vec y, arma::mat X, 
                   Rcpp::Function linkinv, Rcpp::Function variance, 
                   arma::vec beta_new, std::string corstr, double maxclsz, 
                   bool scalefix, double scalevalue);

// [[Rcpp::export]]
Rcpp::List gee_NR_kn(double N, // Number of subjects
                     arma::vec nt, // number of obs per subject
                     arma::vec y,
                     arma::mat X,
                     double nx, // Number of covariates (ncol(X))
                     Rcpp::Function linkinv,
                     Rcpp::Function mueta,
                     Rcpp::Function variance,
                     arma::vec beta_new,
                     arma::cube Rhat, // estimated working correlation matrix (Ehat from cor_gee.R)
                     double fihat, // estimated scale parameter (fi from cor_gee.R)
                     double lambda,
                     double a,
                     double alpha1, // for original features
                     double alpha2, // for knockoff features
                     double ncov,
                     double wcov,
                     double eps=1e-6,
                     double muu=1e6){
  
  // Calculate number of original compositional features
  unsigned int p1 = (nx - ncov) / 2; // Number of original compositional features
  
  arma::vec aindex=cumsum(nt);
  arma::vec index;
  index.zeros(N);
  index.subvec(1,N-1) = aindex.subvec(0,N-2);
  
  arma::vec eta = X * beta_new;
  arma::vec mu = Rcpp::as<arma::vec>(linkinv(eta));
  
  arma::mat E1;
  arma::vec E2;  // Changed from mat to vec to support subvec() method
  
  if (a <= 1){
    
    if (ncov == 0){
      
      E1 = diagmat(lambda*(a + (1-a)*abs(beta_new))/(abs(beta_new)+eps));
      E2.ones(nx);
      // Separate constraints for original and knockoff features
      double constraint1 = muu * (accu(beta_new.subvec(0,p1-1)) + alpha1);
      double constraint2 = muu * (accu(beta_new.subvec(p1,nx-1)) + alpha2);
      E2.subvec(0,p1-1).fill(constraint1);
      E2.subvec(p1,nx-1).fill(constraint2);
      
    }else if (ncov > 0){
      
      // This is E on Wang et al.(2012)
      
      arma::vec lw;
      lw.ones(size(beta_new));
      lw.subvec(0,ncov-1) = lw.subvec(0,ncov-1)*wcov;
      
      E1 = diagmat((lambda*lw % (a + (1-a)*abs(beta_new)))/(abs(beta_new)+eps));
      E2.ones(nx);
      arma::vec muuvec(size(beta_new));
      muuvec.fill(muu);
      muuvec.subvec(0,ncov-1).zeros();
      
      // Separate constraints for original and knockoff features
      arma::vec beta_orig = beta_new.subvec(ncov,ncov+p1-1);
      arma::vec beta_knock = beta_new.subvec(ncov+p1,nx-1);
      
      double constraint1 = muu * (accu(beta_orig) + alpha1);
      double constraint2 = muu * (accu(beta_knock) + alpha2);
      
      E2.subvec(0,ncov-1).zeros(); // Covariates have no constraint
      E2.subvec(ncov,ncov+p1-1).fill(constraint1);
      E2.subvec(ncov+p1,nx-1).fill(constraint2);
      
    }
    
  }else if (a > 0){
    
    if (ncov == 0){
      
      arma::vec beta_abs = abs(beta_new);
      arma::vec b1 = zeros(size(beta_new));
      arma::uvec pos1 = find(beta_abs > lambda);
      b1(pos1).ones();
      
      arma::vec b2 = zeros(size(beta_new));
      arma::uvec pos2 = find(beta_abs < (lambda*a));
      b2(pos2).ones();
      
      E1 = diagmat((lambda*(1-b1)+b1%((lambda*a)-beta_new)%b2/(a-1))/(abs(beta_new)+eps));
      E2.ones(nx);
      // Separate constraints for original and knockoff features
      double constraint1 = muu * (accu(beta_new.subvec(0,p1-1)) + alpha1);
      double constraint2 = muu * (accu(beta_new.subvec(p1,nx-1)) + alpha2);
      E2.subvec(0,p1-1).fill(constraint1);
      E2.subvec(p1,nx-1).fill(constraint2);
      
    }else if (ncov > 0){
      
      // This is E on Wang et al.(2012)
      
      arma::vec beta_abs = abs(beta_new);
      arma::vec b1 = zeros(size(beta_new));
      arma::uvec pos1 = find(beta_abs > lambda);
      b1(pos1).ones();
      
      arma::vec b2 = zeros(size(beta_new));
      arma::uvec pos2 = find(beta_abs < (lambda*a));
      b2(pos2).ones();
      
      arma::vec lw;
      lw.ones(size(beta_new));
      lw.subvec(0,ncov-1) = lw.subvec(0,ncov-1)*wcov;
      
      E1 = diagmat((lw % (lambda*(1-b1)+b1%((lambda*a)-beta_new)%b2/(a-1)))/(abs(beta_new)+eps));
      E2.ones(nx);
      arma::vec muuvec(size(beta_new));
      muuvec.fill(muu);
      muuvec.subvec(0,ncov-1).zeros();
      
      // Separate constraints for original and knockoff features
      arma::vec beta_orig = beta_new.subvec(ncov,ncov+p1-1);
      arma::vec beta_knock = beta_new.subvec(ncov+p1,nx-1);
      
      double constraint1 = muu * (accu(beta_orig) + alpha1);
      double constraint2 = muu * (accu(beta_knock) + alpha2);
      
      E2.subvec(0,ncov-1).zeros(); // Covariates have no constraint
      E2.subvec(ncov,ncov+p1-1).fill(constraint1);
      E2.subvec(ncov+p1,nx-1).fill(constraint2);
      
    }
    
  }
  
  
  arma::vec sum201; //<-matrix(0,nx,1)      //gradient:S
  sum201.zeros(nx);
  
  arma::mat sum301; //<-matrix(0,nx,nx)     //naive variance:H
  sum301.zeros(nx,nx);
  
  arma::mat sum401; //<-matrix(0,nx,nx)     //naive variance:H
  sum401.zeros(nx,nx);
  
  arma::vec varimu = Rcpp::as<arma::vec>(variance(mu));
  arma::vec mee = Rcpp::as<arma::vec>(mueta(eta));
  
  for (int i=0; i<N; ++i){
    
    arma::vec ym;
    ym.zeros(nt(i));
    
    arma::mat bigD;
    bigD.zeros(nt(i),nx);
    
    arma::mat bigA;
    bigA.zeros(nt(i),nt(i));
    
    for (int j=0; j<nt(i); ++j){
      
      ym(j) = y(j+index(i)) - mu(j+index(i));
      
      bigA(j,j) = varimu(j+index(i));
      
      for (int k=0; k<nx; ++k){
        
        bigD(j,k) = mee(j+index(i)) * X(j+index(i),k);
        
      }
    }
    
    //working covariance matrix
    
    arma::mat Rmat = Rhat.slice(i);
    arma::mat bigV = pow(bigA,0.5) * Rmat.submat(0,0,size(nt(i),nt(i))) * pow(bigA,0.5);
    arma::mat nugmat = mat(size(bigV));
    nugmat.diag() += 1e-6;
    
    ////This is S in Wang et al.(2012)
    arma::vec sum200 = bigD.t() * pinv(bigV+nugmat) * ym;      //this is like gradient
    sum201 = sum201 + sum200;
    
    ////This is H in Wang et al.(2012)
    arma::mat sum300 = bigD.t() * pinv(bigV+nugmat) * bigD;    //this is for information matrix.
    sum301 = sum301 + sum300;
    
  }
  
  arma::vec S = fihat*sum201;
  arma::mat H = fihat*sum301;
  
  Rcpp::List ret;
  ret["E"] = E1;
  ret["C"] = E2;
  ret["S"] = S;
  ret["H"] = H;
  
  return ret;
}

// [[Rcpp::export]]
Rcpp::List gee_fit_kn(arma::vec y,
                      arma::mat X,
                      arma::vec nt, // number of obs per subject
                      Rcpp::Function linkinv,
                      Rcpp::Function mueta,
                      Rcpp::Function variance,
                      std::string corstr,
                      arma::vec lambda,
                      double a,
                      double ncov,
                      double wcov,
                      double tol=1e-3,
                      double eps=1e-6,
                      double muu=1e6,
                      int maxiter1=100,
                      int maxiter2=10,
                      bool scalefix=false, // indicator of fixed scale parameter
                      double scalevalue=1, //Value of the scale parameter (if fixed)
                      bool display_progress=true
){
  
  double N = nt.n_elem;
  double nx = X.n_cols;
  double maxclsz = nt.max();
  
  // Calculate number of original compositional features
  unsigned int p1 = (nx - ncov) / 2; // Number of original compositional features
  
  double len = lambda.n_elem;
  arma::mat beta_mat;
  beta_mat.zeros(nx,len);
  
  Progress prog(len, display_progress);
  
  double diff; 
  arma::vec iters = vec(len);
  arma::vec diffs = vec(len);
  
  for (int i=0; i<len; ++i){
    
    prog.increment();
    
    double l = lambda(i);
    int k = 0;
    
    arma::vec betai = beta_mat.col(i);
    if (i >= 1) betai = beta_mat.col(i-1);
    arma::vec beta0;
    
    diff = 1; // Initialize diff > tol for each lambda to start the while loop
    
    double alpha1 = 0; // for original features
    double alpha2 = 0; // for knockoff features
    
    while (abs(diff) > tol and k < maxiter1){
      
      k = k + 1;
      beta0 = betai;
      int k1 = 0;
      double diff1 = 1;
      arma::vec beta00 = beta0;
      
      if (muu > 0){
        
        while (abs(diff1) > tol and k1 < maxiter2){
          
          k1 = k1 + 1;
          beta00 = betai;
          
          Rcpp::List cor_obj = gee_cor(N,
                                       nt,
                                       y,
                                       X,
                                       linkinv,
                                       variance,
                                       betai,
                                       corstr,
                                       maxclsz, 
                                       scalefix,
                                       scalevalue
          );
          
          arma::cube Rhat = cor_obj["Rhat"];
          double fihat = cor_obj["fihat"];
          
          Rcpp::List NR_obj=gee_NR_kn(N,
                                      nt,
                                      y,
                                      X,
                                      nx,
                                      linkinv,
                                      mueta,
                                      variance,
                                      betai,
                                      Rhat,
                                      fihat,
                                      l,
                                      a,
                                      alpha1,
                                      alpha2,
                                      ncov,
                                      wcov,
                                      eps,
                                      muu);
          
          arma::vec S = NR_obj["S"];
          arma::mat H = NR_obj["H"];
          arma::mat E = NR_obj["E"];
          arma::vec C = NR_obj["C"];
          
          arma::mat Nmat = mat(size(E),fill::value(N));
          arma::mat muumat = mat(size(E),fill::zeros); // this should be changed to assign muu only to the constrained part
          muumat.submat(ncov,ncov,nx-1,nx-1).fill(muu);
          arma::mat nugmat = mat(size(E));
          nugmat.diag() += 1e-10;
          
          betai = betai + pinv(H+Nmat%E+muumat+nugmat) * (S-((Nmat%E)*betai)-C);
          
          arma::vec diffvec1 = beta00-betai; 
          diff1 = max(abs(diffvec1));
          
        }
        
        // Update alpha1 and alpha2 separately
        if (ncov > 0){
          alpha1 = alpha1 + accu(betai.subvec(ncov,ncov+p1-1)); // for original features
          alpha2 = alpha2 + accu(betai.subvec(ncov+p1,nx-1)); // for knockoff features
        }else{
          alpha1 = alpha1 + accu(betai.subvec(0,p1-1)); // for original features
          alpha2 = alpha2 + accu(betai.subvec(p1,nx-1)); // for knockoff features
        }
        
      }else if (muu == 0){
        
        Rcpp::List cor_obj = gee_cor(N,
                                     nt,
                                     y,
                                     X,
                                     linkinv,
                                     variance,
                                     betai,
                                     corstr,
                                     maxclsz, 
                                     scalefix,
                                     scalevalue
        );
        
        arma::cube Rhat = cor_obj["Rhat"];
        double fihat = cor_obj["fihat"];
        
        Rcpp::List NR_obj=gee_NR_kn(N,
                                    nt,
                                    y,
                                    X,
                                    nx,
                                    linkinv,
                                    mueta,
                                    variance,
                                    betai,
                                    Rhat,
                                    fihat,
                                    l,
                                    a,
                                    alpha1,
                                    alpha2,
                                    ncov,
                                    wcov,
                                    eps,
                                    muu);
        
        arma::vec S = NR_obj["S"];
        arma::mat H = NR_obj["H"];
        arma::mat E = NR_obj["E"];
        arma::vec C = NR_obj["C"];
        
        arma::mat Nmat = mat(size(E),fill::value(N));
        arma::mat muumat = mat(size(E),fill::zeros); // this should be changed to assign muu only to the constrained part
        arma::mat nugmat = mat(size(E));
        nugmat.diag() += 1e-6;
        
        betai = betai + pinv(H+Nmat%E+muumat+nugmat) * (S-((Nmat%E)*betai)-C);
        
      }
      
      arma::vec diffvec = beta0-betai; 
      diff = max(abs(diffvec));
      
    }
    
    beta_mat.col(i) = betai;
    iters(i) = k;
    diffs(i) = diff;
    
  }
  
  Rcpp::List ret;
  
  ret["beta"] = beta_mat;
  ret["lambda"] = lambda;
  ret["iters"] = iters;
  ret["tol"] = diffs;
  
  return ret;
  
}

