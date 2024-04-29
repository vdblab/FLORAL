// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; indent-tabs-mode: nil; -*-

// [[Rcpp::depends(ast2ast,RcppArmadillo)]]
#include <RcppArmadillo.h>
#include <Rcpp.h>
// [[Rcpp::plugins(cpp17)]]
using namespace arma;
using namespace std;
using namespace Rcpp;
// [[Rcpp::depends(RcppProgress)]]
#include <progress.hpp>
#include <progress_bar.hpp>
#include <etr.hpp>
using namespace etr;

// [[Rcpp::export]]
Rcpp::List gee_NR(double N, // Number of subjects
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
                  double eps=1e-6,
                  double muu=1e6){
  
  // Rcpp::Rcout << "Flag1" << endl;
  
  arma::vec aindex=cumsum(nt);
  arma::vec index;
  index.zeros(N);
  index.subvec(1,N-1) = aindex.subvec(0,N-2);
  
  // Rcpp::Rcout << "Flag2" << endl;
  
  arma::vec eta = X * beta_new;
  sexp linkinveta = linkinv(eta);
  arma::vec mu = linkinveta;
  
  // Rcpp::Rcout << "mu: " << mu << endl;
  
  // This is E on Wang et al.(2012)
  
  arma::mat E1 = diagmat(lambda/(abs(beta_new)+eps));
  arma::vec E2;
  E2.ones(nx);
  E2 = E2 * (muu * accu(beta_new));
  
  // if(is.null(pindex)==TRUE) {
  //   E<-E1 
  // } else 
  //   if(is.null(pindex)!=TRUE) {
  //     E1[,pindex]<-0
  //     E<-E1
  //   }
  //   
  arma::vec sum201; //<-matrix(0,nx,1)      //gradient:S
  sum201.zeros(nx);
  
  arma::mat sum301; //<-matrix(0,nx,nx)     //naive variance:H
  sum301.zeros(nx,nx);
  
  arma::mat sum401; //<-matrix(0,nx,nx)     //naive variance:H
  sum401.zeros(nx,nx);
  
  sexp variancemu = variance(mu);
  arma::vec varimu = variancemu;
  
  sexp muetaeta = mueta(eta);
  arma::vec mee = muetaeta;
  
  // Rcpp::Rcout << "Flag3" << endl;
  
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
    
    ////This is S in Wang et al.(2012)
    arma::vec sum200 = bigD.t() * pinv(bigV) * ym;      //this is like gradient
    sum201 = sum201 + sum200;
    
    ////This is H in Wang et al.(2012)
    arma::mat sum300 = bigD.t() * pinv(bigV) * bigD;    //this is for information matrix.
    sum301 = sum301 + sum300;
    
    //Speed up the code////
    arma::mat SSA = pow(pinv(bigA),2);
    arma::mat SRhat = pinv(Rmat.submat(0,0,size(nt(i),nt(i))));
    arma::vec SSAym = SSA * ym;
    
    arma::mat sum400 = bigD.t() * SSA * SRhat * (SSAym * SSAym.t()) * SRhat * SSA * bigD;
    sum401 = sum401 + sum400;
    
  }
  
  // Rcpp::Rcout << "Flag4" << endl;
  
  arma::mat S = fihat*sum201;
  arma::mat H = fihat*sum301;
  arma::mat M = fihat*sum401;
  
  Rcpp::List ret;
  // ret["eta"] = eta;
  // ret["mu"] = mu;
  ret["E"] = E1;
  ret["C"] = E2;
  ret["S"] = S;
  ret["H"] = H;
  ret["M"] = M;
  
  return ret;
}

// [[Rcpp::export]]
Rcpp::List gee_cor(double N, // Number of subjects
                   arma::vec nt, // number of obs per subject
                   arma::vec y,
                   arma::mat X,
                   Rcpp::Function linkinv,
                   Rcpp::Function variance,
                   arma::vec beta_new,
                   std::string corstr,
                   double maxclsz, // max number of obs
                   bool scalefix, // indicator of fixed scale parameter
                   double scalevalue=1 //Value of the scale parameter (if fixed)
){
  
  // Rcpp::Rcout << "Flag1" << endl;
  
  arma::vec eta = X * beta_new;
  sexp linkinveta = linkinv(eta);
  arma::vec mu = linkinveta;
  
  sexp variancemu = variance(mu);
  arma::vec varimu = variancemu;
  arma::vec sd = pow(varimu,0.5);
  arma::vec res = (y-mu)/sd;
  
  double fi;
  
  if(scalefix) {
    fi = accu(pow(res,2))/(accu(nt));
  }else{
    fi = scalevalue;
  }
  
  arma::vec aindex=cumsum(nt);
  arma::vec index;
  index.zeros(N);
  index.subvec(1,N-1) = aindex.subvec(0,N-2);
  
  double alfa_hat;
  arma::mat alfa_mat;
  
  if (corstr=="independence") {
    
    alfa_hat = 0;
    
  }else if (corstr == "exchangeable"){
    
    double sum1 = 0;
    double sum3 = 0;
    
    for (int i=0; i<N; ++i){
      for (int j=0; j<nt(i); ++j){
        for (int jj=0; jj<nt(i); ++jj){
          if (j != jj){
            double sum2 = res(j+index(i)) * res(jj+index(i));
            sum1 = sum1 + sum2;
          }
        }
      }
      double sum4 = nt(i)*(nt(i)-1);
      sum3 = sum3 + sum4;
    }
    alfa_hat = sum1 / (sum3*fi);
    
  }else if (corstr=="AR-1") { 
    double sum5 = 0;
    double sum6 = 0;
    
    for (int i=0; i<N; ++i){
      for (int j=0; j<nt(i); ++j){
        for (int jj=0; jj<nt(i); ++jj){
          if (j > jj && abs(j-jj)==1){
            double sum7 = res(j+index(i)) * res(jj+index(i));
            sum5 = sum5 + sum7;
          }
        }
      }
      double sum8 = (nt(i)-1);
      sum6 = sum6 + sum8;
    }
    alfa_hat = sum5 / (sum6*fi);
    
  }else if (corstr=="unstructured"){
    
    alfa_mat.zeros(nt(1),nt(1)); //not allowed for unequal number of cluster sizes.
    
    for (int j=0; j<nt(1); ++j){
      for (int jj=0; jj<nt(1); ++jj){
        
        double sum20 = 0;
        
        for (int i=0; i<N; ++i){
          double sum21 = res(j+index(i)) * res(jj+index(i));
          sum20 = sum21 + sum20;    
        }
        
        alfa_mat(j,jj) = sum20/(N*fi); 
        
      }
    }
  }
  
  arma::cube Ehat;
  Ehat.zeros(maxclsz,maxclsz,N); // For each subject, the corr matrix is maxclsz by maxclsz
  
  
  for (int i=0; i<N; ++i){
    
    arma::mat cor1;
    cor1.zeros(nt(i),nt(i));
    
    if (corstr=="independence") {
      
      cor1.eye(nt(i),nt(i));
      
    }else if (corstr == "exchangeable"){
      
      for (int j=0; j<nt(i); ++j){
        for (int jj=0; jj<nt(i); ++jj){
          if (j != jj){
            cor1(j,jj) = alfa_hat;
          }else{
            cor1(j,jj) = 1;
          }
        }
      }
      
    }else if (corstr=="AR-1") { 
      
      for (int j=0; j<nt(i); ++j){
        for (int jj=0; jj<nt(i); ++jj){
          cor1(j,jj) = pow(alfa_hat,abs(j-jj));
        }
      }
      
    }else if (corstr=="unstructured"){
      
      cor1 = alfa_mat + alfa_mat.t();
      cor1.diag().ones();
      
    }
    
    Ehat.slice(i).submat(0,0,nt(i)-1,nt(i)-1) = cor1; 
    
  }
  
  Rcpp::List ret;
  
  ret["Rhat"] = Ehat;
  ret["fihat"] = fi;
  
  return ret;
}