// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; indent-tabs-mode: nil; -*-

// [[Rcpp::plugins(cpp17)]]
// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include <Rcpp.h>
using namespace arma;
// [[Rcpp::depends(ast2ast)]]
#include <etr.hpp>
using namespace etr;
using namespace std;
using namespace Rcpp;
// [[Rcpp::depends(RcppProgress)]]
#include <progress.hpp>
#include <progress_bar.hpp>

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
                  double a,
                  double alpha,
                  double ncov,
                  double wcov,
                  double eps=1e-6,
                  double muu=1e6){
  
  // Rcpp::Rcout << "Flag1" << endl;
  
  arma::vec aindex=cumsum(nt);
  arma::vec index;
  index.zeros(N);
  index.subvec(1,N-1) = aindex.subvec(0,N-2);
  
  arma::vec eta = X * beta_new;
  sexp linkinveta = linkinv(eta);
  arma::vec mu = linkinveta;
  
  // Rcpp::Rcout << "mu: " << mu << endl;
  
  arma::mat E1;
  arma::mat E2;
  
  if (a <= 1){
    
    if (ncov == 0){
      
      E1 = diagmat(lambda*(a + (1-a)*abs(beta_new))/(abs(beta_new)+eps));
      E2.ones(nx);
      E2 = E2 * (muu * (accu(beta_new) + alpha));
      
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
      arma::vec beta_new_sub = beta_new.subvec(ncov,beta_new.n_elem-1);
      E2 = E2 % (muuvec * (accu(beta_new_sub) + alpha));
      
    }
    
  }else if (a > 0){
    
    if (ncov == 0){
      
      arma::vec beta_abs = abs(beta_new);
      arma::vec b1 = zeros(size(beta_new));
      arma::uvec pos1 = find(beta_abs > lambda);
      b1(pos1).ones();
      
      // Rcpp::Rcout << "b1: " << 1-b1 << endl;
      
      arma::vec b2 = zeros(size(beta_new));
      arma::uvec pos2 = find(beta_abs < (lambda*a));
      b2(pos2).ones();
      
      // Rcpp::Rcout << "b2: " << 1-b2 << endl;
      
      E1 = diagmat((lambda*(1-b1)+b1%((lambda*a)-beta_new)%b2/(a-1))/(abs(beta_new)+eps));
      E2.ones(nx);
      E2 = E2 * (muu * (accu(beta_new) + alpha));
      
    }else if (ncov > 0){
      
      // This is E on Wang et al.(2012)
      
      arma::vec beta_abs = abs(beta_new);
      arma::vec b1 = zeros(size(beta_new));
      arma::uvec pos1 = find(beta_abs > lambda);
      b1(pos1).ones();
      
      // Rcpp::Rcout << "b1: " << b1 << endl;
      
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
      arma::vec beta_new_sub = beta_new.subvec(ncov,beta_new.n_elem-1);
      E2 = E2 % (muuvec * (accu(beta_new_sub) + alpha));
      
    }
    
  }
  
  
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
    arma::mat nugmat = mat(size(bigV));
    nugmat.diag() += 1e-10;
    
    ////This is S in Wang et al.(2012)
    arma::vec sum200 = bigD.t() * pinv(bigV+nugmat) * ym;      //this is like gradient
    sum201 = sum201 + sum200;
    
    ////This is H in Wang et al.(2012)
    arma::mat sum300 = bigD.t() * pinv(bigV+nugmat) * bigD;    //this is for information matrix.
    sum301 = sum301 + sum300;
    
    //Speed up the code////
    arma::mat nugmat2 = mat(size(bigA));
    nugmat2.diag() += 1e-10;
    
    arma::mat SSA = pow(pinv(bigA+nugmat2),2);
    
    Rmat.submat(0,0,size(nt(i),nt(i))).diag() += 1e-10;
    arma::mat SRhat = pinv(Rmat.submat(0,0,size(nt(i),nt(i))));
    arma::vec SSAym = SSA * ym;
    
    arma::mat sum400 = bigD.t() * SSA * SRhat * (SSAym * SSAym.t()) * SRhat * SSA * bigD;
    sum401 = sum401 + sum400;
    
  }
  
  // Rcpp::Rcout << "Flag4" << endl;
  
  arma::vec S = fihat*sum201;
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

// [[Rcpp::export]]
Rcpp::List gee_fit(arma::vec y,
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
                   int maxiter=100,
                   bool scalefix=false, // indicator of fixed scale parameter
                   double scalevalue=1, //Value of the scale parameter (if fixed)
                   bool display_progress=true
){
  
  double N = nt.n_elem;
  double nx = X.n_cols;
  double maxclsz = nt.max();
  
  double len = lambda.n_elem;
  arma::mat beta_mat;
  beta_mat.zeros(nx,len);
  
  Progress prog(len, display_progress);
  
  double diff; 
  
  for (int i=0; i<len; ++i){
    
    prog.increment();
    
    double l = lambda(i);
    int k = 0;
    
    arma::vec betai = beta_mat.col(i);
    if (i >= 1) betai = beta_mat.col(i-1);
    arma::vec beta0;
    // double loss0;
    
    diff = 1; // Initialize diff > tol for each lambda to start the while loop
    
    // arma::vec eta = X * betai;
    // sexp linkinveta = linkinv(eta);
    // arma::vec mu = linkinveta;
    // arma::vec wt;
    // wt.ones(size(mu));
    // sexp resi = devresids(y,mu,wt);
    // arma::vec resid = resi;
    // double loss = accu(resid);
    
    double alpha = 0;
    
    while (abs(diff) > tol and k < maxiter){
      
      k = k + 1;
      beta0 = betai;
      // loss0 = loss;
      int k1 = 0;
      double diff1 = 1;
      arma::vec beta00 = beta0;
      
      
      if (muu > 0){
        
        while (abs(diff1) > tol and k1 < maxiter){
          
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
          
          Rcpp::List NR_obj=gee_NR(N,
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
                                   alpha,
                                   ncov,
                                   wcov,
                                   eps,
                                   muu);
          
          arma::vec S = NR_obj["S"];
          arma::mat H = NR_obj["H"];
          arma::mat E = NR_obj["E"];
          arma::vec C = NR_obj["C"];
          
          arma::mat Nmat = mat(size(E),fill::value(N));
          arma::mat muumat = mat(size(E),fill::value(muu));
          arma::mat nugmat = mat(size(E));
          nugmat.diag() += 1e-10;
          
          betai = betai + pinv(H+Nmat%E+muumat+nugmat) * (S-((Nmat%E)*betai)-C);
          
          // Rcpp::Rcout << "Flag440" << endl;
          
          arma::vec diffvec1 = beta00-betai; 
          diff1 = max(abs(diffvec1));
          
          // arma::vec eta1 = X * betai;
          // sexp linkinveta1 = linkinv(eta1);
          // arma::vec mu1 = linkinveta1;
          // arma::vec wt1;
          // wt1.ones(size(mu1));
          // sexp resi1 = devresids(y,mu1,wt1);
          // arma::vec resid1 = resi1;
          // loss = accu(resid1);
          
          // if (display_progress){
          // 
          //   Rcpp::Rcout << "Lambda:" << l << ",Iter:" << k1 << ",Diff:" << diff << endl;
          //   // Rcpp::Rcout << "Lambda:" << l << ",Iter:" << k1 << ",Diff loss:" << loss-loss0 << endl;
          // 
          // }
        }
        
        alpha = alpha + accu(betai.subvec(ncov,betai.n_elem-1));
        
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
        
        Rcpp::List NR_obj=gee_NR(N,
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
                                 alpha,
                                 ncov,
                                 wcov,
                                 eps,
                                 muu);
        
        arma::vec S = NR_obj["S"];
        arma::mat H = NR_obj["H"];
        arma::mat E = NR_obj["E"];
        arma::vec C = NR_obj["C"];
        
        arma::mat Nmat = mat(size(E),fill::value(N));
        arma::mat muumat = mat(size(E),fill::value(muu));
        arma::mat nugmat = mat(size(E));
        nugmat.diag() += 1e-10;
        
        betai = betai + pinv(H+Nmat%E+muumat+nugmat) * (S-((Nmat%E)*betai)-C);
        
      }
      
      arma::vec diffvec = beta0-betai; 
      diff = max(abs(diffvec));
      
    }
    
    beta_mat.col(i) = betai;
    
  }
  
  Rcpp::List ret;
  
  ret["beta"] = beta_mat;
  ret["lambda"] = lambda;
  
  return ret;
  
}