#include <iostream>
#include <cmath>
#include <RcppArmadillo.h>
#include "diagnostics.h"

// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp ;
using namespace arma;

// Estimate the y hat 
mat y_hat_sar(const mat& X, const mat& betas, double rho, const sp_mat& W ){
  
  mat Xb = X*betas.t();
  int n = X.n_rows;
  sp_mat I_sp = speye<sp_mat>(n,n);
  
  //mat invA = spsolve(A,I_nn, "lapack");
  sp_mat invA = I_sp+rho*W+pow(rho,2)*(W*W)+pow(rho,3)*(W*W*W);
  
  mat yhat = invA*Xb;

  return yhat;
}

// Log likelihood function of SAR
double SAR_loglikelihood(const mat& X, const mat& y, double rho, 
                                const mat& betas, double sigma2e, 
                                const mat& detval, const sp_mat& W){
  
  int n = X.n_rows;
  
  //find log( det(A) )
  uvec indt = find(detval.col(0) > rho, 1, "first");
  double logdetA = detval(indt[0],1);
  
  sp_mat I_sp = speye<sp_mat>(n,n);
  sp_mat A = I_sp - rho * W;
    
  mat Ay = A*y;
  
  mat Xb = X*trans(betas);
  mat AymXbAymXb= trans(Ay-Xb)*(Ay-Xb);
  
  double log_lik_sar = (-n/2)*(log(2*datum::pi)+log( pow(sigma2e,2))) 
  + logdetA - AymXbAymXb(0,0)/(2*pow(sigma2e,2));
  
  return log_lik_sar;
}
  
// draw rho-values using griddy Gibbs and inversion
double SAR_draw_rho(const mat& detval,const mat& e0e0,const mat& eded,const mat& e0ed,double sig){
  
    double rho_draw = 0;
    int nrho = detval.n_rows;
    
    vec rho_exist = detval.col(0);
    vec log_detrho = detval.col(1);
    
    vec iota = ones<vec>(nrho);
    
    // Calculate Log-S(e(rho*)) given rho*
    vec S_rho = e0e0(0,0)*iota + pow(rho_exist, 2) *eded(0,0) - 2*rho_exist*e0ed(0,0);
    
    // Calculate the Log-density
    vec log_den = log_detrho - (1.0/(2.0*sig)) * S_rho; 
    double adj = log_den.max();
    vec t = rep(adj, nrho);
    log_den = log_den - t ;
    
    // the density
    vec den = exp(log_den);
    // the interval separating rho is h=0.001
    double h = 0.001;
    
   
    // Integration to calculate the normalized constant
    // using the  trapezoid rule
    double ISUM = h* ( sum(den)-den[0]/2.0-den[nrho-1]/2.0  );
    vec norm_den = (1.0/ISUM)*den;
    
    // cumulative density
    vec cumu_den = cumsum(norm_den);
    
    // Inverse sampling
    double rnd = runif(1,0,1)(0)*sum(norm_den);
    uvec ind = find(cumu_den <= rnd) ;
    int idraw = max(ind);
    
    if(idraw >= 0 && idraw < nrho) 
        rho_draw = rho_exist(idraw);
    
    return rho_draw;
  }

// [[Rcpp::export]]

List sar_cpp_arma( arma::mat X, arma::vec y, arma::sp_mat W, arma::mat detval, 
                   int burnin, int Nsim){
  
  //Starting the MCMC SET UP
  
  int n = X.n_rows;
  int p = X.n_cols;
    
  //Prior distribution specifications
  //For Betas
  vec M0 = zeros(p);
  
  mat T0 = 100.0 * eye<mat>(p,p);
  
  //completely non-informative priors
  int c0(0.01);
  int d0(0.01);
  
  //Store MCMC results
  //int Nsim = 10000;
  //int burnin = 5000;
  
  mat Betas = zeros(Nsim, p);
    
  vec sigma2e = zeros(Nsim);
  
  vec rho = zeros(Nsim);
  
  //initial values for model parameters
  sigma2e[0] = 2;
  rho[0] = 0.5;
  
  int ce( n/2 + c0 );
  
  //Fixed matrix manipulations during the MCMC loops
  
  mat XTX = trans(X) * X;
  mat invT0 = inv(T0);
  
  mat T0M0 = invT0 * M0;
  mat tX = trans(X);
  
  //some fixed values when updating rho
  mat beta0 = solve(X,y);
  
  mat e0 = y-X*beta0;
  mat e0e0 = trans(e0) * e0;
  
  vec Wy = W *y;
  mat betad = solve(X,Wy);
  mat ed = Wy - X * betad;
  
  mat eded = trans(ed) *ed;
  mat e0ed = trans(e0) *ed;
  
  // initialise A
  sp_mat I_sp = speye<sp_mat>(n,n);
  sp_mat A = I_sp - rho[0] * W;
 
  // MCMC updating  
  
  for(int i=1;i<Nsim;i++) {
    // Gibbs sampler for updating Betas
    mat VV = (1.0/sigma2e[i-1]) *XTX + invT0;
    // We use Cholesky decomption to inverse the covariance matrix
    mat vBetas = inv_sympd(VV) ;//chol2inv(chol(VV))
    
    // Define A=I-rho*W 
    A = I_sp - rho[i-1] * W;
    
    vec Ay = A*y;
    mat mBetas = vBetas * (tX * (1.0/sigma2e[i-1])*Ay+T0M0);
    
    // When the number of independent variables is large, drawing from multivariate 
    // distribution could be time-comsuming
    mat cholV = trans(chol(vBetas));
    // draw form p-dimensitional independent norm distribution
    mat betas = rnorm(p);
    betas = mBetas + cholV * betas;
    Betas.row(i) = trans(betas);
  
    //### Gibbs sampler for updating sigma2e
    mat Xb = X*betas;
    mat e = Ay - Xb;
    mat de = 0.5 * trans(e) * e + d0;
    
    sigma2e(i) = 1/Rf_rgamma(ce,1/de(0,0));
    rho(i) = SAR_draw_rho(detval, e0e0, eded, e0ed, sigma2e[i]);    
  }
  
  // evaluate diagnostics
  vec log_lik_samples = zeros(Nsim-burnin+1);
  
  for(int i=burnin;i<Nsim;i++)
  {
  log_lik_samples[i-burnin] = SAR_loglikelihood( X, y, rho[i], Betas.row(i), 
                                    sigma2e[i], detval, W );
  }
  double log_lik_mean_theta = SAR_loglikelihood( X, y, mean( rho.subvec(burnin-1,Nsim-1) ),
                              mean( Betas.submat( burnin-1,0,Nsim-1,p-1  ) ), 
                             mean( sigma2e.subvec(burnin-1,Nsim-1) ), detval, 
                             W = W);

  double dic, pd;
  diagnostic_dic_pd(log_lik_samples,log_lik_mean_theta, dic, pd);
  
  double r2 = diagnostic_Rsquared(y, y_hat_sar(X, mean( Betas.submat( burnin-1,0,Nsim-1,p-1  ) ), mean( rho.subvec(burnin-1,Nsim-1) ),W ));

  mat direct, indirect, total;
  diagnostic_impacts( mean( Betas.submat( burnin-1,0,Nsim-1,p-1  ) ), mean( rho.subvec(burnin-1,Nsim-1) ),W ,
            direct, indirect, total);

  return List ::create( Named("Mbetas")= mean( Betas.submat( burnin-1,0,Nsim-1,p-1  ) ), 
                             Named("SDbetas") = stddev( Betas.tail_rows( Nsim-burnin )  ),
                             Named("Mrho")= mean( rho.subvec(burnin-1,Nsim-1) ), 
                             Named("SDrho") = stddev( rho.subvec(burnin-1,Nsim-1) ),
                             Named("Msigma2e")= mean( sigma2e.subvec(burnin-1,Nsim-1) ), 
                             Named("SDsigma2e") = stddev( sigma2e.subvec(burnin-1,Nsim-1) ),
                             Named("DIC") = dic,
                             Named("pD") = pd,
                             Named("Log_Likelihood") = log_lik_mean_theta,
                             Named("R_Squared") = r2,
                             Named("impact_direct") = direct,
                             Named("impact_idirect") = indirect,
                             Named("impact_total") = total
                             ); 
    
}



