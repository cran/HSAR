#include <iostream>
#include <RcppArmadillo.h>
#include "diagnostics.h"
#include "gibbs_method.h"

// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
using namespace arma;
using namespace std;

// Estimate the y hat 
mat y_hat_hsar_rho_0(const mat& X, const mat& betas, const sp_mat& Z,const mat& Mus  ){
  
  mat Xb = X*betas.t();
  
  mat Zg = Z* Mus.t();
  
  mat yhat = Xb+Zg;

  return yhat;
}

// Log likelihood function of HSAR
double HSAR_loglikelihood_rho_0(const mat& X, const mat& y,  
                                const mat& betas, const mat& us, 
                                const vec& Unum, int Utotal, double sigma2e){
  
  int n = X.n_rows;
  
  mat Xb = X*trans(betas);
  
  mat Zu;
  mat nus = trans(us);
  for(int j=0;j<Utotal;j++){
    mat temp_add = repmat(nus.row(j),Unum[j], 1 );
      
    Zu.insert_rows( Zu.n_rows, temp_add);
  }
  
  mat crosprod_AymXbmZu= trans(y-Xb-Zu)*(y-Xb-Zu);
  
  double log_lik_sar = (-n/2)*(log(2*datum::pi)+log( pow(sigma2e,2))) 
  - crosprod_AymXbmZu(0,0)/(2*pow(sigma2e,2));
  
  return log_lik_sar;
}

// [[Rcpp::export]]

List hsar_cpp_arma_rho_0( arma::mat X, arma::vec y, arma::sp_mat M, 
                      arma::sp_mat Z, arma::mat detvalM, 
                      arma::vec Unum){
  
  //Starting the MCMC SET UP
  //arma_rng::set_seed(124);
  
  int n = X.n_rows;
  int p = X.n_cols;
  int Utotal = M.n_rows;
  
  //Prior distribution specifications
  //For Betas
  vec M0 = zeros(p);
  
  mat T0 = 100.0 * eye<mat>(p,p);
  
  //completely non-informative priors
  int c0(0.01);
  int d0(0.01);
  int a0(0.01);
  int b0(0.01);
  
  //Store MCMC results
  int Nsim = 10000;
  int burnin = 5000;
  
  mat Betas = zeros(Nsim, p);
  mat Us = zeros(Nsim, Utotal);
    
  vec sigma2e = zeros(Nsim);
  vec sigma2u = zeros(Nsim);
  
  vec lambda = zeros(Nsim);
  
  //initial values for model parameters
  sigma2e[0] = 2;
  sigma2u[0] = 2;
  
  lambda[0] = 0.5;
  
  int ce( n/2 + c0 );
  int au( Utotal/2 + a0 );
  
  //Fixed matrix manipulations during the MCMC loops
  
  mat XTX = trans(X) * X;
  mat invT0 = inv(T0);
  
  mat T0M0 = invT0 * M0;
  mat tX = trans(X);
  mat Zfull(Z);
  mat tZ = trans(Zfull);
 
 //  initialise B
  sp_mat I_sp_B = speye<sp_mat>(Utotal,Utotal);
  sp_mat B = I_sp_B - lambda[0] * M;

  // MCMC updating  
  
  for(int i=1;i<Nsim;i++) {
    // Gibbs sampler for updating Betas
    mat VV = (1.0/sigma2e[i-1]) *XTX + invT0;
    // We use Cholesky decomption to inverse the covariance matrix
    mat vBetas = inv_sympd(VV) ;//chol2inv(chol(VV))
    
    vec ZUs = Z*trans(Us.row(i-1)) ; //ZUs <- as.numeric(Z%*%Us[i-1,])
    mat mBetas = vBetas * (tX * (1.0/sigma2e[i-1])*(y-ZUs)+T0M0);    
    
    // When the number of independent variables is large, drawing from multivariate 
    // distribution could be time-comsuming
    mat cholV = trans(chol(vBetas));
    // draw form p-dimensitional independent norm distribution
    mat betas = Rcpp::rnorm(p);
    betas = mBetas + cholV * betas;
    Betas.row(i) = trans(betas);
  
    // Gibbs sampler for updating U. Now, us are spatially dependent.

    // Define and update B=I-lambda*M
    B = I_sp_B - lambda[i-1] * M;
    mat vU = tZ * (1.0/sigma2e[i-1])*Z+ trans(B) * (1.0/sigma2u[i-1])*B;
    vU = inv_sympd(vU); //vU <- chol2inv(chol(vU))

    vec Xb = X*betas;
    mat mU = vU * (tZ * (1.0/sigma2e[i-1])*(y-Xb ));
   
    // When the number of higher level units is large, drawing from multivariate 
    // distribution could be time-comsuming
    cholV = trans(chol(vU));

    // draw form J-dimensitional independent norm distribution
    mat us = Rcpp::rnorm(Utotal);
    us = mU + cholV * us;
    Us.row(i) = trans(us);  

    // Gibbs sampler for updating sigma2e
    mat Zu;
    for(int j=0;j<Utotal;j++){
      mat temp_add = repmat(us.row(j),Unum[j], 1 );
      
      Zu.insert_rows( Zu.n_rows, temp_add);
    }
    
    mat e = y - Zu -Xb;
    mat de = 0.5 * trans(e) * e + d0;
    
    sigma2e(i) = 1/Rf_rgamma(ce,1/de(0,0));
    
   // Gibbs sampler for updating sigma2u
   vec Bus = B*us;
   mat bu = 0.5 * trans(Bus) * Bus + b0;
   sigma2u(i) = 1/Rf_rgamma(au,1/bu(0,0));
    
    // Giddy Gibbs integration and inverse sampling for lambda
        
    mat uu = trans(us) *us;
    mat uMu = trans(us) * M * us;
    mat Mu = M * us;
    mat uMMu = trans(Mu) * Mu;
   
    lambda(i) = HSAR_draw_lambda(detvalM, uu, uMu, uMMu, sigma2u[i]);
    
  }
  
  // Diagnostics
  vec log_lik_samples = zeros(Nsim-burnin+1);
  
  for(int i=burnin;i<Nsim;i++)
  {
  log_lik_samples[i-burnin] = HSAR_loglikelihood_rho_0( X, y, Betas.row(i), 
                                    Us.row(i),Unum,Utotal,
                                    sigma2e[i] );
  }
  double log_lik_mean_theta = HSAR_loglikelihood_rho_0( X, y, 
                              mean( Betas.submat( burnin-1,0,Nsim-1,p-1  ) ), 
                              mean( Us.tail_rows( Nsim-burnin ) ), Unum,Utotal,
                             mean( sigma2e.subvec(burnin-1,Nsim-1) ));
                             
  double dic, pd;
  diagnostic_dic_pd(log_lik_samples,log_lik_mean_theta, dic, pd);
  
  double r2 = diagnostic_Rsquared(y, y_hat_hsar_rho_0(X, mean( Betas.submat( burnin-1,0,Nsim-1,p-1  ) ), Z ,mean( Us.tail_rows( Nsim-burnin ) ) ));
            
  return List ::create( Named("Mbetas")= mean( Betas.tail_rows( Nsim-burnin ) ), 
                             Named("SDbetas") = stddev( Betas.tail_rows( Nsim-burnin )  ),
                             Named("Mlambda")= mean( lambda.subvec(burnin-1,Nsim-1) ), 
                             Named("SDlambda") = stddev( lambda.subvec(burnin-1,Nsim-1) ),
                             Named("Msigma2e")= mean( sigma2e.subvec(burnin-1,Nsim-1) ), 
                             Named("SDsigma2e") = stddev( sigma2e.subvec(burnin-1,Nsim-1)),
                             Named("Msigma2u")= mean( sigma2u.subvec(burnin-1,Nsim-1) ), 
                             Named("SDsigma2u") = stddev( sigma2u.subvec(burnin-1,Nsim-1)),
                             Named("Mus")= mean( Us.tail_rows( Nsim-burnin ) ), 
                             Named("SDus") = stddev( Us.tail_rows( Nsim-burnin )  ),
                             Named("DIC") = dic,
                             Named("pD") = pd,
                             Named("Log_Likelihood") = log_lik_mean_theta,
                             Named("R_Squared") = r2
                             ); 
    
}


