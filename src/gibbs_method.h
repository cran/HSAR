#include <RcppArmadillo.h>

// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp ;
using namespace arma;

double HSAR_draw_rho(const mat& detval,const mat& e0e0,
                     const mat& eded, const mat& eueu,
                     const mat& e0ed,const mat& e0eu,
                     const mat& edeu, double sig);

double HSAR_draw_lambda(const mat& detvalM, const mat& uu,
                        const mat& uMu, const mat& uMMu, double sig);

