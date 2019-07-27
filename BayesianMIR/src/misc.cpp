//------- Source from BMIR.v3.0.cpp: do not edit by hand
// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
using namespace Rcpp;
using namespace arma;


arma::mat mvrnormArma(int n, arma::vec mu, arma::mat sigma) {
  int ncols = sigma.n_cols;
  mat Y = randn(n, ncols);
  return repmat(mu, 1, n).t() + Y * chol(sigma);
}

double log_sum_exp(const arma::vec& x) {
  unsigned int maxi = x.index_max();
  double maxv = x(maxi);
  if (!(maxv > -datum::inf)) {
    return -datum::inf;
  }
  double cumsum = 0.0;
  for (unsigned int i = 0; i < x.n_elem; i++) {
    if ((i != maxi) & (x(i) > -datum::inf)) {
      cumsum += exp(x(i) - maxv);
    }
  }
  return maxv + log1p(cumsum);
}

Rcpp::IntegerVector oneMultinom(arma::vec probs) {
  int k = probs.size();
  IntegerVector ans(k);
  rmultinom(1, probs.begin(), k, ans.begin());
  return ans ;
}

Rcpp::IntegerVector FindEveryOne(Rcpp::IntegerVector delta){
  int n = delta.length();
  int jj = 0;
  IntegerVector idx;
  for(int ii = 0; ii < n; ii++){
    if(delta[ii] == 1){
      idx.push_back(ii);
      jj++;
    }
  }
  return idx;
}

