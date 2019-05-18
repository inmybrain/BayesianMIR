// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
using namespace Rcpp;
using namespace arma;


mat mvrnormArma(int n, vec mu, mat sigma) {
  int ncols = sigma.n_cols;
  mat Y = randn(n, ncols);
  return repmat(mu, 1, n).t() + Y * chol(sigma);
}

double log_sum_exp(const vec& x) {
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

IntegerVector oneMultinom(vec probs) {
  int k = probs.size();
  IntegerVector ans(k);
  rmultinom(1, probs.begin(), k, ans.begin());
  return(ans);
}

IntegerVector FindEveryOne(IntegerVector delta){
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

// [[Rcpp::export]]
Rcpp::List Run1Gibbs_cpp(
    Rcpp::IntegerVector ninst,
    arma::colvec Y,
    arma::mat X1,
    
    arma::colvec hp_mu_coef, // mu such that beta ~ MVN(mu, xxx)
    double hp_a, // a : sigma_eps^2 ~ IG(a,b)
    double hp_b, // b : sigma_eps^2 ~ IG(a,b)
    double hp_g_coef, 
    double hp_sig2_y,
    arma::colvec hp_pi,
    arma::mat hp_inv_Sig_coef, // Sig^{-1} such that beta ~ MVN(mu, hp_sig2_y/g_beta * Sig)
    
    arma::colvec coef,
    double sig2_error,
    Rcpp::IntegerVector delta
){
  int p = X1.n_cols;
  int n = Y.n_rows;
  
  // Update beta
  int pos = 0;
  mat Z(n, p);
  for(int ii = 0; ii < n; ii++) {
    for(int jj = 0; jj < ninst[ii]; jj++){
      if(delta[pos + jj] == 1){
        Z.row(ii) = X1.row(pos + jj);
      }
    }
    pos = pos + ninst[ii];
  }
  mat inv_cov = inv(Z.t() * Z / sig2_error + hp_inv_Sig_coef * hp_g_coef / hp_sig2_y);
  colvec mu1 = inv_cov * (Z.t() * (Y) / sig2_error + hp_g_coef / hp_sig2_y * hp_inv_Sig_coef * hp_mu_coef);
  coef = mvrnormArma(1, mu1, inv_cov).t();
  
  // Update sig2_error
  sig2_error = 1.0 / R::rgamma(n / 2.0 + hp_a,
                               1.0 / (sum((Y - Z * coef) % (Y - Z * coef)) / 2.0 + hp_b));
  
  // Update delta
  pos = 0;
  for(int ii = 0; ii < n; ii++) {
    vec hp_pi_tilde = log(hp_pi.rows(pos, pos + ninst[ii] - 1)) - (Y[ii] - X1.rows(pos, pos + ninst[ii] - 1) * coef) % (Y[ii] - X1.rows(pos, pos + ninst[ii] - 1) * coef) / (2 * sig2_error);
    hp_pi_tilde = exp(hp_pi_tilde - log_sum_exp(hp_pi_tilde));
    delta[seq(pos, pos + ninst[ii] - 1)] = oneMultinom(hp_pi_tilde);
    pos = pos + ninst[ii];
  }
  return List::create(
    Named("coef") = coef,
    Named("sig2_error") = sig2_error,
    Named("delta") = delta                    
  );
}


// [[Rcpp::export]]
Rcpp::List BMIR_cpp(int ntotal,
              int nwarm,
              int nthin,
              
              Rcpp::IntegerVector ninst,
              arma::colvec Y,
              arma::mat X1,
              
              arma::colvec hp_mu_coef, // mu such that beta ~ MVN(mu, xxx)
              double hp_a, // a : sigma_eps^2 ~ IG(a,b)
              double hp_b, // b : sigma_eps^2 ~ IG(a,b)
              double hp_g_coef,
              double hp_sig2_y,
              arma::colvec hp_pi,
              arma::mat hp_inv_Sig_coef, // Sig^{-1} such that beta ~ MVN(mu, hp_sig2_y/g_beta * Sig)
              
              arma::colvec coef,
              double sig2_error,
              Rcpp::IntegerVector delta,
              
              bool return_delta = false
){
  int iter;
  double tick = 0.2;
  printf("=============================================================\n");
  printf("Start warming up %d MCMC samples!\n", nwarm);
  printf("Progress :");
  for(iter = 1; iter <= nwarm; iter++){
    if(any(iter == regspace(1, int(nwarm * tick), nwarm))){
      printf("%3.f%%...", 100 * double(iter) / double(nwarm));
    }
    List res_gibbs = Run1Gibbs_cpp(ninst,
                                   Y, X1,
                                   hp_mu_coef, hp_a, hp_b, hp_g_coef, hp_sig2_y, hp_pi, hp_inv_Sig_coef,
                                   coef, sig2_error, delta
    );
    coef = as<colvec>(res_gibbs["coef"]);
    sig2_error = as<double>(res_gibbs["sig2_error"]);
    delta = as<IntegerVector>(res_gibbs["delta"]);
  } // nwarm
  printf("\n");
  printf("Finish warming up!\n");
  printf("-------------------------------------------------------------\n");
  
  int niter = ntotal - nwarm;
  int nsave = 1 + floor((niter - 1.0) / (double) nthin);
  // mat res_coef(niter/nthin, coef.n_rows);
  // mat res_sig2_error(niter/nthin, 1);
  sp_mat res_delta(nsave, delta.length());
  mat res_coef(niter, coef.n_rows);
  mat res_sig2_error(niter, 1);
  // sp_mat res_delta(niter, delta.length());
  NumericVector pip(X1.n_rows);
  printf("Start extracting %d MCMC samples!\n", niter);
  printf("Progress :");
  for(iter = 1; iter <= niter; iter++){
    if(any(iter == regspace(1, int(niter * tick), niter))){
      printf("%3.f%%...", 100 * double(iter) / double(niter));
    }
    List res_gibbs = Run1Gibbs_cpp(ninst,
                                   Y, X1,
                                   hp_mu_coef, hp_a, hp_b, hp_g_coef, hp_sig2_y, hp_pi, hp_inv_Sig_coef,
                                   coef, sig2_error, delta
    );
    coef = as<colvec>(res_gibbs["coef"]);
    sig2_error = as<double>(res_gibbs["sig2_error"]);
    delta = as<IntegerVector>(res_gibbs["delta"]);
    if(any(iter == regspace(1, nthin, niter))){
      if(return_delta){
        res_delta.row((iter - 1)/nthin) = as<rowvec>(delta);
      }
      pip += as<NumericVector>(delta);
    }
    res_coef.row(iter - 1) = coef.t();
    res_sig2_error.row(iter - 1) = sig2_error;
  }
  pip = pip / nsave;
  printf("\n");
  printf("Finish MCMC sampling!\n");
  printf("=============================================================\n");
  
  return List::create(
    Named("coef") = res_coef,
    Named("sig2_error") = res_sig2_error,
    Named("pip") = pip,
    Named("delta") = res_delta
  );
}