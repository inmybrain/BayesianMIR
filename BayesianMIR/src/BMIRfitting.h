//------- Source from BMIR.v3.0.cpp: do not edit by hand
#ifndef BMIRfitting
#define BMIRfitting

Rcpp::List Run1Gibbs_cpp(
    Rcpp::IntegerVector ninst,
    arma::colvec Y,
    arma::mat X1,
    
    arma::colvec hp_mu_coef, // mu such that beta ~ MVN(mu, xxx)
;
Rcpp::List BMIR_cpp(int ntotal,
                    int nwarm,
                    int nthin,
                    
                    Rcpp::IntegerVector ninst,
                    arma::colvec Y,
                    arma::mat X1,
                    
                    arma::colvec hp_mu_coef, // mu such that beta ~ MVN(mu, xxx)
;
#endif
