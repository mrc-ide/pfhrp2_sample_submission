#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
SEXP loglike(Rcpp::NumericVector params, Rcpp::List data, Rcpp::List misc) {
  
  // extract parameters
  double A = params["A"];
  double B = params["B"];
  double K = params["K"];
  double w = params["w"];
  
  // unpack data
  std::vector<int> N = Rcpp::as< std::vector<int> >(data["N"]);
  std::vector<int> num = Rcpp::as< std::vector<int> >(data["num"]);
  std::vector<int> denom = Rcpp::as< std::vector<int> >(data["denom"]);
  
  // calculate loglikelihood
  double ret = 0.0;
  for (int i = 0; i < N.size(); ++i) {
    
    // predict power
    double power = A + (1 - A)*B / (1 + pow(K / N[i], w));
    
    // binomial likelihood
    ret += R::dbinom(num[i], denom[i], power, true);
  }
  
  return Rcpp::wrap(ret);
}

// [[Rcpp::export]]
SEXP logprior(Rcpp::NumericVector params, Rcpp::List misc) {
  return Rcpp::wrap(0.0);
}


// NOTE: Do not edit this function name
// [[Rcpp::export]]  
SEXP create_xptr(std::string function_name) {  
  typedef SEXP (*funcPtr_likelihood)(Rcpp::NumericVector params, Rcpp::List data, Rcpp::List misc);  
  typedef SEXP (*funcPtr_prior)(Rcpp::NumericVector params, Rcpp::List misc);  
  
  // NOTE: If your loglikelihood function is not called "loglike" please edit:
  if (function_name == "loglike"){
    return(Rcpp::XPtr<funcPtr_likelihood>(new funcPtr_likelihood(&loglike)));
  } 
  // NOTE: If your logprior function is not called "logprior" please edit:
  if (function_name == "logprior"){
    return(Rcpp::XPtr<funcPtr_prior>(new funcPtr_prior(&logprior)));
  } 
  
  stop("cpp function %i not found", function_name);
}
