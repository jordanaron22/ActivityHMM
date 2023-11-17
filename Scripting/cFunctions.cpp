#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
NumericVector vectorEqBool(NumericVector vec, double lod) {
	NumericVector vecBool;
	for(NumericVector::iterator i = vec.begin(); i != vec.end(); ++i) {
		if (*i == lod) {
			vecBool.push_back(1);
		} else {
			vecBool.push_back(0);
		}
  }
  return vecBool;
}

// [[Rcpp::export]]
NumericVector logClassificationC(int current_state, NumericVector act_obs, double mu, double sig, double act_binom, double lod) {
	NumericVector temp;
	if (current_state == 0) {
		temp = Rcpp::dnorm( act_obs, mu, sig, true );
	} else {
		NumericVector vec_eq = vectorEqBool(act_obs, lod);
		temp = log(((1-act_binom) * (1 - vec_eq) * Rcpp::dnorm( act_obs, mu, sig, false ))+(act_binom*vec_eq));
	}
	return temp;
}

