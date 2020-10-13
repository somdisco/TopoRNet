#include "trn_utils_dm_prune.hpp"

//' Density of Dirichlet-Multinomial Distribution
//' 
//' @description Computes the density of a Dirichlet-Multinomial (DM) distribution, at a given 
//' count vector \code{x}, with a user-specified prior count vector \code{prior}. 
//' 
//' @param x the count vector at which to evaluate the DM density 
//' @param prior the prior count vector, must have \code{length = length(x)}
//' @param log_form whether to return the log-density, default = TRUE
//' @param normalize whether to normalize the calculated density by \code{sum(x)}. 
//' Default = TRUE. 
//' 
//' @return The evaluated DM density (a number)
//' 
//' @export
//[[Rcpp::export]]
double dDirMult(const arma::vec& x, const arma::vec& prior, bool log_form = true, bool normalize = true) {
  return TRN_UTILS_DMPRUNE::dDirMult(x, prior, log_form, normalize);
}