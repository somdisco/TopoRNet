#ifndef TRN_UTILS_TOPOPRES_HPP
#include "trn_utils_topo_pres.hpp"
#endif

//' Topographic Product of a TRN
//' 
//' @description The Topographic Product of a TRN which is imbued with an output topology.  
//' 
//' @param input_coords the coordinates of the vertices in the Input space (i.e, R^d) of the TRN
//' @param output_coords the coordinates of the vertices in the Output space of the TRN 
//' @param parallel, whether to compute in parallel. Default = T.
//' 
//' @return The Topographic Product (a number)
//' 
//' @references 
//' \insertRef{BauerPawelzikGeisel1992}{TopoRNet}
//' @export
// [[Rcpp::export]]
double Topographic_Product(const arma::mat& input_coords, const arma::mat& output_coords, bool parallel = true) {
  return TRN_UTILS_TOPOPRES::Topographic_Product(input_coords, output_coords, parallel); 
}


//' Topographic Functions of a TRN 
//' 
//' @description The Topograph Functions of a TRN measure topology preservation across a range of "folding" lengths, 
//' which are defined as the length (geodesic distance) by which edges in one space (input / output) must "fold" in order 
//' to be represented in the other.  
//' 
//' The functions computed here are the TF (Topographic Function) of Villmann et al. and the 
//' DTF (Differential Topographic Function ) & WDTF (Weighted DTF) of Zhang et al. 
//' 
//' @param inputADJ the adjacency matrix of vertices in the Input space of the TRN 
//' @param outputADJ the adjacency matrix of vertices in the Output sapce of the TRN 
//' @param parallel, whether to compute in parallel. Default = T. 
//' 
//' @return a data frame with columns: 
//' \itemize{
//' \item \code{k} list of folding lengths
//' \item \code{TF} the TF computed at each folding length 
//' \item \code{DTF} the DTF computed at each folding length 
//' \item \code{WDTF} the WDTF computed at each folding length 
//' }
//' 
//' \insertRef{VillmannDerHerrmannMartinetz1997}{TopoRNet}
//' \insertRef{ZhangMerenyi2006}{TopoRNet}
//' @export
// [[Rcpp::export]]
Rcpp::DataFrame Topographic_Functions(const arma::umat& inputADJ, const arma::umat& outputADJ, bool parallel = true) {
  return TRN_UTILS_TOPOPRES::Topographic_Functions(inputADJ, outputADJ, parallel); 
}
