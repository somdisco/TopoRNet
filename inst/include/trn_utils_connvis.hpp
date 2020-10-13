#ifndef TRN_UTILS_CONNVIS_HPP
#define TRN_UTILS_CONNVIS_HPP

#ifndef RcppArmadillo_H
#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
#endif

#ifndef RcppParallel_H
#include <RcppParallel.h>
// [[Rcpp::depends(RcppParallel)]]
#endif

// [[Rcpp::plugins(cpp11)]]


namespace TRN_UTILS_CONNVIS{

// Parallel worker to assign CONNvis local ranks 
struct CONNvis_lrank_worker : public RcppParallel::Worker {
  
  // Inputs 
  const arma::umat& EL; // edge list
  const arma::vec& ew; // edge weights, must have length = nrow(EL)

  // Internal variables 
  arma::uvec unq_EL1; // unique values in first column of edge list ("from" vertices)
  unsigned int nunq_EL1; // length(unq_EL1), controls parallel processing 
  
  // output container
  arma::uvec lrank;
  
  // Constructor 
  CONNvis_lrank_worker(const arma::umat& EL, const arma::vec& ew)
    : EL(EL), ew(ew)
  {
    // Make sure the nrows(EL) = ew
    if(EL.n_rows != ew.n_elem) Rcpp::stop("nrow(edge list) != length(edge weights)");
    
    // Compute unique "from" vertex list 
    unq_EL1 = arma::sort(arma::unique(EL.col(0))); 
    nunq_EL1 = unq_EL1.n_elem; 
    
    lrank.set_size(EL.n_rows); 
  }
  
  // Find the rank of a single row 
  void lrank_from_vertj(unsigned int j) {
    
    // Get rows of edge list where from vert = j 
    arma::uvec these_edges = arma::find(EL.col(0) == unq_EL1[j]); 
    
    // Strip out edge weights for these edges and sort them 
    arma::vec x = ew.elem(these_edges); 
    x = arma::unique(x); 
    x = arma::sort(x, "descend");
    
    for(unsigned int k=0; k<these_edges.n_elem; ++k) {
      lrank(these_edges[k]) = std::distance(x.begin(), std::find(x.begin(), x.end(), ew(these_edges[k]) ) ); 
      lrank(these_edges[k]) += 1; // + 1 b/c distance from first element is 0,1,2,3 instead of 1,2,3,4, etc
    }
    
  }
  
  
  // Parallel operator - find BMU of each row of X in parallel
  void operator()(std::size_t begin, std::size_t end) {
    for(unsigned int j = begin; j < end; j++) {
      lrank_from_vertj(j);
    }
  }
  
  // Parallel call method
  void calc_parallel() {
    RcppParallel::parallelFor(0, nunq_EL1, *this);
  }
  
  // Non-parallel call method
  void calc_serial() {
    // Find BMU of each row of X
    for(unsigned int j=0; j<nunq_EL1; ++j) {
      lrank_from_vertj(j);
    }
  }
};


inline arma::uvec CONNvis_grank(const arma::vec& CADJ, const arma::vec& lrank_means) {
  // Inputs: CADJ is vectorized CADJ values 
  //         lrank_means are vector of lrank means, sorted in descending order 
  
  // Strip out the local rank means, append infinity to beginning, -infinity to end
  std::vector<double> lrank_means_ = arma::conv_to<std::vector<double>>::from(lrank_means);
  lrank_means_.insert(lrank_means_.begin(), std::numeric_limits<double>::max());
  lrank_means_.push_back(std::numeric_limits<double>::lowest());
  
  // Assign global ranks
  arma::uvec grank(CADJ.n_elem); 
  unsigned int grank_counter = 1;
  for(unsigned int j=0; j<(lrank_means_.size()-1); ++j) {
    
    // Find values in between, if none skip
    arma::uvec in_cur_rank = arma::find(CADJ < lrank_means_[j] && CADJ >= lrank_means_[j+1]);
    if(in_cur_rank.n_elem == 0) continue;
    
    // Set them to the current global rank, and increment it
    grank.elem(in_cur_rank).fill(grank_counter);
    grank_counter++;
  }
  
  return grank; 
}




} // close namespace 
#endif
