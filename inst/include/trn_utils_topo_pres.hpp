#ifndef RcppArmadillo_H
#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
#endif

#ifndef RcppParallel_H
#include <RcppParallel.h>
// [[Rcpp::depends(RcppParallel)]]
#endif

// [[Rcpp::plugins(cpp11)]]

#ifndef SOM_UTILS_DIST_HPP
#include "som_utils_dist.hpp"
#endif


#ifndef TRN_UTILS_TOPOPRES_HPP
#define TRN_UTILS_TOPOPRES_HPP

namespace TRN_UTILS_TOPOPRES {

struct topo_product_worker : public RcppParallel::Worker {
  
  // Inputs
  const arma::mat& W; // prototypes 
  const arma::mat& nu_xy; // vertex (x,y) coordinates 
  

  // Internal variables 
  unsigned int nW; 
  
  // Output 
  arma::vec sumlogP3; // stores sum_{k=1}^{nW-1} log(P3) calculated from every prototype


  // Constructor
  topo_product_worker(const arma::mat& W, const arma::mat& nu_xy)
    : W(W), nu_xy(nu_xy)
  {
    
    // Check that W and nu_xy have same # of rows 
    if(W.n_rows != nu_xy.n_rows) Rcpp::stop("nrow(W) != nrow(nu_xy)"); 
    if(nu_xy.n_cols != 2) Rcpp::stop("ncol(nu_xy) != 2");
    
    // Setup output container 
    nW = W.n_rows; 
    sumlogP3.set_size(nW); 
    sumlogP3.zeros(); 
  }
  

  // Find TP contribution from single prototype 
  void topo_product_from_j(unsigned int j) {
    // Subscripts "M" denote manifold (Rd) distances, "L" denote lattice distances 
    
    // Compute Euclidean dist from this proto to all others, and this neuron to all others 
    arma::vec dM(nW); dM.zeros(); 
    arma::vec dL(nW); dL.zeros(); 
    for(unsigned int k=0; k<nW; ++k) {
      if(k == j) continue; 
      // Manifold distance
      dM[k] = SOM_UTILS_DIST::dist_L22(W.row(j), W.row(k)); 
      // Latttice distance 
      dL[k] = SOM_UTILS_DIST::dist_L22(nu_xy.row(j), nu_xy.row(k)); 
    }
    
    // Find the index of the nearest neighbors of j, in both M and L 
    arma::uvec nnM = arma::sort_index(dM, "ascend");
    arma::uvec nnL = arma::sort_index(dL, "ascend"); 
    
    // Loop over the distances to the other N-1 protos, 
    // compute Q1 and Q2
    arma::vec logQ1(nW-1), logQ2(nW-1); 
    for(unsigned int k=1; k<nW; ++k) {
      logQ1[k-1] = std::log(dM[nnL[k]] / dM[nnM[k]]); 
      logQ2[k-1] = std::log(dL[nnL[k]] / dL[nnM[k]]);
    }
    
    // Compute logP3
    arma::vec logP3 = arma::cumsum(logQ1 + logQ2); 
    for(unsigned int k=1; k<nW; ++k) {
      logP3[k-1] = logP3[k-1] / (2.0*double(k)); 
    }
    
    // Store its sum 
    sumlogP3[j] = arma::accu(logP3); 
    
    
    return;
  }
  
  
  // Parallel operator - find BMU of each row of X in parallel
  void operator()(std::size_t begin, std::size_t end) {
    for(unsigned int j = begin; j < end; j++) {
      topo_product_from_j(j);
    }
  }
  
  // Parallel call method
  void calc_parallel() {
    RcppParallel::parallelFor(0, nW, *this);
  }
  
  // Non-parallel call method
  void calc_serial() {
    for(unsigned int j=0; j<nW; ++j) {
      topo_product_from_j(j);
    }
  }
  
  // Return the Topographic Product 
  double TP() {
    return arma::accu(sumlogP3) / (double(nW) * double(nW-1)); 
  }
};

inline double Topographic_Product(const arma::mat& input_coords, const arma::mat& output_coords, bool parallel) {
  
  topo_product_worker wkr(input_coords, output_coords); 
  if(parallel) wkr.calc_parallel(); else wkr.calc_serial(); 
  return wkr.TP(); 
  
} 


// *** Topographic Function *** 
struct topo_function_worker : public RcppParallel::Worker {
  
  // Inputs
  const arma::umat& inputADJ; 
  const arma::umat& outputADJ; 
  
  // Internal variables 
  unsigned int nW; 
  unsigned int maxK;
  arma::umat inGDIST; 
  arma::umat outGDIST; 
  
  // Output 
  arma::ivec K; 
  arma::vec TF; 
  arma::vec DTF; 
  arma::vec WDTF;
  
  // Constructor
  topo_function_worker(const arma::umat& inputADJ, const arma::umat& outputADJ)
    : inputADJ(inputADJ), outputADJ(outputADJ)
  {
    // Check that in & out topoADJ have same size 
    if(inputADJ.n_rows != outputADJ.n_rows) Rcpp::stop("nrow(inputADJ) != nrow(outputADJ)");
    if(inputADJ.n_cols != outputADJ.n_cols) Rcpp::stop("ncol(inputADJ) != ncol(outputADJ)");
    
    nW = inputADJ.n_rows; 
    
    // Compute the geodesic distances of the graphs
    inGDIST = SOM_UTILS_DIST::geodesicdist(inputADJ, false, false); // (adjacency, is_weighted, is_directed)
    outGDIST = SOM_UTILS_DIST::geodesicdist(outputADJ, false, false);

    maxK = std::max(inGDIST.max(), outGDIST.max());

    // Setup output containers
    K = arma::regspace<arma::ivec>(-maxK, maxK); 
    TF.set_size(K.n_elem); TF.zeros(); 
    DTF.set_size(K.n_elem); DTF.zeros(); 
    WDTF.set_size(K.n_elem); WDTF.zeros(); 
  }
  
  
  // Find TP contribution from single prototype 
  // Only consider the positive folding range here 
  void dtf_length_k(unsigned int k) {
    // Subscripts "M" denote manifold (Rd), "L" denote lattice 
    
    if(k==1) return; 
    
    arma::uvec find_posK = arma::find(K == k);
    arma::uvec find_negK = arma::find(K == -k);
    
    // Forward DTF
    arma::uvec fwd_fold = arma::find(inGDIST==1 && outGDIST==k);
    DTF[find_posK[0]] = double(fwd_fold.n_elem) / double(nW);
      
    // Forward WDTF
    WDTF[find_posK[0]] = double(arma::accu(inputADJ(fwd_fold))) / double(arma::accu(inputADJ)); 

    // Backward DTF
    arma::uvec bwd_fold = arma::find(inGDIST==k && outGDIST==1);
    DTF[find_negK[0]] = double(bwd_fold.n_elem) / double(nW);
    
    return;
  }
  
  // Parallel operator - find BMU of each row of X in parallel
  void operator()(std::size_t begin, std::size_t end) {
    for(unsigned int j = begin; j < end; j++) {
      dtf_length_k(j);
    }
  }

  // Set the TF values, must be done after calculation of DTF 
  void set_TF() {
    
    // Find where K=0, will use this to separate the postive & negative parts of DTF for summing 
    arma::uvec find_zeroK = arma::find(K == 0);
    unsigned int idx_zeroK = find_zeroK[0];
    
    // Extract the positive & negative folding lenghts 
    arma::vec posDTF = DTF.tail(maxK-1);
    posDTF = arma::reverse(posDTF);
    arma::vec negDTF = DTF.head(maxK-1); 
    
    arma::vec cumsum_posDTF = arma::cumsum(posDTF); 
    arma::vec cumsum_negDTF = arma::cumsum(negDTF); 
    
    TF.subvec((idx_zeroK+1), arma::size(cumsum_posDTF)) = arma::reverse(cumsum_posDTF);
    TF.subvec(1, arma::size(cumsum_negDTF)) = cumsum_negDTF;
    
    // Set the K=0 value
    TF[idx_zeroK] = TF[idx_zeroK+1] + TF[idx_zeroK-1];
  }
  
  // Remove empty folding lengths from all containers 
  void process_TFs() {
    std::vector<unsigned int> keep_these; 
    for(unsigned int k=0; k<K.n_elem; ++k) {
      if(TF[k]>0 || DTF[k]>0 || WDTF[k]>0) {
        keep_these.push_back(k); 
      }
    }
    
    arma::uvec keep_these_uvec = arma::conv_to<arma::uvec>::from(keep_these); 
    this->K = this->K.elem(keep_these_uvec); 
    this->TF = this->TF.elem(keep_these_uvec);
    this->DTF = this->DTF.elem(keep_these_uvec); 
    this->WDTF = this->WDTF.elem(keep_these_uvec); 
  }
  
  // Parallel call method
  void calc_parallel() {
    RcppParallel::parallelFor(2, maxK+1, *this);
    this->set_TF(); 
    this->process_TFs(); 
  }
  
  // Non-parallel call method
  void calc_serial() {
    for(unsigned int j=2; j<=maxK; ++j) {
      dtf_length_k(j);
    }
    this->set_TF(); 
    this->process_TFs(); 
  }
  
  // Return all TFs in a data frame 
  Rcpp::DataFrame get_all_TFs() {
    
    Rcpp::DataFrame out = Rcpp::DataFrame::create(
      Rcpp::Named("k") = Rcpp::wrap(K),
      Rcpp::Named("TF") = Rcpp::wrap(TF),
      Rcpp::Named("DTF") = Rcpp::wrap(DTF),
      Rcpp::Named("WDTF") = Rcpp::wrap(WDTF)); 
    
    return out;
  }
};

inline Rcpp::DataFrame Topographic_Functions(const arma::umat& inputADJ, const arma::umat& outputADJ, bool parallel) {
  topo_function_worker wkr(inputADJ, outputADJ); 
  if(parallel) wkr.calc_parallel(); else wkr.calc_serial(); 
  return wkr.get_all_TFs();  
}

} // close namespace 

#endif


