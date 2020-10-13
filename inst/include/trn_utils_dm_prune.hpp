#ifndef TRN_UTILS_DMPRUNE_HPP
#define TRN_UTILS_DMPRUNE_HPP

#ifndef RcppArmadillo_H
#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
#endif

#ifndef RcppParallel_H
#include <RcppParallel.h>
// [[Rcpp::depends(RcppParallel)]]
#endif

// [[Rcpp::plugins(cpp11)]]

//#ifndef SOM_UTILS_DIST_HPP
//#include "som_utils_dist.hpp"
//#endif



namespace TRN_UTILS_DMPRUNE {


// Density of a Dirichlet Multinomial distribution
inline double dDirMult(const arma::vec& x, const arma::vec& prior, bool log_form = true, bool normalize = true) {
  // x = vector of observed bin counts
  // prior = vector of prior alpha values
  
  double n = arma::accu(x); 
  if(!(n > 0)) {
    return 0.0; 
  }
  double a0 = arma::accu(prior);
  
  double out = 0.0; 
  
  out += std::lgamma(a0); 
  out += std::lgamma(n + 1);
  out -= std::lgamma(n + a0);
  
  for(int i=0; i<x.size(); ++i) {
    out += std::lgamma(x[i] + prior[i]);
    out -= std::lgamma(prior[i]);
    out -= std::lgamma(x[i] + 1);
  }
  
  // out += arma::accu(arma::lgamma(x + prior)); 
  // out -= arma::accu(arma::lgamma(prior)); 
  // out -= arma::accu(arma::lgamma(x + 1)); 
  // 
  
  
  
  // for(int i=0; i<x.size(); ++i) {
  //   out += std::lgamma(x[i] + prior[i]);
  //   out -= std::lgamma(prior[i]);
  //   out -= std::lgamma(x[i] + 1);
  // }
  
  
  
  if(normalize) out /= n;
  
  if(!log_form) out = std::exp(out);
  
  return out;
  
  
  // double n = arma::accu(x); 
  // if(!(n > 0)) {
  //   if(log_form) return min_log_probability; else return 0.0; 
  // }
  // double a0 = arma::accu(prior);
  // 
  // double out = 0.0; 
  // 
  // out += std::lgamma(a0); 
  // out += std::lgamma(n + 1);
  // out -= std::lgamma(n + a0);
  // 
  // for(int i=0; i<x.size(); ++i) {
  //   out += std::lgamma(x[i] + prior[i]);
  //   out -= std::lgamma(prior[i]);
  //   out -= std::lgamma(x[i] + 1);
  // }
  // 
  // if(normalize) out /= n;
  // 
  // if(!log_form) out = std::exp(out);
  // 
}

 
// Internal function to compute the length portion of the score score
inline arma::vec DMP_length_score(arma::uvec lengths, unsigned int max_length, unsigned int TP_radius) {

  double y_at_tp_radius = double(max_length - TP_radius + 1) / double(max_length);
  double y_at_end = 1.0 / double(max_length);

  arma::vec expy = y_at_tp_radius*arma::exp(std::log(y_at_end / y_at_tp_radius) * (arma::conv_to<arma::vec>::from(lengths) - double(TP_radius)));
  arma::vec liny = (double(max_length) - arma::conv_to<arma::vec>::from(lengths) + 1.0) / double(max_length);

  arma::uvec replace_these = arma::find(liny < expy);
  if(replace_these.size() > 0) {
    expy.elem(replace_these) = liny.elem(replace_these);
  }

  return expy;
}


struct DMPrune_Lambda_worker : public RcppParallel::Worker {

  // Inputs
  const arma::uvec& mu_rank; 
  const arma::vec& CADJ; 
  const arma::vec& CADJ_prior; 
  
  // Internals 
  unsigned int nranks;
  
  // Outputs 
  arma::uvec pruneStep; 
  arma::vec Lambda; 

  
  // Constructor
  DMPrune_Lambda_worker(const arma::uvec& mu_rank, const arma::vec& CADJ, const arma::vec& CADJ_prior)
    : mu_rank(mu_rank), CADJ(CADJ), CADJ_prior(CADJ_prior)
  {
    
    pruneStep = arma::sort(arma::unique(mu_rank)); 
    nranks = pruneStep.max(); 
    
    pruneStep -= 1; 
    // Don't compute Lambda at the end step, there is no more data in the model 
    Lambda.set_size(nranks); 
  }


  // Find TP contribution from single prototype
  void prune_step_j(unsigned int j) {
    
    arma::vec tmpCADJ = CADJ; 
    arma::uvec prune_these = arma::find(mu_rank <= j);
    if(prune_these.n_elem > 0) {
      tmpCADJ.elem(prune_these).zeros();   
    }
    
    Lambda[j] = dDirMult(tmpCADJ, CADJ_prior, true, true); 

    return;
  }


  // Parallel operator - find BMU of each row of X in parallel
  void operator()(std::size_t begin, std::size_t end) {
    for(unsigned int j = begin; j < end; j++) {
      prune_step_j(j);
    }
  }

  // Parallel call method
  void calc_parallel() {
    RcppParallel::parallelFor(0, nranks, *this);
  }

  // Non-parallel call method
  void calc_serial() {
    for(unsigned int j=1; j<nranks; ++j) {
      prune_step_j(j);
    }
  }
};

struct DMPrune_Lambda_errorband_worker : public RcppParallel::Worker {
  
  // Inputs
  const arma::uvec& mu_rank; 
  const arma::vec& CADJ; 
  const arma::vec& CADJ_prior; 
  unsigned int nperms; 
  double qlo, qhi; 
  
  // Internals 
  unsigned int nranks;
  
  // Outputs 
  arma::uvec pruneStep; 
  arma::vec lo, hi; 
  
  
  // Constructor
  DMPrune_Lambda_errorband_worker(const arma::uvec& mu_rank, const arma::vec& CADJ, const arma::vec& CADJ_prior, 
                                  unsigned int nperms_, double qlo_, double qhi_)
    : mu_rank(mu_rank), CADJ(CADJ), CADJ_prior(CADJ_prior)
  {
    
    nperms = nperms_; 
    if(!(qlo_ >= 0.0 && qlo_ <= 1.0)) Rcpp::stop("qlo must be in [0,1]");
    if(!(qhi_ >= 0.0 && qhi_ <= 1.0)) Rcpp::stop("qhi must be in [0,1]");
    if(!(qlo_ < qhi_)) Rcpp::stop("qlo must be < qhi");
    qlo = qlo_; 
    qhi = qhi_; 
    
    pruneStep = arma::sort(arma::unique(mu_rank)); 
    nranks = pruneStep.max(); 
    pruneStep -= 1; 
    
    
    // Don't compute Lambda at the end step, there is no more data in the model 
    lo.set_size(nranks); 
    hi.set_size(nranks); 
  }
  
  
  // Find TP contribution from single prototype
  void sample_prune_step_j(unsigned int j) {
    // j assumed 0-indexed 
    
    if(j == 0) {
      lo[j] = dDirMult(CADJ, CADJ_prior, true, true);
      hi[j] = lo[j]; 
    }
    
    arma::uvec sample_idx = arma::ones<arma::uvec>(mu_rank.n_elem); 
    sample_idx.elem(arma::find(mu_rank <= j)).zeros(); 
    
    arma::vec permLambda(nperms); 
    for(unsigned int perm=0; perm < nperms; ++perm) {
      arma::vec tmpCADJ = CADJ; 
      sample_idx = sample_idx.elem(arma::randperm(mu_rank.n_elem)); 
      arma::uvec prune_these = arma::find(sample_idx == 0); 
      tmpCADJ.elem(prune_these).zeros(); 
      
      permLambda[perm] = dDirMult(tmpCADJ, CADJ_prior, true, true); 
    }
    
    arma::vec qprobs(2); qprobs[0] = qlo; qprobs[1] = qhi; 
    arma::vec quants = arma::quantile(permLambda, qprobs); 
    
    lo[j] = quants[0]; 
    hi[j] = quants[1]; 

    return;
  }
  
  
  // Parallel operator - find BMU of each row of X in parallel
  void operator()(std::size_t begin, std::size_t end) {
    for(unsigned int j = begin; j < end; j++) {
      sample_prune_step_j(j);
    }
  }
  
  // Parallel call method
  void calc_parallel() {
    RcppParallel::parallelFor(0, nranks, *this);
  }
  
  // Non-parallel call method
  void calc_serial() {
    for(unsigned int j=1; j<nranks; ++j) {
      sample_prune_step_j(j);
    }
  }
};


struct DMPrune_Lambda_loo_worker : public RcppParallel::Worker {
  
  // Inputs
  const arma::vec& CADJ; 
  const arma::vec& CADJ_prior; 
  
  // Internals 
  arma::vec sparseCADJ; 
  arma::uvec loo_order; 
  arma::vec loo_llk; 
  arma::vec loo_select; 
  unsigned int n;
  
  // Constructor
  DMPrune_Lambda_loo_worker(const arma::vec& CADJ, const arma::vec& CADJ_prior)
    : CADJ(CADJ), CADJ_prior(CADJ_prior)
  {
    
    n = CADJ.n_elem;
    loo_order.set_size(n); loo_order.zeros(); // 0 indicates it has not yet been removed 
    loo_llk.set_size(n); loo_llk.zeros(); 
    loo_select.set_size(n); // will be initialized before every run 
    sparseCADJ = CADJ; 
    
  }
  
  
  // Find TP contribution from single prototype
  void test_loo_j(unsigned int j) {
    
    if(loo_order[j] > 0) return; 
    
    arma::vec tmpCADJ = sparseCADJ; 
    tmpCADJ[j] = 0; 
  
    loo_select[j] = dDirMult(tmpCADJ, CADJ_prior, true, true); 
    
    return;
  }
  
  
  // Parallel operator - find BMU of each row of X in parallel
  void operator()(std::size_t begin, std::size_t end) {
    for(unsigned int j = begin; j < end; j++) {
      test_loo_j(j);
    }
  }
  
  // Parallel call method
  void calc_parallel() {
    
    unsigned int print_counter = 0; 
    double llk0 = dDirMult(CADJ, CADJ_prior, true, true); 
    
    for(unsigned int j=1; j<n; ++j) {
      // Reset select 
      loo_select.fill(-1e10); 
      
      // Test all in parallel 
      RcppParallel::parallelFor(0, n, *this);
      
      // Find winner 
      //unsigned int remove_this = loo_select.index_max();
      unsigned int remove_this = arma::abs(loo_select - llk0).index_min();
      
      // Mark it as removed 
      loo_order[remove_this] = j; 
      loo_llk[remove_this] = loo_select[remove_this];
      sparseCADJ[remove_this] = 0; 
      //cur_llk = loo_select[remove_this];
      
      // Display 
      if(j % int(0.05*n) == 0) {
        Rcpp::Rcout << int( std::round(double(j)/double(n)*100) ) << "%\t";
        //Rcpp::Rcout << "\t" << j; 
        print_counter++; 
        if(print_counter % 5 == 0 || j == (n-1)) Rcpp::Rcout << std::endl; 
      }
    }
    
  }
  
  // // Non-parallel call method
  // void calc_serial() {
  //   for(unsigned int j=1; j<nranks; ++j) {
  //     prune_step_j(j);
  //   }
  // }
};



// // *** Topographic Function *** 
// struct topo_function_worker : public RcppParallel::Worker {
//   
//   // Inputs
//   const arma::umat& in_topoADJ; 
//   const arma::umat& out_topoADJ; 
//   
//   // Internal variables 
//   unsigned int nW; 
//   unsigned int maxK;
//   arma::umat in_geodist; 
//   arma::umat out_geodist; 
//   
//   // Output 
//   arma::ivec K; 
//   arma::vec TF; 
//   arma::vec DTF; 
//   arma::vec WDTF;
//   
//   //arma::uvec posK; 
//   //arma::ivec negK; 
//   
//   //arma::vec posTF; 
//   //arma::vec negTF;
//   //arma::vec posDTF; 
//   //arma::vec negDTF; 
//   //arma::vec posWDTF; 
//   //arma::vec negWDTF; 
// 
//   
//   // Constructor
//   topo_function_worker(const arma::umat& in_topoADJ, const arma::umat& out_topoADJ)
//     : in_topoADJ(in_topoADJ), out_topoADJ(out_topoADJ)
//   {
//     // Check that in & out topoADJ have same size 
//     if(in_topoADJ.n_rows != out_topoADJ.n_rows) Rcpp::stop("nrow(input_topoADJ) != nrow(output_topoADJ)");
//     if(in_topoADJ.n_cols != out_topoADJ.n_cols) Rcpp::stop("ncol(input_topoADJ) != ncol(output_topoADJ)");
//     
//     nW = in_topoADJ.n_rows; 
//     
//     // Compute the geodesic distances of the graphs
//     arma::umat in_binary_topoADJ = in_topoADJ + in_topoADJ.t();
//     in_geodist = SOM_UTILS_DIST::geodesicdist(in_binary_topoADJ, false, false); // (adjacency, is_weighted, is_directed)
// 
//     arma::umat out_binary_topoADJ = out_topoADJ + out_topoADJ.t();
//     out_geodist = SOM_UTILS_DIST::geodesicdist(out_binary_topoADJ, false, false);
// 
//     maxK = std::max(in_geodist.max(), out_geodist.max());
// 
//     // Setup output containers
//     K = arma::regspace<arma::ivec>(-int(maxK), int(maxK)); 
//     //posK = arma::regspace<arma::uvec>(1,  int(maxK));
//     //negK = arma::regspace<arma::ivec>(-int(maxK), -1); 
//     
//     //posTF.set_size(maxK); posTF.zeros(); 
//     //negTF.set_size(maxK); negTF.zeros(); 
//     //posDTF.set_size(maxK); posDTF.zeros(); 
//     //negDTF.set_size(maxK); negDTF.zeros(); 
//     //posWDTF.set_size(maxK); posWDTF.zeros(); 
//     //negWDTF.set_size(maxK); negWDTF.zeros(); 
//     
//     TF.set_size(K.n_elem); TF.zeros(); 
//     DTF.set_size(K.n_elem); DTF.zeros(); 
//     WDTF.set_size(K.n_elem); WDTF.zeros(); 
//   }
//   
//   
//   // Find TP contribution from single prototype 
//   // Only consider the positive folding range here 
//   void topo_function_length_k(unsigned int k) {
//     // Subscripts "M" denote manifold (Rd), "L" denote lattice 
//     
//     if(k==1) return; 
//     
//     arma::uvec find_posK = arma::find(K == int(k));
//     unsigned int idx_posK = find_posK[0];
//     arma::uvec find_negK = arma::find(K == -1*int(k));
//     unsigned int idx_negK = find_negK[0];
// 
//     // Forward DTF
//     arma::uvec fwd_fold = arma::find(out_geodist==1 && in_geodist==k);
//     DTF[idx_posK] = double(fwd_fold.n_elem) / double(nW);
//       
//     // Forward WDTF
//     WDTF[idx_posK] = arma::accu(in_topoADJ.elem(fwd_fold)) / double(arma::accu(in_topoADJ)); 
// 
//     // Backward DTF
//     arma::uvec bwd_fold = arma::find(in_geodist==1 && out_geodist==k);
//     DTF[idx_negK] = double(bwd_fold.n_elem) / double(nW);
//     
//     return;
//   }
//   
//   void topo_function_length_0() {
//     // Subscripts "M" denote manifold (Rd), "L" denote lattice 
//     
//     arma::uvec find_zeroK = arma::find(K == 0);
//     unsigned int idx_zeroK = find_zeroK[0];
// 
//     // Forward
//     arma::uvec find_Kpos1 = arma::find(K == 1);
//     arma::uvec find_Kneg1 = arma::find(K == -1);
//     TF[idx_zeroK] = TF[find_Kpos1[0]] + TF[find_Kneg1[0]];
//     
//     return;
//   }
//   
//   
//   // Parallel operator - find BMU of each row of X in parallel
//   void operator()(std::size_t begin, std::size_t end) {
//     for(unsigned int j = begin; j < end; j++) {
//       topo_function_length_k(j);
//     }
//   }
//   
//   // Parallel call method
//   void calc_parallel() {
//     RcppParallel::parallelFor(1, maxK+1, *this);
//     this->Build_TF(); 
//     this->topo_function_length_0(); 
//   }
//   
//   // Non-parallel call method
//   void calc_serial() {
//     for(unsigned int j=1; j<=maxK; ++j) {
//       topo_function_length_k(j);
//     }
//     this->Build_TF(); 
//     this->topo_function_length_0();
//   }
//   
//   void Build_TF() {
//     
//     arma::uvec find_zeroK = arma::find(K == 0);
//     unsigned int idx_zeroK = find_zeroK[0];
//     
//     for(unsigned int j=idx_zeroK+1; j < K.n_elem; ++j) {
//       for(unsigned int jj=(j+1); jj < K.n_elem; ++jj) {
//          TF[j] += DTF[jj];
//       }
//     }
//     
//     for(unsigned int j=idx_zeroK-1; j>=0; --j) {
//       for(unsigned int jj=(j-1); jj>=0; --jj) {
//         TF[j] += DTF[jj]; 
//       }
//     }
//     
//     
//   }
//   
//   Rcpp::DataFrame AllTopoFunctions() {
//     Rcpp::DataFrame out = Rcpp::DataFrame::create(
//       Rcpp::Named("k") = Rcpp::wrap(K), 
//       Rcpp::Named("TF") = Rcpp::wrap(TF), 
//       Rcpp::Named("DTF") = Rcpp::wrap(DTF), 
//       Rcpp::Named("WDTF") = Rcpp::wrap(WDTF)
//     ); 
//     return out; 
//   }
// };
// 
// inline Rcpp::DataFrame Discrete_Topographic_Function(const arma::umat& inputADJ, const arma::umat& outputADJ, bool parallel) {
//   topo_function_worker wkr(inputADJ, outputADJ); 
//   if(parallel) wkr.calc_parallel(); else wkr.calc_serial(); 
//   return wkr.AllTopoFunctions(); 
// }

} // close namespace 

#endif


