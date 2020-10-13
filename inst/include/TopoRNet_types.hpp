#ifndef TOPORNET_TRNOBJ_HPP
#define TOPORNET_TRNOBJ_HPP

#ifndef RcppArmadillo_H
#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
#endif

#ifndef RcppParallel_H
#include <RcppParallel.h>
// [[Rcpp::depends(RcppParallel)]]
#endif

// [[Rcpp::plugins(cpp11)]]

#include <string>
#include <random>
#include <algorithm>


#include "som_utils_dist.hpp"
#include "trn_utils_groupstats.hpp"
#include "trn_utils_topo_pres.hpp"
#include "trn_utils_connvis.hpp"
#include "trn_utils_dm_prune.hpp"

// RcppParallel is defining some variable named FALSE, 
// which causes errors when TRN_MODULE is exported 
#ifdef FALSE
#undef FALSE
#endif


// Forward declaration for RCPP_MODULES 
class TRN; 
RCPP_EXPOSED_CLASS(TRN);



class TRN{
  
public:
  bool parallel; 
  void set_parallel(bool parallel_); 
  
  // Graphs stored internally
  unsigned int nV;      // # of vertices in the graph 
  unsigned int CADJ_nE; // # of CADJ edges 
  unsigned int CONN_nE; // # of CONN edges
  unsigned int OT_nE;   // # of lattice edges 
  
  arma::umat CADJ_EL;   // 1-indexed CADJ edge list
  arma::umat CONN_EL;   // 1-indexed CONN edge list 
  arma::umat OT_EL;     // 1-indexed lattice edge list 
  
  arma::vec CADJ;       // CADJ edge weights for each edge in CADJ_EL
  arma::vec CONN;       // CONN edge weights for each edge in CONN_EL 
  
  arma::uvec CADJ_deg;  // Unweighted vertex degrees on CADJ graph 
  arma::uvec CONN_deg;  // Unweighted vertex degrees on CONN graph 
  arma::uvec OT_deg;    // Unweighted vertex degrees on output topology graph 
  arma::vec CADJ_wdeg;  // Weighted vertex degrees on CADJ graph 
  arma::vec CONN_wdeg;  // Weighted vertex degrees on CONN graph 
  
  arma::uvec CADJ_fwdflen; // Forward folding length of each CADJ edge, length = nrow(CADJ_EL)
  arma::uvec CONN_fwdflen; // Forward folding length of each CONN edge, length = nrow(CONN_EL)
  arma::uvec CADJ_bwdflen; // Backward folding length of each lattice edge, length = nrow(OT_EL)
  arma::uvec CONN_bwdflen; // Backward folding length of each lattice edge, length = nrow(OT_EL)
  
  unsigned int CADJ_TP_radius; // Topology Preserving Radius according to CADJ 
  unsigned int CONN_TP_radius; // Topology Preserving Radius according to CONN 
  unsigned int calc_TP_radius(arma::uvec vertex_degrees); 
  
  arma::umat OT_geodesic_dist; // The geodesic distance matrix of the output topology 
  unsigned int OT_max_dist;    // The maximum geodesic distance between any two TRN vertices, measured according to its output topology
  
  arma::umat OT_nhb_sizes;     // an (nV x OT_max_dist) matrix storing the # of lattice neighbors of each neuron at each distance in [1,OT_max_dist]
  void calc_OT_nhb_sizes(); 
  
  void set_CADJ(const arma::umat& CADJ_); 
  void set_CONN(); 
  
  // Constructor
  TRN(); 
  
  // Initialize with input & output topologies, defined by their adjacency matrices 
  void initialize_TRN(arma::umat CADJ_, arma::umat OTADJ_);
  
  // CONNvis fields and methods
  arma::uvec CADJ_lrank;
  arma::uvec CADJ_grank;
  TRN_UTILS_GROUPSTATS::groupstat_container CADJ_CVstats_lrank;
  TRN_UTILS_GROUPSTATS::groupstat_container CADJ_CVstats_grank;
  TRN_UTILS_GROUPSTATS::groupstat_container CADJ_CVstats_length;
  void calc_CADJvis_stats();
  Rcpp::DataFrame get_CADJvis_stats(std::string which_stat); 
  std::vector<std::string> get_CADJvis_colors(); 
  arma::vec get_CADJvis_widths(); 
  
  arma::uvec CONN_lrank;
  arma::uvec CONN_grank;
  TRN_UTILS_GROUPSTATS::groupstat_container CONN_CVstats_lrank;
  TRN_UTILS_GROUPSTATS::groupstat_container CONN_CVstats_grank;
  TRN_UTILS_GROUPSTATS::groupstat_container CONN_CVstats_length;
  void calc_CONNvis_stats();
  Rcpp::DataFrame get_CONNvis_stats(std::string which_stat); 
  std::vector<std::string> get_CONNvis_colors(); 
  arma::vec get_CONNvis_widths(); 

  // Pruning
  arma::uvec CADJ_active; 
  arma::uvec CONN_active; 
  void restore_CADJ_edges(); 
  void restore_CONN_edges(); 
  
  void prune_CADJ_edge_list(const arma::umat& edge_list);
  void prune_CADJ_edge_id(const arma::uvec& edge_ids);
  void prune_CADJ_vertex_id(const arma::uvec& vertex_ids);
  void prune_CADJ_edge_weight(double min_weight); // tv equivalent, by value
  void prune_CADJ_lrank(unsigned int max_lrank); // tn equivalent
  void prune_CADJ_grank(unsigned int max_grank); // tv equivalent, but by rank instead of value
  void prune_CADJ_length(unsigned int max_len);  // tl equivalent
  void prune_CADJ_wdeg(double min_weight); // tr equivalent
  arma::mat get_CADJ(); 

  void prune_CONN_edge_list(const arma::umat& edge_list);
  void prune_CONN_edge_id(const arma::uvec& edge_ids);
  void prune_CONN_vertex_id(const arma::uvec& vertex_ids);
  void prune_CONN_edge_weight(double min_weight); // tv equivalent, by value
  void prune_CONN_lrank(unsigned int max_lrank); // tn equivalent
  void prune_CONN_grank(unsigned int max_grank); // tv equivalent, but by rank instead of value
  void prune_CONN_length(unsigned int max_len);  // tl equivalent
  void prune_CONN_wdeg(double min_weight); // tr equivalent
  void prune_CONN_CADJ(); 
  arma::mat get_CONN(); 
  
  arma::mat get_OTADJ(); 
  
  arma::uvec CADJ_active_verts();
  arma::uvec CONN_active_verts(); 
  
  // Topology Preservation Measures 
  double TopoProd; 
  Rcpp::DataFrame CADJ_TopoFxns; 
  Rcpp::DataFrame CONN_TopoFxns; 
  void calc_TopoMeasures(const arma::mat& inV, const arma::mat& outV); 
  
  // DM-Prune 
  arma::vec CADJ_BootSig; 
  arma::vec CONN_BootSig; 
  void set_BootSig(const arma::mat& BSADJ); 
  
  arma::vec DMP_prior; 
  double DMP_gamma_G, DMP_gamma_N, DMP_gamma_L, DMP_gamma_S; 
  arma::vec DMP_mu_G, DMP_mu_N, DMP_mu_L, DMP_mu_S, DMP_mu;
  arma::uvec DMP_mu_rank; 
  arma::uvec DMP_pruneStep; 
  arma::vec DMP_Lambda; 
  void calc_DMPrune_LambdaPath(const arma::mat& priorADJ); 
  
  // These are exposed but currently undocumented!! 
  unsigned int DMP_Lambda_error_nperms; 
  double DMP_Lambda_error_qlo, DMP_Lambda_error_qhi; 
  arma::mat DMP_Lambda_error; 
  void calc_DMPrune_LambdaPath_error(); 
  
  arma::uvec DMP_loo_order; 
  arma::vec DMP_loo_llk; 
  void calc_DMPrune_LambdaPath_loo(const arma::mat& UPADJ); 
  
  void DMPrune_CADJ_step(unsigned int min_step); 
  
  // Flags 
  bool isinit; 
  bool isset_TPM; 
  bool isset_DMPrune; 
    
    
  // IO 
  void save(std::string rdsfile); 
  void load(std::string rdsfile); 
  Rcpp::List as_list(); 
  void load_list(Rcpp::List TRNList); 
  
};


// ***** Constructor 
inline TRN::TRN() {
  
  // Defaults 
  this->parallel = true; 
  
  this->DMP_gamma_G = 1.0; 
  this->DMP_gamma_N = 1.0; 
  this->DMP_gamma_L = 1.0; 
  this->DMP_gamma_S = 1.0; 
  this->DMP_Lambda_error_nperms = 1000; 
  this->DMP_Lambda_error_qlo = 0.025; 
  this->DMP_Lambda_error_qhi = 0.975; 
  
  //this->isset_in_topo = false; 
  //this->isset_out_topo = false; 
  
  //this->isset_CV = false; 
  this->isinit = false; 
  this->isset_TPM = false;
  
  this->isset_DMPrune = false; 
  
}

inline void TRN::set_parallel(bool parallel_) {
  this->parallel = parallel_; 
}

inline void TRN::calc_OT_nhb_sizes() {
  
  this->OT_nhb_sizes.set_size(this->nV, this->OT_max_dist); 
  
  // Loop over each vertex 
  for(unsigned int neuroni = 0; neuroni < this->nV; ++neuroni) {
    
    // Loop over each distance 
    for(unsigned int distj = 1; distj <= this->OT_max_dist; ++distj) {
      arma::uvec num_neurons_at_distj = arma::find(this->OT_geodesic_dist.row(neuroni) == distj); 
      this->OT_nhb_sizes(neuroni, distj-1) = num_neurons_at_distj.n_elem; 
    }
    
    // Compute cumsum of neighbors at each distance 
    this->OT_nhb_sizes.row(neuroni) = arma::cumsum(this->OT_nhb_sizes.row(neuroni));
  }
  
}

inline unsigned int TRN::calc_TP_radius(arma::uvec vertex_degrees) {
  
  unsigned int TP_radius = 0; 
  
  for(unsigned int neuroni=0; neuroni<this->nV; ++neuroni) {
    for(unsigned int distj = 1; distj <= this->OT_max_dist; ++distj) {
      
      if(vertex_degrees[neuroni] <= this->OT_nhb_sizes(neuroni, distj-1)) {
        TP_radius = std::max(TP_radius, distj); 
        break; 
      }
      
    }
  }
  
  return TP_radius; 
}

inline void TRN::set_CADJ(const arma::umat& CADJ_) {
  
  // *** Storing CADJ Topology 
  Rcpp::Rcout << "Storing CADJ Topology:" << std::endl; 
  
  Rcpp::Rcout << "++ edge list ... "; 
  arma::uvec CADJ_nonzero = arma::find(CADJ_ > 0);
  this->CADJ_EL = arma::ind2sub(arma::size(CADJ_), CADJ_nonzero).t() + 1; // store edgelist as 1-based index
  Rcpp::Rcout << "done" << std::endl; 
  
  this->CADJ_nE = CADJ_nonzero.n_elem;
  Rcpp::Rcout << "   " << this->CADJ_nE << " edges stored" << std::endl;
  
  Rcpp::Rcout << "++ weights ... "; 
  this->CADJ = arma::conv_to<arma::vec>::from(CADJ_(CADJ_nonzero));
  Rcpp::Rcout << "done" << std::endl; 
  
  Rcpp::Rcout << "++ vertex degrees ... ";
  this->CADJ_wdeg = arma::conv_to<arma::vec>::from(arma::sum(CADJ_, 1));
  
  arma::umat binaryCADJ = 0*CADJ_;
  binaryCADJ.elem(CADJ_nonzero).ones();
  this->CADJ_deg = arma::sum(binaryCADJ, 1);
  Rcpp::Rcout << " done" << std::endl;
  
  Rcpp::Rcout << "++ forward folding lengths ... ";
  this->CADJ_fwdflen = this->OT_geodesic_dist.elem(arma::sub2ind(arma::size(CADJ_), this->CADJ_EL.t()-1));
  Rcpp::Rcout << "done" << std::endl; 
  
  Rcpp::Rcout << "++ backward folding lengths ... ";
  arma::umat CADJ_geo_dist = SOM_UTILS_DIST::geodesicdist(binaryCADJ, false, false);
  this->CADJ_bwdflen = CADJ_geo_dist.elem(arma::sub2ind(arma::size(CADJ_), this->OT_EL.t()-1));
  Rcpp::Rcout << "done" << std::endl;
  
  Rcpp::Rcout << "++ topology preserving radius ... "; 
  this->CADJ_TP_radius = this->calc_TP_radius(this->CADJ_deg); 
  Rcpp::Rcout << "done" << std::endl;
  Rcpp::Rcout << "   CADJ_TP_radius = " << this->CADJ_TP_radius << std::endl; 
  
  Rcpp::Rcout << "++ activating all edges ... "; 
  this->restore_CADJ_edges(); 
  Rcpp::Rcout << "done" << std::endl; 
  
  Rcpp::Rcout << "----------------------------------------------------------------" << std::endl;
  
  // *** CADJVis Stats
  this->calc_CADJvis_stats();
  Rcpp::Rcout << "Call $get_CADJvis_stats to view" << std::endl;
}

inline void TRN::set_CONN() {
  
  // *** Storing CONN Topology 
  Rcpp::Rcout << "Storing CONN Topology:" << std::endl; 
  
  arma::umat CONN_ = arma::conv_to<arma::umat>::from(this->get_CADJ()); 
  CONN_+= CONN_.t(); 
  
  Rcpp::Rcout << "++ edge list ... "; 
  arma::uvec CONN_nonzero = arma::find(CONN_ > 0);
  this->CONN_EL = arma::ind2sub(arma::size(CONN_), CONN_nonzero).t() + 1; // store edgelist as 1-based index
  Rcpp::Rcout << "done" << std::endl; 
  
  this->CONN_nE = CONN_nonzero.n_elem;
  Rcpp::Rcout << "   " << this->CONN_nE << " edges stored" << std::endl;
  
  Rcpp::Rcout << "++ weights ... "; 
  this->CONN = arma::conv_to<arma::vec>::from(CONN_(CONN_nonzero));
  Rcpp::Rcout << "done" << std::endl; 
  
  Rcpp::Rcout << "++ vertex degrees ... ";
  this->CONN_wdeg = arma::conv_to<arma::vec>::from(arma::sum(CONN_, 1));
  
  arma::umat binaryCONN = 0*CONN_;
  binaryCONN.elem(CONN_nonzero).ones();
  this->CONN_deg = arma::sum(binaryCONN, 1);
  Rcpp::Rcout << " done" << std::endl;
  
  Rcpp::Rcout << "++ forward folding lengths ... ";
  this->CONN_fwdflen = this->OT_geodesic_dist.elem(arma::sub2ind(arma::size(CONN_), this->CONN_EL.t()-1));
  Rcpp::Rcout << "done" << std::endl; 
  
  Rcpp::Rcout << "++ backward folding lengths ... ";
  arma::umat CONN_geo_dist = SOM_UTILS_DIST::geodesicdist(binaryCONN, false, false);
  this->CONN_bwdflen = CONN_geo_dist.elem(arma::sub2ind(arma::size(CONN_), this->OT_EL.t()-1));
  Rcpp::Rcout << "done" << std::endl;
  
  Rcpp::Rcout << "++ topology preserving radius ... "; 
  this->CONN_TP_radius = this->calc_TP_radius(this->CONN_deg); 
  Rcpp::Rcout << "done" << std::endl;
  Rcpp::Rcout << "   CONN_TP_radius = " << this->CONN_TP_radius << std::endl; 
  
  Rcpp::Rcout << "++ activating all edges ... "; 
  this->restore_CONN_edges(); 
  Rcpp::Rcout << "done" << std::endl; 
  
  Rcpp::Rcout << "----------------------------------------------------------------" << std::endl;
  
  
  // *** CONNVis Stats
  this->calc_CONNvis_stats();
  Rcpp::Rcout << "Call $get_CONNvis_stats to view" << std::endl;
}

inline void TRN::initialize_TRN(arma::umat CADJ_, arma::umat OTADJ_) {
  
  // Set number of vertices in the graph, make sure both adjacencies are square and the same dimension 
  this->nV = CADJ_.n_rows; 
  if(OTADJ_.n_rows != this->nV) Rcpp::stop("nrow(CADJ) != nrow(OTADJ)");
  if(CADJ_.n_cols != this->nV) Rcpp::stop("CADJ must be square");
  if(OTADJ_.n_cols != this->nV) Rcpp::stop("OTADJ must be square");
  
  Rcpp::Rcout << "Initializing TRN object:" << std::endl; 
  Rcpp::Rcout << "----------------------------------------------------------------" << std::endl;
  
  
  
  // *** Store Output Topology 
  Rcpp::Rcout << "Storing Output Topology:" << std::endl; 
  
  Rcpp::Rcout << "++ edge list ... "; 
  arma::uvec OT_nonzero = arma::find(OTADJ_ > 0);
  this->OT_EL = arma::ind2sub(arma::size(OTADJ_), OT_nonzero).t() + 1; // store edgelist as 1-based index 
  Rcpp::Rcout << "done" << std::endl; 
  
  this->OT_nE = OT_nonzero.n_elem;
  Rcpp::Rcout << "   " << this->OT_nE << " edges stored" << std::endl;
  
  Rcpp::Rcout << "++ vertex degrees ... "; 
  this->OT_deg = arma::sum(OTADJ_, 1); 
  Rcpp::Rcout << " done" << std::endl; 
  
  Rcpp::Rcout << "++ geodesic distances ... "; 
  this->OT_geodesic_dist = SOM_UTILS_DIST::geodesicdist(OTADJ_, false, false);
  this->OT_max_dist = this->OT_geodesic_dist.max(); 
  Rcpp::Rcout << " done" << std::endl; 
  Rcpp::Rcout << "   max dist = " << this->OT_max_dist << std::endl; 
  
  Rcpp::Rcout << "++ neighborhood sizes ... "; 
  this->calc_OT_nhb_sizes(); 
  Rcpp::Rcout << "done" << std::endl; 
  if(arma::min(arma::max(this->OT_nhb_sizes, 1)) < (this->nV - 1)) {
    Rcpp::Rcout << "   warning: Output Topology is not fully connected" << std::endl; 
  }
  
  Rcpp::Rcout << "----------------------------------------------------------------" << std::endl;
  
  
  
  // *** Storing CADJ Topology 
  this->set_CADJ(CADJ_); 
  
  //Rcpp::Rcout << "----------------------------------------------------------------" << std::endl;
  
  // *** CADJVis Stats
  //this->calc_CADJvis_stats();
  //Rcpp::Rcout << "Call $get_CADJvis_stats to view" << std::endl;

  Rcpp::Rcout << "----------------------------------------------------------------" << std::endl;
  
  // *** Storing CONN Topology 
  this->set_CONN(); 
  
  //Rcpp::Rcout << "----------------------------------------------------------------" << std::endl;
  

  // *** CONNVis Stats
  //this->calc_CONNvis_stats();
  //Rcpp::Rcout << "Call $get_CONNvis_stats to view" << std::endl;
  
  Rcpp::Rcout << "----------------------------------------------------------------" << std::endl;

  this->isinit = true; 
}


// ***** CONNvis Local Rank
inline void TRN::calc_CONNvis_stats() {

  Rcpp::Rcout << "Calculating CONNvis statistics" << std::endl;

  
  // *** CONN local ranks & their stats 
  Rcpp::Rcout << "++ CONN local rank stats ... ";
  TRN_UTILS_CONNVIS::CONNvis_lrank_worker cv_lrank_wkr(this->CONN_EL, this->CONN);
  if(parallel) cv_lrank_wkr.calc_parallel(); else cv_lrank_wkr.calc_serial();
  this->CONN_lrank = cv_lrank_wkr.lrank;
  
  TRN_UTILS_GROUPSTATS::groupstat_worker connstat_local_wkr(this->CONN_lrank, this->CONN);
  connstat_local_wkr.set_stat("all");
  if(parallel) connstat_local_wkr.calc_parallel(); else connstat_local_wkr.calc_serial();
  this->CONN_CVstats_lrank = connstat_local_wkr.groupstats; 
  
  Rcpp::Rcout << "done" << std::endl;

  
  // *** CONN global ranks & their stats 
  Rcpp::Rcout << "++ CONN global rank stats ... ";
  this->CONN_grank = TRN_UTILS_CONNVIS::CONNvis_grank(this->CONN, this->CONN_CVstats_lrank.mean); 
  
  TRN_UTILS_GROUPSTATS::groupstat_worker connstat_global_wkr(this->CONN_grank, this->CONN);
  connstat_global_wkr.set_stat("all");
  if(parallel) connstat_global_wkr.calc_parallel(); else connstat_global_wkr.calc_serial();
  this->CONN_CVstats_grank = connstat_global_wkr.groupstats; 
  
  Rcpp::Rcout << "done" << std::endl;

  
  // *** CADJ length stats 
  Rcpp::Rcout << "++ CONN length stats ... ";

  TRN_UTILS_GROUPSTATS::groupstat_worker connstat_length_wkr(this->CONN_fwdflen, this->CONN);
  connstat_length_wkr.set_stat("all");
  if(parallel) connstat_length_wkr.calc_parallel(); else connstat_length_wkr.calc_serial();
  this->CONN_CVstats_length = connstat_length_wkr.groupstats; 
  
  Rcpp::Rcout << "done" << std::endl;
}

inline void TRN::calc_CADJvis_stats() {
  
  Rcpp::Rcout << "Calculating CADJvis statistics" << std::endl;
  
  
  // *** CADJ local ranks & their stats 
  Rcpp::Rcout << "++ CADJ local rank stats ... ";
  TRN_UTILS_CONNVIS::CONNvis_lrank_worker cv_lrank_wkr(this->CADJ_EL, this->CADJ);
  if(parallel) cv_lrank_wkr.calc_parallel(); else cv_lrank_wkr.calc_serial();
  this->CADJ_lrank = cv_lrank_wkr.lrank;
  
  TRN_UTILS_GROUPSTATS::groupstat_worker cadjstat_local_wkr(this->CADJ_lrank, this->CADJ);
  cadjstat_local_wkr.set_stat("all");
  if(parallel) cadjstat_local_wkr.calc_parallel(); else cadjstat_local_wkr.calc_serial();
  this->CADJ_CVstats_lrank = cadjstat_local_wkr.groupstats; 
  
  Rcpp::Rcout << "done" << std::endl;
  
  
  // *** CADJ global ranks & their stats 
  Rcpp::Rcout << "++ CADJ global rank stats ... ";
  this->CADJ_grank = TRN_UTILS_CONNVIS::CONNvis_grank(this->CADJ, this->CADJ_CVstats_lrank.mean); 
  
  TRN_UTILS_GROUPSTATS::groupstat_worker cadjstat_global_wkr(this->CADJ_grank, this->CADJ);
  cadjstat_global_wkr.set_stat("all");
  if(parallel) cadjstat_global_wkr.calc_parallel(); else cadjstat_global_wkr.calc_serial();
  this->CADJ_CVstats_grank = cadjstat_global_wkr.groupstats; 
  
  Rcpp::Rcout << "done" << std::endl;
  
  
  // *** CADJ length stats 
  Rcpp::Rcout << "++ CADJ length stats ... ";
  
  TRN_UTILS_GROUPSTATS::groupstat_worker cadjstat_length_wkr(this->CADJ_fwdflen, this->CADJ);
  cadjstat_length_wkr.set_stat("all");
  if(parallel) cadjstat_length_wkr.calc_parallel(); else cadjstat_length_wkr.calc_serial();
  this->CADJ_CVstats_length = cadjstat_length_wkr.groupstats; 
  
  Rcpp::Rcout << "done" << std::endl;
}

inline Rcpp::DataFrame TRN::get_CADJvis_stats(std::string which_stat) {
  // ** Checks
  if(!this->isinit) Rcpp::stop("Must call $initialize_TRN before calling $get_CADJvis_stats");

  // Initialize empty data frame
  Rcpp::DataFrame out;

  // Build & Return the requested stats
  if(which_stat == "lrank") {

    out = Rcpp::DataFrame::create(
      Rcpp::Named("lrank") = Rcpp::wrap(this->CADJ_CVstats_lrank.group_id),
      Rcpp::Named("count") = Rcpp::wrap(this->CADJ_CVstats_lrank.count),
      Rcpp::Named("pct") = Rcpp::wrap(this->CADJ_CVstats_lrank.pct),
      Rcpp::Named("cumpct") = Rcpp::wrap(this->CADJ_CVstats_lrank.cumpct),
      Rcpp::Named("mean") = Rcpp::wrap(this->CADJ_CVstats_lrank.mean),
      Rcpp::Named("sd") = Rcpp::wrap(this->CADJ_CVstats_lrank.sd),
      Rcpp::Named("q0") = Rcpp::wrap(this->CADJ_CVstats_lrank.q0),
      Rcpp::Named("q25") = Rcpp::wrap(this->CADJ_CVstats_lrank.q25),
      Rcpp::Named("q50") = Rcpp::wrap(this->CADJ_CVstats_lrank.q50),
      Rcpp::Named("q75") = Rcpp::wrap(this->CADJ_CVstats_lrank.q75),
      Rcpp::Named("q100") = Rcpp::wrap(this->CADJ_CVstats_lrank.q100));

  } else if(which_stat == "grank") {

    out = Rcpp::DataFrame::create(
      Rcpp::Named("grank") = Rcpp::wrap(this->CADJ_CVstats_grank.group_id),
      Rcpp::Named("count") = Rcpp::wrap(this->CADJ_CVstats_grank.count),
      Rcpp::Named("pct") = Rcpp::wrap(this->CADJ_CVstats_grank.pct),
      Rcpp::Named("cumpct") = Rcpp::wrap(this->CADJ_CVstats_grank.cumpct),
      Rcpp::Named("mean") = Rcpp::wrap(this->CADJ_CVstats_grank.mean),
      Rcpp::Named("sd") = Rcpp::wrap(this->CADJ_CVstats_grank.sd),
      Rcpp::Named("q0") = Rcpp::wrap(this->CADJ_CVstats_grank.q0),
      Rcpp::Named("q25") = Rcpp::wrap(this->CADJ_CVstats_grank.q25),
      Rcpp::Named("q50") = Rcpp::wrap(this->CADJ_CVstats_grank.q50),
      Rcpp::Named("q75") = Rcpp::wrap(this->CADJ_CVstats_grank.q75),
      Rcpp::Named("q100") = Rcpp::wrap(this->CADJ_CVstats_grank.q100));

  } else if(which_stat == "length") {

    out = Rcpp::DataFrame::create(
      Rcpp::Named("length") = Rcpp::wrap(this->CADJ_CVstats_length.group_id),
      Rcpp::Named("count") = Rcpp::wrap(this->CADJ_CVstats_length.count),
      Rcpp::Named("pct") = Rcpp::wrap(this->CADJ_CVstats_length.pct),
      Rcpp::Named("cumpct") = Rcpp::wrap(this->CADJ_CVstats_length.cumpct),
      Rcpp::Named("mean") = Rcpp::wrap(this->CADJ_CVstats_length.mean),
      Rcpp::Named("sd") = Rcpp::wrap(this->CADJ_CVstats_length.sd),
      Rcpp::Named("q0") = Rcpp::wrap(this->CADJ_CVstats_length.q0),
      Rcpp::Named("q25") = Rcpp::wrap(this->CADJ_CVstats_length.q25),
      Rcpp::Named("q50") = Rcpp::wrap(this->CADJ_CVstats_length.q50),
      Rcpp::Named("q75") = Rcpp::wrap(this->CADJ_CVstats_length.q75),
      Rcpp::Named("q100") = Rcpp::wrap(this->CADJ_CVstats_length.q100));

  } else {
    Rcpp::stop("Unknown stat: options are 'lrank','grank','length'");
  }

  return out;
}

inline Rcpp::DataFrame TRN::get_CONNvis_stats(std::string which_stat) {
  // ** Checks
  if(!this->isinit) Rcpp::stop("Must call $initialize_TRN before calling $get_CONNvis_stats");
  
  // Initialize empty data frame
  Rcpp::DataFrame out;
  
  // Build & Return the requested stats
  if(which_stat == "lrank") {
    
    out = Rcpp::DataFrame::create(
      Rcpp::Named("lrank") = Rcpp::wrap(this->CONN_CVstats_lrank.group_id),
      Rcpp::Named("count") = Rcpp::wrap(this->CONN_CVstats_lrank.count),
      Rcpp::Named("pct") = Rcpp::wrap(this->CONN_CVstats_lrank.pct),
      Rcpp::Named("cumpct") = Rcpp::wrap(this->CONN_CVstats_lrank.cumpct),
      Rcpp::Named("mean") = Rcpp::wrap(this->CONN_CVstats_lrank.mean),
      Rcpp::Named("sd") = Rcpp::wrap(this->CONN_CVstats_lrank.sd),
      Rcpp::Named("q0") = Rcpp::wrap(this->CONN_CVstats_lrank.q0),
      Rcpp::Named("q25") = Rcpp::wrap(this->CONN_CVstats_lrank.q25),
      Rcpp::Named("q50") = Rcpp::wrap(this->CONN_CVstats_lrank.q50),
      Rcpp::Named("q75") = Rcpp::wrap(this->CONN_CVstats_lrank.q75),
      Rcpp::Named("q100") = Rcpp::wrap(this->CONN_CVstats_lrank.q100));
    
  } else if(which_stat == "grank") {
    
    out = Rcpp::DataFrame::create(
      Rcpp::Named("grank") = Rcpp::wrap(this->CONN_CVstats_grank.group_id),
      Rcpp::Named("count") = Rcpp::wrap(this->CONN_CVstats_grank.count),
      Rcpp::Named("pct") = Rcpp::wrap(this->CONN_CVstats_grank.pct),
      Rcpp::Named("cumpct") = Rcpp::wrap(this->CONN_CVstats_grank.cumpct),
      Rcpp::Named("mean") = Rcpp::wrap(this->CONN_CVstats_grank.mean),
      Rcpp::Named("sd") = Rcpp::wrap(this->CONN_CVstats_grank.sd),
      Rcpp::Named("q0") = Rcpp::wrap(this->CONN_CVstats_grank.q0),
      Rcpp::Named("q25") = Rcpp::wrap(this->CONN_CVstats_grank.q25),
      Rcpp::Named("q50") = Rcpp::wrap(this->CONN_CVstats_grank.q50),
      Rcpp::Named("q75") = Rcpp::wrap(this->CONN_CVstats_grank.q75),
      Rcpp::Named("q100") = Rcpp::wrap(this->CONN_CVstats_grank.q100));
    
  } else if(which_stat == "length") {
    
    out = Rcpp::DataFrame::create(
      Rcpp::Named("length") = Rcpp::wrap(this->CONN_CVstats_length.group_id),
      Rcpp::Named("count") = Rcpp::wrap(this->CONN_CVstats_length.count),
      Rcpp::Named("pct") = Rcpp::wrap(this->CONN_CVstats_length.pct),
      Rcpp::Named("cumpct") = Rcpp::wrap(this->CONN_CVstats_length.cumpct),
      Rcpp::Named("mean") = Rcpp::wrap(this->CONN_CVstats_length.mean),
      Rcpp::Named("sd") = Rcpp::wrap(this->CONN_CVstats_length.sd),
      Rcpp::Named("q0") = Rcpp::wrap(this->CONN_CVstats_length.q0),
      Rcpp::Named("q25") = Rcpp::wrap(this->CONN_CVstats_length.q25),
      Rcpp::Named("q50") = Rcpp::wrap(this->CONN_CVstats_length.q50),
      Rcpp::Named("q75") = Rcpp::wrap(this->CONN_CVstats_length.q75),
      Rcpp::Named("q100") = Rcpp::wrap(this->CONN_CVstats_length.q100));
    
  } else {
    Rcpp::stop("Unknown stat: options are 'lrank','grank','length'");
  }
  
  return out;
}

inline std::vector<std::string> TRN::get_CADJvis_colors() {
  if(!this->isinit) Rcpp::stop("Must call $initialize_TRN before calling $get_CADJvis_colors");
  
  unsigned int max_lrank = this->CADJ_lrank.max(); 
  
  std::vector<std::string> color_map;
  
  if(max_lrank >= 1) {
    color_map.push_back("red"); 
    //if(max_lrank == 1) return color_map; 
  }
  
  if(max_lrank >= 2) {
    color_map.push_back("blue"); 
    //if(max_lrank == 2) return color_map; 
  }
  
  if(max_lrank >= 3) {
    color_map.push_back("green"); 
    //if(max_lrank == 3) return color_map; 
  }
  
  if(max_lrank >= 4) {
    color_map.push_back("yellow"); 
    //if(max_lrank == 4) return color_map; 
  }
  
  int remaining_lranks = max_lrank - 4; 
  if(remaining_lranks > 0) {
    arma::vec gray_vals = arma::linspace<arma::vec>(25.0, 75.0, remaining_lranks); 
    for(unsigned int i=0; i<(unsigned int)remaining_lranks; ++i) {
      color_map.push_back("grey" + std::to_string(int(std::round(gray_vals[i]))));
    }
  }
  
  // Fill up return colors 
  std::vector<std::string> colors(this->CADJ_nE); 
  for(unsigned int i=0; i<this->CADJ_nE; ++i) {
    colors[i] = color_map[this->CADJ_lrank[i]-1];
  }
  return colors; 
}

inline std::vector<std::string> TRN::get_CONNvis_colors() {
  if(!this->isinit) Rcpp::stop("Must call $initialize_TRN before calling $get_CONNvis_colors");
  
  unsigned int max_lrank = this->CONN_lrank.max(); 
  
  std::vector<std::string> color_map;
  
  if(max_lrank >= 1) {
    color_map.push_back("red"); 
    //if(max_lrank == 1) return color_map; 
  }
  
  if(max_lrank >= 2) {
    color_map.push_back("blue"); 
    //if(max_lrank == 2) return color_map; 
  }
  
  if(max_lrank >= 3) {
    color_map.push_back("green"); 
    //if(max_lrank == 3) return color_map; 
  }
  
  if(max_lrank >= 4) {
    color_map.push_back("yellow"); 
    //if(max_lrank == 4) return color_map; 
  }
  
  int remaining_lranks = max_lrank - 4; 
  if(remaining_lranks > 0) {
    arma::vec gray_vals = arma::linspace<arma::vec>(25.0, 75.0, remaining_lranks); 
    for(unsigned int i=0; i<(unsigned int)remaining_lranks; ++i) {
      color_map.push_back("grey" + std::to_string(int(std::round(gray_vals[i]))));
    }
  }
  
  // Fill up return colors 
  std::vector<std::string> colors(this->CONN_nE); 
  for(unsigned int i=0; i<this->CONN_nE; ++i) {
    colors[i] = color_map[this->CONN_lrank[i]-1];
  }
  return colors; 
}

inline arma::vec TRN::get_CADJvis_widths() {
  if(!this->isinit) Rcpp::stop("Must call $initialize_TRN before calling $get_CADJvis_widths");
  
  return arma::conv_to<arma::vec>::from(this->CADJ_grank.max() + this->CADJ_grank.min() - this->CADJ_grank); 
}

inline arma::vec TRN::get_CONNvis_widths() {
  if(!this->isinit) Rcpp::stop("Must call $initialize_TRN before calling $get_CONNvis_widths");
  
  return arma::conv_to<arma::vec>::from(this->CONN_grank.max() + this->CONN_grank.min() - this->CONN_grank); 
}


// // ***** Pruning
inline void TRN::restore_CADJ_edges() {
  this->CADJ_active.set_size(this->CADJ_nE);
  this->CADJ_active.fill(1);
}

inline void TRN::restore_CONN_edges() {
  this->CONN_active.set_size(this->CONN_nE);
  this->CONN_active.fill(1);
}


inline void TRN::prune_CADJ_edge_list(const arma::umat& edge_list) {
  // edge_list is assumed to contain (v1,v2) pairs in rows, which will be matched and inactivated
  if(!this->isinit) Rcpp::stop("Must call $initialize_TRN before calling $prune_CADJ_edge_list");

  if(edge_list.n_cols != 2) Rcpp::stop("Input edge_list of wrong format. Assumed to be a matrix with (v1,v2) pairs in rows");

  bool issue_missing_warning = false;
  std::string missing_list = ""; 

  // Loop over given edgelist
  unsigned int nremoved = 0; 
  for(unsigned int i=0; i<edge_list.n_rows; ++i) {

    // Find this edge in the TRG EL. If it doesn't exist, issue warning and skip
    arma::uvec matchidx = arma::find(this->CADJ_EL.col(0)==edge_list(i,0) && this->CADJ_EL.col(1)==edge_list(i,1));
    if(matchidx.n_elem == 0) {
      issue_missing_warning = true;
      missing_list += "(" + std::to_string(edge_list(i,0)) + "," + std::to_string(edge_list(i,1)) + ") "; 
      continue;
    } else {
      // Otherwise inactive the edge
      this->CADJ_active[matchidx[0]] = 0;
      nremoved++; 
    }

  }
  
  Rcpp::Rcout << "Set " << nremoved << " CADJ edges inactive." << std::endl;

  if(issue_missing_warning) {
    std::string warning_msg = "The following edges were not found in TRN$CADJ_EL (ignored):\n" + missing_list; 
    Rcpp::warning(warning_msg);
  } 
  return;
}

inline void TRN::prune_CADJ_edge_id(const arma::uvec& edge_ids) {
  // edge_ids is assumed to contain (1-based) indices to the rows of EL to inactive
  if(!this->isinit) Rcpp::stop("Must call $initialize_TRN before calling $prune_CADJ_edge_id");

  bool issue_missing_warning = false;
  std::string missing_list = ""; 
  
  // Loop over given edge ids
  unsigned int nremoved = 0; 
  for(unsigned int i=0; i<edge_ids.n_elem; ++i) {
    
    // If edge_id in range, set it to 0. 
    // Otherwise, log it as missing 
    if(edge_ids[i] > this->CADJ_nE || edge_ids[i] == 0) {
      issue_missing_warning = true;
      missing_list += std::to_string(edge_ids[i]) + ", ";
    } else {
      this->CADJ_active[edge_ids[i]-1] = 0; 
      nremoved++; 
    }
  }
  
  Rcpp::Rcout << "Set " << nremoved << " CADJ edges inactive." << std::endl;
  
  if(issue_missing_warning) {
    std::string warning_msg = "The following edge ids do not reference rows of TRN$CADJ_EL (ignored):\n" + missing_list; 
    Rcpp::warning(warning_msg);
  } 
  
  return;
}

inline void TRN::prune_CADJ_vertex_id(const arma::uvec& vertex_ids) {
  // vertex_ids is assumed to contain (1-based) vertex ids
  if(!this->isinit) Rcpp::stop("Must call $initialize_TRN before calling $prune_CADJ_vertex_id");

  bool issue_missing_warning = false;
  std::string missing_list = ""; 
  
  // Loop over each given vertex_id
  unsigned int nremoved = 0; 
  for(unsigned int i=0; i<vertex_ids.n_elem; ++i) {
    // Find rows of EL containing this vertex
    arma::uvec idxlist = arma::find(this->CADJ_EL.col(0)==vertex_ids[i] || this->CADJ_EL.col(1)==vertex_ids[i]);
    
    if(idxlist.n_elem == 0) {
      issue_missing_warning = true;
      missing_list += std::to_string(vertex_ids[i]) + ", "; 
    } else {
      this->CADJ_active.elem(idxlist).fill(0);  
      nremoved += idxlist.n_elem; 
    }

  }
  
  Rcpp::Rcout << "Set " << nremoved << " CADJ edges inactive." << std::endl;

  if(issue_missing_warning) {
    std::string warning_msg = "The following vertex ids do not appear in TRN$CADJ_EL (ignored):\n" + missing_list; 
    Rcpp::warning(warning_msg);
  } 

  return;
}

inline void TRN::prune_CADJ_edge_weight(double min_weight) {
  if(!this->isinit) Rcpp::stop("Must call $initialize_TRN before calling $prune_CADJ_edge_weight");

  arma::uvec idxlist = arma::find(this->CADJ < min_weight);
  if(idxlist.n_elem == 0) return;
  this->CADJ_active.elem(idxlist).zeros();
  
  Rcpp::Rcout << "Set " << idxlist.n_elem << " CADJ edges inactive." << std::endl;
}

inline void TRN::prune_CADJ_lrank(unsigned int max_lrank) {
  if(!this->isinit) Rcpp::stop("Must call $initialize_TRN before calling $prune_CADJ_lrank");

  arma::uvec idxlist = arma::find(this->CADJ_lrank > max_lrank);
  if(idxlist.n_elem == 0) return;
  this->CADJ_active.elem(idxlist).zeros();
  
  Rcpp::Rcout << "Set " << idxlist.n_elem << " CADJ edges inactive." << std::endl;
}

inline void TRN::prune_CADJ_grank(unsigned int max_grank) {
  if(!this->isinit) Rcpp::stop("Must call $initialize_TRN before calling $prune_CADJ_grank");

  arma::uvec idxlist = arma::find(this->CADJ_grank > max_grank);
  if(idxlist.n_elem == 0) return;
  this->CADJ_active.elem(idxlist).zeros();
  
  Rcpp::Rcout << "Set " << idxlist.n_elem << " CADJ edges inactive." << std::endl;
}

inline void TRN::prune_CADJ_length(unsigned int max_length) {
  if(!this->isinit) Rcpp::stop("Must call $initialize_TRN before calling $prune_CADJ_length");

  arma::uvec idxlist = arma::find(this->CADJ_fwdflen > max_length);
  if(idxlist.n_elem == 0) return;
  this->CADJ_active.elem(idxlist).zeros();
  
  Rcpp::Rcout << "Set " << idxlist.n_elem << " CADJ edges inactive." << std::endl;
}

inline void TRN::prune_CADJ_wdeg(double min_weight) {
  if(!this->isinit) Rcpp::stop("Must call $initialize_TRN before calling $prune_CADJ_wdeg");

  arma::uvec idxlist = arma::find(this->CADJ_wdeg < min_weight && this->CADJ_wdeg > 0);
  if(idxlist.n_elem == 0) return;
  
  idxlist += 1; 
  this->prune_CADJ_vertex_id(idxlist); 
}

inline arma::mat TRN::get_CADJ() {
  //if(!this->isinit) Rcpp::stop("Must call $initialize_TRN before calling $get_CADJ");

  arma::mat ADJ(this->nV, this->nV);
  ADJ.zeros();

  // Add all weights
  arma::uvec is_active = arma::find(this->CADJ_active);
  ADJ.elem(arma::sub2ind(arma::size(ADJ), this->CADJ_EL.rows(is_active).t()-1)) = this->CADJ.elem(is_active);

  return ADJ;
}


inline void TRN::prune_CONN_edge_list(const arma::umat& edge_list) {
  // edge_list is assumed to contain (v1,v2) pairs in rows, which will be matched and inactivated
  if(!this->isinit) Rcpp::stop("Must call $initialize_TRN before calling $prune_CONN_edge_list");
  
  if(edge_list.n_cols != 2) Rcpp::stop("Input edge_list of wrong format. Assumed to be a matrix with (v1,v2) pairs in rows");
  
  bool issue_missing_warning = false;
  std::string missing_list = ""; 
  
  // Loop over given edgelist
  unsigned int nremoved = 0; 
  for(unsigned int i=0; i<edge_list.n_rows; ++i) {
    
    // Find this edge in the TRG EL. If it doesn't exist, issue warning and skip
    arma::uvec matchidx = arma::find(this->CONN_EL.col(0)==edge_list(i,0) && this->CONN_EL.col(1)==edge_list(i,1));
    if(matchidx.n_elem == 0) {
      issue_missing_warning = true;
      missing_list += "(" + std::to_string(edge_list(i,0)) + "," + std::to_string(edge_list(i,1)) + ") "; 
      continue;
    } else {
      // Otherwise inactive the edge
      this->CONN_active[matchidx[0]] = 0;
      nremoved++;
    }
    
  }
  
  Rcpp::Rcout << "Set " << nremoved << " CONN edges inactive." << std::endl;
  
  if(issue_missing_warning) {
    std::string warning_msg = "The following edges were not found in TRN$CONN_EL (ignored):\n" + missing_list; 
    Rcpp::warning(warning_msg);
  } 
  return;
}

inline void TRN::prune_CONN_edge_id(const arma::uvec& edge_ids) {
  // edge_ids is assumed to contain (1-based) indices to the rows of EL to inactive
  if(!this->isinit) Rcpp::stop("Must call $initialize_TRN before calling $prune_CONN_edge_id");
  
  bool issue_missing_warning = false;
  std::string missing_list = ""; 
  
  // Loop over given edge ids
  unsigned int nremoved = 0; 
  for(unsigned int i=0; i<edge_ids.n_elem; ++i) {
    
    // If edge_id in range, set it to 0. 
    // Otherwise, log it as missing 
    if(edge_ids[i] > this->CONN_nE || edge_ids[i] == 0) {
      issue_missing_warning = true;
      missing_list += std::to_string(edge_ids[i]) + ", ";
    } else {
      this->CONN_active[edge_ids[i]-1] = 0; 
    }
  }
  
  Rcpp::Rcout << "Set " << nremoved << " CONN edges inactive." << std::endl;
  
  if(issue_missing_warning) {
    std::string warning_msg = "The following edge ids do not reference rows of TRN$CONN_EL (ignored):\n" + missing_list; 
    Rcpp::warning(warning_msg);
  } 
  
  return;
}

inline void TRN::prune_CONN_vertex_id(const arma::uvec& vertex_ids) {
  // vertex_ids is assumed to contain (1-based) vertex ids
  if(!this->isinit) Rcpp::stop("Must call $initialize_TRN before calling $prune_CONN_vertex_id");
  
  bool issue_missing_warning = false;
  std::string missing_list = ""; 
  
  // Loop over each given vertex_id
  unsigned int nremoved = 0; 
  for(unsigned int i=0; i<vertex_ids.n_elem; ++i) {
    // Find rows of EL containing this vertex
    arma::uvec idxlist = arma::find(this->CONN_EL.col(0)==vertex_ids[i] || this->CONN_EL.col(1)==vertex_ids[i]);
    
    if(idxlist.n_elem == 0) {
      issue_missing_warning = true;
      missing_list += std::to_string(vertex_ids[i]) + ", "; 
    } else {
      this->CONN_active.elem(idxlist).fill(0);  
    }
    
  }
  
  Rcpp::Rcout << "Set " << nremoved << " CONN edges inactive." << std::endl;
  
  if(issue_missing_warning) {
    std::string warning_msg = "The following vertex ids do not appear in TRN$CONN_EL (ignored):\n" + missing_list; 
    Rcpp::warning(warning_msg);
  } 
  
  return;
}

inline void TRN::prune_CONN_edge_weight(double min_weight) {
  if(!this->isinit) Rcpp::stop("Must call $initialize_TRN before calling $prune_CONN_edge_weight");
  
  arma::uvec idxlist = arma::find(this->CONN < min_weight);
  if(idxlist.n_elem == 0) return;
  this->CONN_active.elem(idxlist).zeros();
  
  Rcpp::Rcout << "Set " << idxlist.n_elem << " CONN edges inactive." << std::endl;
}

inline void TRN::prune_CONN_lrank(unsigned int max_lrank) {
  if(!this->isinit) Rcpp::stop("Must call $initialize_TRN before calling $prune_CONN_lrank");
  
  arma::uvec idxlist = arma::find(this->CONN_lrank > max_lrank);
  if(idxlist.n_elem == 0) return;
  this->CONN_active.elem(idxlist).zeros();
  
  Rcpp::Rcout << "Set " << idxlist.n_elem << " CONN edges inactive." << std::endl;
}

inline void TRN::prune_CONN_grank(unsigned int max_grank) {
  if(!this->isinit) Rcpp::stop("Must call $initialize_TRN before calling $prune_CONN_grank");
  
  arma::uvec idxlist = arma::find(this->CONN_grank > max_grank);
  if(idxlist.n_elem == 0) return;
  this->CONN_active.elem(idxlist).zeros();
  
  Rcpp::Rcout << "Set " << idxlist.n_elem << " CONN edges inactive." << std::endl;
}

inline void TRN::prune_CONN_length(unsigned int max_length) {
  if(!this->isinit) Rcpp::stop("Must call $initialize_TRN before calling $prune_CONN_length");
  
  arma::uvec idxlist = arma::find(this->CONN_fwdflen > max_length);
  if(idxlist.n_elem == 0) return;
  this->CONN_active.elem(idxlist).zeros();
  
  Rcpp::Rcout << "Set " << idxlist.n_elem << " CONN edges inactive." << std::endl;
}

inline void TRN::prune_CONN_wdeg(double min_weight) {
  if(!this->isinit) Rcpp::stop("Must call $initialize_TRN before calling $prune_CONN_wdeg");
  
  arma::uvec idxlist = arma::find(this->CONN_wdeg < min_weight && this->CONN_wdeg > 0);
  if(idxlist.n_elem == 0) return;
  
  idxlist += 1; 
  this->prune_CONN_vertex_id(idxlist); 
}

inline void TRN::prune_CONN_CADJ() {
  // Find inactive CADJ edges
  arma::uvec inactive_CADJ = arma::find(this->CADJ_active == 0); 
  
  for(unsigned int i=0; i<inactive_CADJ.n_elem; ++i) {
    unsigned int iidx = this->CADJ_EL(inactive_CADJ[i], 0); 
    unsigned int jidx = this->CADJ_EL(inactive_CADJ[i], 1); 
    arma::uvec CONN_EL_idx = arma::find(this->CONN_EL.col(0) == iidx && this->CONN_EL.col(1) == jidx);
    if(CONN_EL_idx.n_elem > 0) {
      this->CONN_active[CONN_EL_idx[0]] = 0; 
    }
    
    // If (jidx,iidx) is an ACTIVE CADJ edge, leave the reciprocal CONN edge 
    // Otherwise, remove it 
    bool recip_CADJ_active = arma::any(this->CADJ_EL.col(0) == jidx && this->CADJ_EL.col(1) == iidx && this->CADJ_active == 1);
    
    if(!recip_CADJ_active) {
      CONN_EL_idx = arma::find(this->CONN_EL.col(0) == jidx && this->CONN_EL.col(1) == iidx);
      if(CONN_EL_idx.n_elem > 0) {
        this->CONN_active[CONN_EL_idx[0]] = 0; 
      }
    }
  }
  
}

inline arma::mat TRN::get_CONN() {
  //if(!this->isinit) Rcpp::stop("Must call $initialize_TRN before calling $get_CONN");
  
  arma::mat ADJ(this->nV, this->nV);
  ADJ.zeros();
  
  // Add all weights
  arma::uvec is_active = arma::find(this->CONN_active);
  ADJ.elem(arma::sub2ind(arma::size(ADJ), this->CONN_EL.rows(is_active).t()-1)) = this->CONN.elem(is_active);
  
  return ADJ;
}

inline arma::mat TRN::get_OTADJ() {
  //if(!this->isinit) Rcpp::stop("Must call $initialize_TRN before calling $get_outputADJ");

  arma::mat ADJ(this->nV, this->nV);
  ADJ.zeros();

  // Add all weights
  ADJ.elem(arma::sub2ind(arma::size(ADJ), this->OT_EL.t()-1)).ones();

  return ADJ;
}
 
inline arma::uvec TRN::CADJ_active_verts() {
  if(!this->isinit) Rcpp::stop("Must call $initialize_TRN before calling $CADJ_active_verts");

  arma::mat tmp = this->get_CADJ();
  arma::vec rowsums = arma::sum(tmp, 1);
  arma::rowvec colsums = arma::sum(tmp, 0);

  arma::uvec out(this->nV); out.zeros();
  out.elem(arma::find(rowsums > 0)).ones();
  out.elem(arma::find(colsums > 0)).ones();

  return arma::find(out)+1;
}

inline arma::uvec TRN::CONN_active_verts() {
  if(!this->isinit) Rcpp::stop("Must call $initialize_TRN before calling $CONN_active_verts");
  
  arma::mat tmp = this->get_CONN();
  arma::vec rowsums = arma::sum(tmp, 1);
  arma::rowvec colsums = arma::sum(tmp, 0);
  
  arma::uvec out(this->nV); out.zeros();
  out.elem(arma::find(rowsums > 0)).ones();
  out.elem(arma::find(colsums > 0)).ones();
  
  return arma::find(out)+1;
}


// ***** Topology Preserving Measures
inline void TRN::calc_TopoMeasures(const arma::mat& inV, const arma::mat& outV) {
  if(!this->isinit) Rcpp::stop("Must call $initialize_TRN before calling $calc_TopoMeasures");

  Rcpp::Rcout << "Computing Toplogy Preservation Measures:" << std::endl;

  Rcpp::Rcout << "++ Topographic Product ... ";
  this->TopoProd = TRN_UTILS_TOPOPRES::Topographic_Product(inV, outV, this->parallel);
  Rcpp::Rcout << "done" << std::endl;

  Rcpp::Rcout << "++ CADJ Topographic Functions ... ";
  arma::umat inputADJ = arma::conv_to<arma::umat>::from(this->get_CADJ());
  arma::umat outputADJ = arma::conv_to<arma::umat>::from(this->get_OTADJ());
  this->CADJ_TopoFxns = TRN_UTILS_TOPOPRES::Topographic_Functions(inputADJ, outputADJ, this->parallel);
  Rcpp::Rcout << "done" << std::endl;
  
  Rcpp::Rcout << "++ CONN Topographic Functions ... ";
  inputADJ = arma::conv_to<arma::umat>::from(this->get_CONN());
  this->CONN_TopoFxns = TRN_UTILS_TOPOPRES::Topographic_Functions(inputADJ, outputADJ, this->parallel);
  Rcpp::Rcout << "done" << std::endl;
  

  this->isset_TPM = true;

  return;
}


// ***** DM-Prune 
inline void TRN::set_BootSig(const arma::mat& BSADJ) {
  if(!this->isinit) Rcpp::stop("Must call $initialize_TRN before calling $set_BootSig");
  if(BSADJ.n_rows!=this->nV || BSADJ.n_cols!= this->nV) Rcpp::stop("BSADJ must be square, with nrows = ncols = num. graph vertices.");
  if(arma::any(arma::vectorise(BSADJ) < 0.0) || arma::any(arma::vectorise(BSADJ) > 1.0)) Rcpp::stop("BSADJ must have values in [0,1]");
  
  this->CADJ_BootSig = BSADJ.elem(arma::sub2ind(arma::size(BSADJ), this->CADJ_EL.t()-1));
  
  arma::mat BSADJsym = (BSADJ + BSADJ.t()) / 2; 
  this->CONN_BootSig = BSADJsym.elem(arma::sub2ind(arma::size(BSADJ), this->CONN_EL.t()-1));
}

inline void TRN::calc_DMPrune_LambdaPath(const arma::mat& priorADJ) {
  if(!this->isinit) Rcpp::stop("Must call $initialize_TRN before calling $calc_DMPrune_LambdaPath");
  
  Rcpp::Rcout << "Building DM-Prune Lambda Path:" << std::endl; 
  Rcpp::Rcout << "----------------------------------------------------------------" << std::endl;
  
  // Set the gamma values 
  DMP_gamma_G = 1.0, DMP_gamma_L = 1.0, DMP_gamma_N = 1.0; DMP_gamma_S = 1.0; 
  Rcpp::Rcout << "++ gamma values: G = " << DMP_gamma_G << ", N = " << DMP_gamma_N << ", L = " << DMP_gamma_L << ", S = " << DMP_gamma_S << std::endl; 
  if(this->CADJ_BootSig.n_elem == 0) {
    Rcpp::Rcout << "No bootstrap edge significances found, skipping.\n   Call $set_BootSig to incorporate." << std::endl;   
  }
  
  // Rcpp::Rcout << "----------------------------------------------------------------" << std::endl;
  // // Store the uniform probabilities 
  // Rcpp::Rcout << "Storing Uniform probabilities ... ";
  // this->set_unif_probs(UPADJ); 
  // Rcpp::Rcout << "done" << std::endl; 
  // 
  // Rcpp::Rcout << "Computing Uniform prior ... "; 
  // this->DMP_CADJ_prior = arma::accu(this->CADJ) * this->CADJ_unif_probs; 
  // Rcpp::Rcout << "done" << std::endl; 
  
  Rcpp::Rcout << "Checking & storing prior ... ";
  if(priorADJ.n_rows!=this->nV || priorADJ.n_cols!=this->nV) 
    Rcpp::stop("priorADJ must be square, with nrows = ncols = num. graph vertices.");
  
  arma::uvec extract_these = arma::sub2ind(arma::size(priorADJ), this->CADJ_EL.t()-1);
  if(!arma::all(priorADJ.elem(extract_these) > 0))
    Rcpp::stop("priorADJ must be > 0 wherever CADJ is > 0.");
  this->DMP_prior = priorADJ.elem(extract_these); 
  Rcpp::Rcout << "done" << std::endl; 
  
  // Visualize prior 
  Rcpp::Environment base = Rcpp::Environment::namespace_env("TopoRNet");
  Rcpp::Function visprior = base["vis_DMPrune_prior"];
  visprior(*this); 
  //saveRDS(Rcpp::wrap(TRNList), Rcpp::Named("file", rdsfile));
  
  
  // ***** Computing mu scores
  Rcpp::Rcout << "----------------------------------------------------------------" << std::endl;
  Rcpp::Rcout << "Computing mu scores:" << std::endl; 
  
  // Global score 
  Rcpp::Rcout << "++ mu_G (global) ... ";
  this->DMP_mu_G = this->CADJ / this->CADJ.max(); 
  Rcpp::Rcout << "done" << std::endl; 
  
  // Local score 
  Rcpp::Rcout << "++ mu_N (local) ... ";
  arma::mat tmpCADJ = arma::zeros<arma::mat>(this->nV, this->nV);
  tmpCADJ.elem(arma::sub2ind(arma::size(tmpCADJ), this->CADJ_EL.t()-1)) = this->CADJ; 
  arma::vec local_maxrank = arma::max(tmpCADJ, 1); 
  this->DMP_mu_N = this->CADJ / local_maxrank.elem(this->CADJ_EL.col(0)-1); 
  Rcpp::Rcout << "done" << std::endl; 
  
  // Length score 
  Rcpp::Rcout << "++ mu_L (length) ... ";
  this->DMP_mu_L = TRN_UTILS_DMPRUNE::DMP_length_score(this->CADJ_fwdflen, this->OT_max_dist, this->CADJ_TP_radius); 
  Rcpp::Rcout << "done" << std::endl; 
  
  // Bootstrap score 
  if(this->CADJ_BootSig.n_elem == 0) {
    this->DMP_mu_S = arma::zeros<arma::vec>(this->CADJ_nE); 
    DMP_gamma_S = 0.0; 
  } else {
    Rcpp::Rcout << "++ mu_S (bootstrap significance) ... ";
    this->DMP_mu_S = this->CADJ_BootSig;
    Rcpp::Rcout << "done" << std::endl; 
  }
  
  
  // *** Compute aggregate score 
  Rcpp::Rcout << "++ aggregate mu ... ";
  this->DMP_mu = arma::pow(DMP_mu_G, DMP_gamma_G) % arma::pow(DMP_mu_N, DMP_gamma_N) % arma::pow(DMP_mu_L, DMP_gamma_L) % arma::pow(DMP_mu_S, DMP_gamma_S); 
  this->DMP_mu = arma::pow(this->DMP_mu, 1.0 / (DMP_gamma_G + DMP_gamma_N + DMP_gamma_L + DMP_gamma_S)); 
  Rcpp::Rcout << "done" << std::endl; 
  
  // *** Rank 
  Rcpp::Rcout << "++ ranking mu ... "; 
  arma::vec unq_mu = arma::sort(arma::unique(this->DMP_mu)); 
  std::map<double, unsigned int> unq_mu_lookup; 
  for(unsigned int i=0; i<unq_mu.n_elem; ++i) {
    unq_mu_lookup[unq_mu[i]] = i+1; 
  }
  this->DMP_mu_rank.set_size(this->CADJ_nE); 
  for(unsigned int i=0; i<this->CADJ_nE; ++i) {
    DMP_mu_rank[i] = unq_mu_lookup[DMP_mu[i]];
  }
  Rcpp::Rcout << "done" << std::endl; 
  
  
  
  // ***** Computing Lambda
  Rcpp::Rcout << "----------------------------------------------------------------" << std::endl;
  Rcpp::Rcout << "Computing Lambda Path ... "; 
  TRN_UTILS_DMPRUNE::DMPrune_Lambda_worker Lambda_wkr(this->DMP_mu_rank, this->CADJ, this->DMP_prior); 
  if(this->parallel) Lambda_wkr.calc_parallel(); else Lambda_wkr.calc_serial(); 
  this->DMP_pruneStep = Lambda_wkr.pruneStep; 
  this->DMP_Lambda = Lambda_wkr.Lambda; 
  Rcpp::Rcout << "done" << std::endl; 
  
  this->isset_DMPrune = true; 
  
}

inline void TRN::calc_DMPrune_LambdaPath_loo(const arma::mat& priorADJ) {
  if(!this->isinit) Rcpp::stop("Must call $initialize_TRN before calling $calc_DMPrune_LambdaPath_loo");
  
  Rcpp::Rcout << "Building DM-Prune Lambda Path:" << std::endl; 
  Rcpp::Rcout << "----------------------------------------------------------------" << std::endl;
  
  Rcpp::Rcout << "Checking & storing prior ... ";
  if(priorADJ.n_rows!=this->nV || priorADJ.n_cols!=this->nV) 
    Rcpp::stop("priorADJ must be square, with nrows = ncols = num. graph vertices.");
  
  arma::uvec extract_these = arma::sub2ind(arma::size(priorADJ), this->CADJ_EL.t()-1);
  if(!arma::all(priorADJ.elem(extract_these) > 0))
    Rcpp::stop("priorADJ must be > 0 wherever CADJ is > 0.");
  this->DMP_prior = priorADJ.elem(extract_these); 
  Rcpp::Rcout << "done" << std::endl; 
  
  // ***** Computing Lambda
  Rcpp::Rcout << "----------------------------------------------------------------" << std::endl;
  Rcpp::Rcout << "Computing Leave-One-Out Lambda Path ... " << std::endl; 
  TRN_UTILS_DMPRUNE::DMPrune_Lambda_loo_worker wkr(this->CADJ, this->DMP_prior); 
  wkr.calc_parallel(); 
  this->DMP_loo_llk = wkr.loo_llk; 
  this->DMP_loo_order = wkr.loo_order; 
  
  this->isset_DMPrune = true; 
  
}

inline void TRN::calc_DMPrune_LambdaPath_error() {
  if(!this->isinit) Rcpp::stop("Must call $initialize_TRN before calling $calc_DMPrune_LambdaPath_error");
  if(!this->isset_DMPrune) Rcpp::stop("Must call $calc_DMPrune_LambdaPath before calling $calc_DMPrune_LambdaPath_error");
  
  TRN_UTILS_DMPRUNE::DMPrune_Lambda_errorband_worker wkr(this->DMP_mu_rank, this->CADJ, this->DMP_prior, this->DMP_Lambda_error_nperms, this->DMP_Lambda_error_qlo, this->DMP_Lambda_error_qhi);
  if(this->parallel) wkr.calc_parallel(); else wkr.calc_serial(); 
  
  this->DMP_Lambda_error.set_size(this->DMP_Lambda.n_elem, 2); 
  this->DMP_Lambda_error.col(0) = wkr.lo; 
  this->DMP_Lambda_error.col(1) = wkr.hi; 
  
  return; 
}

inline void TRN::DMPrune_CADJ_step(unsigned int min_step) {
  if(!this->isinit) Rcpp::stop("Must call $initialize_TRN before calling $DMPrune_CADJ_step");
  if(!this->isset_DMPrune) Rcpp::stop("Must call $calc_DMPrune_LambdaPath first");
  
  // Check step 
  if(min_step > this->DMP_pruneStep.max()) Rcpp::stop("min_step must be <= max(DMP_pruneStep)");
  
  // Set all edges whose mu_rank <  min_step to INACTIVE 
  arma::uvec remove_these = arma::find(this->DMP_mu_rank < min_step);
  this->CADJ_active.elem(remove_these).zeros(); 
  
  Rcpp::Rcout << "Set " << remove_these.n_elem << " CADJ edges inactive." << std::endl; 
}


// // ***** Saving & Loading ***** 
inline void TRN::save(std::string rdsfile) {

  std::string rdsend = ".trn";
  if(rdsend.size() > rdsfile.size() || !std::equal(rdsend.rbegin(), rdsend.rend(), rdsfile.rbegin()))
    Rcpp::stop("Save file string must end with .trn");

  Rcpp::List TRNList = this->as_list(); 

  // Save 
  Rcpp::Environment base = Rcpp::Environment::namespace_env("base");
  Rcpp::Function saveRDS = base["saveRDS"];
  saveRDS(Rcpp::wrap(TRNList), Rcpp::Named("file", rdsfile));

  return;
}

inline void TRN::load(std::string rdsfile) {

  // Read file 
  Rcpp::Environment base = Rcpp::Environment::namespace_env("base");
  Rcpp::Function readRDS = base["readRDS"];
  Rcpp::List TRNList = readRDS(Rcpp::Named("file", rdsfile));
  
  this->load_list(TRNList); 

  return;

}

inline Rcpp::List TRN::as_list() {
  
  Rcpp::List TRNList;
  
  TRNList["parallel"] = this->parallel;
  
  // Graphs stored internally
  TRNList["nV"] = this->nV;
  TRNList["CADJ_nE"] = this->CADJ_nE; 
  TRNList["CONN_nE"] = this->CONN_nE; 
  TRNList["OT_nE"] = this->OT_nE; 
  TRNList["CADJ_EL"] = this->CADJ_EL;
  TRNList["CONN_EL"] = this->CONN_EL;
  TRNList["OT_EL"] = this->OT_EL;
  TRNList["CADJ"] = this->CADJ;
  TRNList["CONN"] = this->CONN;
  TRNList["CADJ_deg"] = this->CADJ_deg;
  TRNList["CONN_deg"] = this->CONN_deg;
  TRNList["OT_deg"] = this->OT_deg;
  TRNList["CADJ_wdeg"] = this->CADJ_wdeg;
  TRNList["CONN_wdeg"] = this->CONN_wdeg;
  TRNList["CADJ_fwdflen"] = this->CADJ_fwdflen;
  TRNList["CONN_fwdflen"] = this->CONN_fwdflen;
  TRNList["CADJ_bwdflen"] = this->CADJ_bwdflen;
  TRNList["CONN_bwdflen"] = this->CONN_bwdflen;
  TRNList["CADJ_TP_radius"] = this->CADJ_TP_radius;
  TRNList["CONN_TP_radius"] = this->CONN_TP_radius;
  
  TRNList["OT_geodesic_dist"] = this->OT_geodesic_dist; 
  TRNList["OT_max_dist"] = this->OT_max_dist;
  TRNList["OT_nhb_sizes"] = this->OT_nhb_sizes; 
  
  
  // CONNvis 
  TRNList["CADJ_lrank"] = this->CADJ_lrank;
  TRNList["CADJ_grank"] = this->CADJ_grank;
  TRNList["CADJ_CVstats_lrank"] = this->get_CADJvis_stats("lrank");
  TRNList["CADJ_CVstats_grank"] = this->get_CADJvis_stats("grank");
  TRNList["CADJ_CVstats_length"] = this->get_CADJvis_stats("length");
  
  TRNList["CONN_lrank"] = this->CONN_lrank;
  TRNList["CONN_grank"] = this->CONN_grank;
  TRNList["CONN_CVstats_lrank"] = this->get_CONNvis_stats("lrank");
  TRNList["CONN_CVstats_grank"] = this->get_CONNvis_stats("grank");
  TRNList["CONN_CVstats_length"] = this->get_CONNvis_stats("length");
  
  // Pruning 
  TRNList["CADJ_active"] = this->CADJ_active;
  TRNList["CONN_active"] = this->CONN_active;
  
  // Topology Preservation Measures 
  TRNList["TopoProd"] = this->TopoProd;
  TRNList["CADJ_TopoFxns"] = this->CADJ_TopoFxns;
  TRNList["CONN_TopoFxns"] = this->CONN_TopoFxns;
  
  
  // DM-Prune 
  TRNList["CADJ_BootSig"] = this->CADJ_BootSig;
  TRNList["CONN_BootSig"] = this->CONN_BootSig;
  TRNList["DMP_prior"] = this->DMP_prior;
  
  TRNList["DMP_gamma_G"] = this->DMP_gamma_G;
  TRNList["DMP_gamma_N"] = this->DMP_gamma_N;
  TRNList["DMP_gamma_L"] = this->DMP_gamma_L;
  TRNList["DMP_gamma_S"] = this->DMP_gamma_S;
  
  TRNList["DMP_mu_G"] = this->DMP_mu_G;
  TRNList["DMP_mu_N"] = this->DMP_mu_N;
  TRNList["DMP_mu_L"] = this->DMP_mu_L;
  TRNList["DMP_mu_S"] = this->DMP_mu_S;
  TRNList["DMP_mu"] = this->DMP_mu;
  TRNList["DMP_mu_rank"] = this->DMP_mu_rank;
  
  TRNList["DMP_pruneStep"] = this->DMP_pruneStep;
  TRNList["DMP_Lambda"] = this->DMP_Lambda;
  
  TRNList["DMP_Lambda_error_nperms"] = this->DMP_Lambda_error_nperms;
  TRNList["DMP_Lambda_error_qlo"] = this->DMP_Lambda_error_qlo;
  TRNList["DMP_Lambda_error_qhi"] = this->DMP_Lambda_error_qhi;
  TRNList["DMP_Lambda_error"] = this->DMP_Lambda_error;
  
  TRNList["DMP_loo_order"] = this->DMP_loo_order;
  TRNList["DMP_loo_llk"] = this->DMP_loo_llk;
  
  // Flags 
  TRNList["isinit"] = this->isinit; 
  TRNList["isset_TPM"] = this->isset_TPM;
  TRNList["isset_DMPrune"] = this->isset_DMPrune;
  
  return TRNList; 
}

inline void TRN::load_list(Rcpp::List TRNList) {
  
  // Graph storage 
  this->parallel = TRNList["parallel"];
  this->nV = TRNList["nV"];
  this->CADJ_nE = TRNList["CADJ_nE"];
  this->CONN_nE = TRNList["CONN_nE"];
  this->OT_nE = TRNList["OT_nE"];
  this->CADJ_EL = Rcpp::as<arma::umat>(TRNList["CADJ_EL"]);
  this->CONN_EL = Rcpp::as<arma::umat>(TRNList["CONN_EL"]);
  this->OT_EL = Rcpp::as<arma::umat>(TRNList["OT_EL"]);
  this->CADJ = Rcpp::as<arma::vec>(TRNList["CADJ"]);
  this->CONN = Rcpp::as<arma::vec>(TRNList["CONN"]);
  this->CADJ_deg = Rcpp::as<arma::uvec>(TRNList["CADJ_deg"]);
  this->CONN_deg = Rcpp::as<arma::uvec>(TRNList["CONN_deg"]);
  this->OT_deg = Rcpp::as<arma::uvec>(TRNList["OT_deg"]);
  this->CADJ_wdeg = Rcpp::as<arma::vec>(TRNList["CADJ_wdeg"]);
  this->CONN_wdeg = Rcpp::as<arma::vec>(TRNList["CONN_wdeg"]);
  this->CADJ_fwdflen = Rcpp::as<arma::uvec>(TRNList["CADJ_fwdflen"]);
  this->CONN_fwdflen = Rcpp::as<arma::uvec>(TRNList["CONN_fwdflen"]);
  this->CADJ_bwdflen = Rcpp::as<arma::uvec>(TRNList["CADJ_bwdflen"]);
  this->CONN_bwdflen = Rcpp::as<arma::uvec>(TRNList["CONN_bwdflen"]);
  this->CADJ_TP_radius = TRNList["CADJ_TP_radius"];
  this->CONN_TP_radius = TRNList["CONN_TP_radius"];
  
  this->OT_geodesic_dist = Rcpp::as<arma::umat>(TRNList["OT_geodesic_dist"]); 
  this->OT_max_dist = TRNList["OT_max_dist"];
  this->OT_nhb_sizes = Rcpp::as<arma::umat>(TRNList["OT_nhb_sizes"]);
  
  // CONNvis 
  this->CADJ_lrank = Rcpp::as<arma::uvec>(TRNList["CADJ_lrank"]);
  this->CADJ_grank = Rcpp::as<arma::uvec>(TRNList["CADJ_grank"]);
  
  Rcpp::DataFrame tmp_CADJ_CVstats_lrank = Rcpp::as<Rcpp::DataFrame>(TRNList["CADJ_CVstats_lrank"]);
  this->CADJ_CVstats_lrank.group_id = Rcpp::as<arma::uvec>(tmp_CADJ_CVstats_lrank["lrank"]);
  this->CADJ_CVstats_lrank.count = Rcpp::as<arma::uvec>(tmp_CADJ_CVstats_lrank["count"]);
  this->CADJ_CVstats_lrank.pct = Rcpp::as<arma::vec>(tmp_CADJ_CVstats_lrank["pct"]);
  this->CADJ_CVstats_lrank.cumpct = Rcpp::as<arma::vec>(tmp_CADJ_CVstats_lrank["cumpct"]);
  this->CADJ_CVstats_lrank.mean = Rcpp::as<arma::vec>(tmp_CADJ_CVstats_lrank["mean"]);
  this->CADJ_CVstats_lrank.sum = this->CADJ_CVstats_lrank.mean;
  this->CADJ_CVstats_lrank.sum.each_col() %= arma::conv_to<arma::vec>::from(this->CADJ_CVstats_lrank.count);
  this->CADJ_CVstats_lrank.sd = Rcpp::as<arma::vec>(tmp_CADJ_CVstats_lrank["sd"]);
  this->CADJ_CVstats_lrank.q0 = Rcpp::as<arma::vec>(tmp_CADJ_CVstats_lrank["q0"]);
  this->CADJ_CVstats_lrank.q25 = Rcpp::as<arma::vec>(tmp_CADJ_CVstats_lrank["q25"]);
  this->CADJ_CVstats_lrank.q50 = Rcpp::as<arma::vec>(tmp_CADJ_CVstats_lrank["q50"]);
  this->CADJ_CVstats_lrank.q75 = Rcpp::as<arma::vec>(tmp_CADJ_CVstats_lrank["q75"]);
  this->CADJ_CVstats_lrank.q100 = Rcpp::as<arma::vec>(tmp_CADJ_CVstats_lrank["q100"]);
  this->CADJ_CVstats_lrank.ngroups = this->CADJ_CVstats_lrank.group_id.n_elem;
  
  Rcpp::DataFrame tmp_CADJ_CVstats_grank = Rcpp::as<Rcpp::DataFrame>(TRNList["CADJ_CVstats_grank"]);
  this->CADJ_CVstats_grank.group_id = Rcpp::as<arma::uvec>(tmp_CADJ_CVstats_grank["grank"]);
  this->CADJ_CVstats_grank.count = Rcpp::as<arma::uvec>(tmp_CADJ_CVstats_grank["count"]);
  this->CADJ_CVstats_grank.pct = Rcpp::as<arma::vec>(tmp_CADJ_CVstats_grank["pct"]);
  this->CADJ_CVstats_grank.cumpct = Rcpp::as<arma::vec>(tmp_CADJ_CVstats_grank["cumpct"]);
  this->CADJ_CVstats_grank.mean = Rcpp::as<arma::vec>(tmp_CADJ_CVstats_grank["mean"]);
  this->CADJ_CVstats_grank.sum = this->CADJ_CVstats_grank.mean;
  this->CADJ_CVstats_grank.sum.each_col() %= arma::conv_to<arma::vec>::from(this->CADJ_CVstats_grank.count);
  this->CADJ_CVstats_grank.sd = Rcpp::as<arma::vec>(tmp_CADJ_CVstats_grank["sd"]);
  this->CADJ_CVstats_grank.q0 = Rcpp::as<arma::vec>(tmp_CADJ_CVstats_grank["q0"]);
  this->CADJ_CVstats_grank.q25 = Rcpp::as<arma::vec>(tmp_CADJ_CVstats_grank["q25"]);
  this->CADJ_CVstats_grank.q50 = Rcpp::as<arma::vec>(tmp_CADJ_CVstats_grank["q50"]);
  this->CADJ_CVstats_grank.q75 = Rcpp::as<arma::vec>(tmp_CADJ_CVstats_grank["q75"]);
  this->CADJ_CVstats_grank.q100 = Rcpp::as<arma::vec>(tmp_CADJ_CVstats_grank["q100"]);
  this->CADJ_CVstats_grank.ngroups = this->CADJ_CVstats_grank.group_id.n_elem;
  
  Rcpp::DataFrame tmp_CADJ_CVstats_length = Rcpp::as<Rcpp::DataFrame>(TRNList["CADJ_CVstats_length"]);
  this->CADJ_CVstats_length.group_id = Rcpp::as<arma::uvec>(tmp_CADJ_CVstats_length["length"]);
  this->CADJ_CVstats_length.count = Rcpp::as<arma::uvec>(tmp_CADJ_CVstats_length["count"]);
  this->CADJ_CVstats_length.pct = Rcpp::as<arma::vec>(tmp_CADJ_CVstats_length["pct"]);
  this->CADJ_CVstats_length.cumpct = Rcpp::as<arma::vec>(tmp_CADJ_CVstats_length["cumpct"]);
  this->CADJ_CVstats_length.mean = Rcpp::as<arma::vec>(tmp_CADJ_CVstats_length["mean"]);
  this->CADJ_CVstats_length.sum = this->CADJ_CVstats_length.mean;
  this->CADJ_CVstats_length.sum.each_col() %= arma::conv_to<arma::vec>::from(this->CADJ_CVstats_length.count);
  this->CADJ_CVstats_length.sd = Rcpp::as<arma::vec>(tmp_CADJ_CVstats_length["sd"]);
  this->CADJ_CVstats_length.q0 = Rcpp::as<arma::vec>(tmp_CADJ_CVstats_length["q0"]);
  this->CADJ_CVstats_length.q25 = Rcpp::as<arma::vec>(tmp_CADJ_CVstats_length["q25"]);
  this->CADJ_CVstats_length.q50 = Rcpp::as<arma::vec>(tmp_CADJ_CVstats_length["q50"]);
  this->CADJ_CVstats_length.q75 = Rcpp::as<arma::vec>(tmp_CADJ_CVstats_length["q75"]);
  this->CADJ_CVstats_length.q100 = Rcpp::as<arma::vec>(tmp_CADJ_CVstats_length["q100"]);
  this->CADJ_CVstats_length.ngroups = this->CADJ_CVstats_length.group_id.n_elem;
  
  
  this->CONN_lrank = Rcpp::as<arma::uvec>(TRNList["CONN_lrank"]);
  this->CONN_grank = Rcpp::as<arma::uvec>(TRNList["CONN_grank"]);
  
  Rcpp::DataFrame tmp_CONN_CVstats_lrank = Rcpp::as<Rcpp::DataFrame>(TRNList["CONN_CVstats_lrank"]);
  this->CONN_CVstats_lrank.group_id = Rcpp::as<arma::uvec>(tmp_CONN_CVstats_lrank["lrank"]);
  this->CONN_CVstats_lrank.count = Rcpp::as<arma::uvec>(tmp_CONN_CVstats_lrank["count"]);
  this->CONN_CVstats_lrank.pct = Rcpp::as<arma::vec>(tmp_CONN_CVstats_lrank["pct"]);
  this->CONN_CVstats_lrank.cumpct = Rcpp::as<arma::vec>(tmp_CONN_CVstats_lrank["cumpct"]);
  this->CONN_CVstats_lrank.mean = Rcpp::as<arma::vec>(tmp_CONN_CVstats_lrank["mean"]);
  this->CONN_CVstats_lrank.sum = this->CONN_CVstats_lrank.mean;
  this->CONN_CVstats_lrank.sum.each_col() %= arma::conv_to<arma::vec>::from(this->CONN_CVstats_lrank.count);
  this->CONN_CVstats_lrank.sd = Rcpp::as<arma::vec>(tmp_CONN_CVstats_lrank["sd"]);
  this->CONN_CVstats_lrank.q0 = Rcpp::as<arma::vec>(tmp_CONN_CVstats_lrank["q0"]);
  this->CONN_CVstats_lrank.q25 = Rcpp::as<arma::vec>(tmp_CONN_CVstats_lrank["q25"]);
  this->CONN_CVstats_lrank.q50 = Rcpp::as<arma::vec>(tmp_CONN_CVstats_lrank["q50"]);
  this->CONN_CVstats_lrank.q75 = Rcpp::as<arma::vec>(tmp_CONN_CVstats_lrank["q75"]);
  this->CONN_CVstats_lrank.q100 = Rcpp::as<arma::vec>(tmp_CONN_CVstats_lrank["q100"]);
  this->CONN_CVstats_lrank.ngroups = this->CONN_CVstats_lrank.group_id.n_elem;
  
  Rcpp::DataFrame tmp_CONN_CVstats_grank = Rcpp::as<Rcpp::DataFrame>(TRNList["CONN_CVstats_grank"]);
  this->CONN_CVstats_grank.group_id = Rcpp::as<arma::uvec>(tmp_CONN_CVstats_grank["grank"]);
  this->CONN_CVstats_grank.count = Rcpp::as<arma::uvec>(tmp_CONN_CVstats_grank["count"]);
  this->CONN_CVstats_grank.pct = Rcpp::as<arma::vec>(tmp_CONN_CVstats_grank["pct"]);
  this->CONN_CVstats_grank.cumpct = Rcpp::as<arma::vec>(tmp_CONN_CVstats_grank["cumpct"]);
  this->CONN_CVstats_grank.mean = Rcpp::as<arma::vec>(tmp_CONN_CVstats_grank["mean"]);
  this->CONN_CVstats_grank.sum = this->CONN_CVstats_grank.mean;
  this->CONN_CVstats_grank.sum.each_col() %= arma::conv_to<arma::vec>::from(this->CONN_CVstats_grank.count);
  this->CONN_CVstats_grank.sd = Rcpp::as<arma::vec>(tmp_CONN_CVstats_grank["sd"]);
  this->CONN_CVstats_grank.q0 = Rcpp::as<arma::vec>(tmp_CONN_CVstats_grank["q0"]);
  this->CONN_CVstats_grank.q25 = Rcpp::as<arma::vec>(tmp_CONN_CVstats_grank["q25"]);
  this->CONN_CVstats_grank.q50 = Rcpp::as<arma::vec>(tmp_CONN_CVstats_grank["q50"]);
  this->CONN_CVstats_grank.q75 = Rcpp::as<arma::vec>(tmp_CONN_CVstats_grank["q75"]);
  this->CONN_CVstats_grank.q100 = Rcpp::as<arma::vec>(tmp_CONN_CVstats_grank["q100"]);
  this->CONN_CVstats_grank.ngroups = this->CONN_CVstats_grank.group_id.n_elem;
  
  Rcpp::DataFrame tmp_CONN_CVstats_length = Rcpp::as<Rcpp::DataFrame>(TRNList["CONN_CVstats_length"]);
  this->CONN_CVstats_length.group_id = Rcpp::as<arma::uvec>(tmp_CONN_CVstats_length["length"]);
  this->CONN_CVstats_length.count = Rcpp::as<arma::uvec>(tmp_CONN_CVstats_length["count"]);
  this->CONN_CVstats_length.pct = Rcpp::as<arma::vec>(tmp_CONN_CVstats_length["pct"]);
  this->CONN_CVstats_length.cumpct = Rcpp::as<arma::vec>(tmp_CONN_CVstats_length["cumpct"]);
  this->CONN_CVstats_length.mean = Rcpp::as<arma::vec>(tmp_CONN_CVstats_length["mean"]);
  this->CONN_CVstats_length.sum = this->CONN_CVstats_length.mean;
  this->CONN_CVstats_length.sum.each_col() %= arma::conv_to<arma::vec>::from(this->CONN_CVstats_length.count);
  this->CONN_CVstats_length.sd = Rcpp::as<arma::vec>(tmp_CONN_CVstats_length["sd"]);
  this->CONN_CVstats_length.q0 = Rcpp::as<arma::vec>(tmp_CONN_CVstats_length["q0"]);
  this->CONN_CVstats_length.q25 = Rcpp::as<arma::vec>(tmp_CONN_CVstats_length["q25"]);
  this->CONN_CVstats_length.q50 = Rcpp::as<arma::vec>(tmp_CONN_CVstats_length["q50"]);
  this->CONN_CVstats_length.q75 = Rcpp::as<arma::vec>(tmp_CONN_CVstats_length["q75"]);
  this->CONN_CVstats_length.q100 = Rcpp::as<arma::vec>(tmp_CONN_CVstats_length["q100"]);
  this->CONN_CVstats_length.ngroups = this->CONN_CVstats_length.group_id.n_elem;
  
  
  // Pruning
  this->CADJ_active = Rcpp::as<arma::uvec>(TRNList["CADJ_active"]);
  this->CONN_active = Rcpp::as<arma::uvec>(TRNList["CONN_active"]);
  
  
  // Topology Preserving Measures
  this->TopoProd = TRNList["TopoProd"];
  this->CADJ_TopoFxns = Rcpp::as<Rcpp::DataFrame>(TRNList["CADJ_TopoFxns"]);
  this->CONN_TopoFxns = Rcpp::as<Rcpp::DataFrame>(TRNList["CONN_TopoFxns"]);
  
  
  // DMPrune
  this->CADJ_BootSig = Rcpp::as<arma::vec>(TRNList["CADJ_BootSig"]);
  this->CONN_BootSig = Rcpp::as<arma::vec>(TRNList["CONN_BootSig"]);
  
  this->DMP_prior = Rcpp::as<arma::vec>(TRNList["DMP_prior"]);
  
  this->DMP_gamma_G = TRNList["DMP_gamma_G"];
  this->DMP_gamma_N = TRNList["DMP_gamma_N"];
  this->DMP_gamma_L = TRNList["DMP_gamma_L"];
  this->DMP_gamma_S = TRNList["DMP_gamma_S"];
  
  this->DMP_mu_G = Rcpp::as<arma::vec>(TRNList["DMP_mu_G"]);
  this->DMP_mu_N = Rcpp::as<arma::vec>(TRNList["DMP_mu_N"]);
  this->DMP_mu_L = Rcpp::as<arma::vec>(TRNList["DMP_mu_L"]);
  this->DMP_mu_S = Rcpp::as<arma::vec>(TRNList["DMP_mu_S"]);
  this->DMP_mu = Rcpp::as<arma::vec>(TRNList["DMP_mu"]);
  this->DMP_mu_rank = Rcpp::as<arma::uvec>(TRNList["DMP_mu_rank"]);
  
  this->DMP_pruneStep = Rcpp::as<arma::uvec>(TRNList["DMP_pruneStep"]);
  this->DMP_Lambda = Rcpp::as<arma::vec>(TRNList["DMP_Lambda"]);
  
  this->DMP_Lambda_error_nperms = TRNList["DMP_Lambda_error_nperms"];
  this->DMP_Lambda_error_qlo = TRNList["DMP_Lambda_error_qlo"];
  this->DMP_Lambda_error_qhi = TRNList["DMP_Lambda_error_qhi"];
  this->DMP_Lambda_error = Rcpp::as<arma::mat>(TRNList["DMP_Lambda_error"]);
  
  this->DMP_loo_order = Rcpp::as<arma::uvec>(TRNList["DMP_loo_order"]);
  this->DMP_loo_llk = Rcpp::as<arma::vec>(TRNList["DMP_loo_llk"]);
  
  
  // Flags
  this->isinit = TRNList["isinit"];
  this->isset_TPM = TRNList["isset_TPM"];
  this->isset_DMPrune = TRNList["isset_DMPrune"];
  
  return; 
}



#endif

