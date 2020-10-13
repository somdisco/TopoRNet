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


#ifndef TRN_UTILS_GROUPSTATS_HPP
#define TRN_UTILS_GROUPSTATS_HPP

namespace TRN_UTILS_GROUPSTATS {

struct groupstat_container {
  // Output containers 
  // Each has length or nrows = ngroups. 
  // All except "group" and "count" have ncols = ncol
  arma::uvec group_id; // unique group IDs found in vals_group, sorted 
  arma::uvec count;    // number of vals in each group 
  arma::vec pct;       // % of number of vals in each group 
  arma::vec cumpct;    // cum % of number of vals in each group 
  
  // Following are statistics in each group, by dimension. dim(stat matrix) = ngroups x dim; 
  arma::mat sum; 
  arma::mat mean; 
  arma::mat sd; 
  arma::mat q0; 
  arma::mat q25; 
  arma::mat q50; 
  arma::mat q75; 
  arma::mat q100; 
  
  unsigned int ngroups; 
  
  groupstat_container() {} 
};


struct groupstat_worker : public RcppParallel::Worker {
  
  // Inputs
  const arma::uvec& group_id; // group id of each row of vals matrix 
  const arma::mat& values; // values whose summaries will be computed, nX rows x p columns
  
  // Internal variables
  //arma::uvec unq_group_id; 
  //unsigned int ngroups; // number of unique group IDs found in vals_group
  unsigned int dim; // ncols(vals matrix)
  std::string stat; 
  std::vector<std::string> stat_options; // list of possible stat names that can be computed with this worker 
  
  // Output containers 
  groupstat_container groupstats; 
  
  // // Each has length or nrows = ngroups. 
  // // All except "group" and "count" have ncols = ncol
  // arma::uvec group; // unique group IDs found in vals_group, sorted 
  // arma::uvec count; // number of vals in each group 
  // arma::vec pct; // % of number of vals in each group 
  // arma::vec cumpct; // cum % of number of vals in each group 
  // // Following are statistics in each group, by dimension. dim(stat matrix) = ngroups x dim; 
  // arma::mat sum; 
  // arma::mat mean; 
  // arma::mat sd; 
  // arma::mat q0; 
  // arma::mat q25; 
  // arma::mat q50; 
  // arma::mat q75; 
  // arma::mat q100; 
  
  // Constructor
  groupstat_worker(const arma::uvec& group_id, const arma::mat& values)
    : group_id(group_id), values(values)
  {
    
    // Check the group_id vector has same number of elements as values
    if(group_id.n_elem != values.n_rows) Rcpp::stop("Input values should have nrows = length(group_id)");
    
    // Initialize stat to empty 
    stat = ""; 
    
    // Compute the unique group IDs
    groupstats.group_id = arma::sort(arma::unique(group_id));
    groupstats.ngroups = groupstats.group_id.n_elem; 
    
    dim = values.n_cols; 
    
    // Populate a vector of the possible statistics that can be computed with this worker 
    stat_options.resize(10);
    stat_options[0] = "count"; // computing count also computed pct and cumpct 
    stat_options[1] = "sum";
    stat_options[2] = "mean";
    stat_options[3] = "sd";
    stat_options[4] = "q0";
    stat_options[5] = "q25";
    stat_options[6] = "q50";
    stat_options[7] = "q75";
    stat_options[8] = "q100";
    stat_options[9] = "all";
  }
  
  void set_stat(std::string stat_) {
    
    // Check that input stat is one of the possible options 
    std::vector<std::string>::const_iterator it = std::find(stat_options.begin(), stat_options.end(), stat_);
    if(it == stat_options.end()) {
      Rcpp::stop("Input 'stat' must be one of count, sum, mean, sd, q0, q25, q50, q75, q100, all");
    } else {
      this->stat = stat_;
    }
    
    // Initialize containers, as needed 
    if(this->stat == "all" || this->stat == "count") {
      this->groupstats.count.set_size(this->groupstats.ngroups); 
      this->groupstats.count.fill(NA_REAL); 
      
      this->groupstats.pct.set_size(this->groupstats.ngroups);
      this->groupstats.pct.fill(NA_REAL);
      
      this->groupstats.cumpct.set_size(this->groupstats.ngroups); 
      this->groupstats.cumpct.fill(NA_REAL); 
    }
    
    if(this->stat == "all" || this->stat == "sum") {
      this->groupstats.sum.set_size(this->groupstats.ngroups, dim); 
      this->groupstats.sum.fill(NA_REAL); 
    }
    
    if(this->stat == "all" || this->stat == "mean") {
      this->groupstats.mean.set_size(this->groupstats.ngroups, dim); 
      this->groupstats.mean.fill(NA_REAL);
    }
    
    if(this->stat == "all" || this->stat == "sd") {
      this->groupstats.sd.set_size(this->groupstats.ngroups, dim); 
      this->groupstats.sd.fill(NA_REAL);
    }
    
    if(this->stat == "all" || this->stat == "q0") {
      this->groupstats.q0.set_size(this->groupstats.ngroups, dim); 
      this->groupstats.q0.fill(NA_REAL);
    }
    
    if(this->stat == "all" || this->stat == "q25") {
      this->groupstats.q25.set_size(this->groupstats.ngroups, dim); 
      this->groupstats.q25.fill(NA_REAL);
    }
    
    if(this->stat == "all" || this->stat == "q50") {
      this->groupstats.q50.set_size(this->groupstats.ngroups, dim); 
      this->groupstats.q50.fill(NA_REAL);
    }
    
    if(this->stat == "all" || this->stat == "q75") {
      this->groupstats.q75.set_size(this->groupstats.ngroups, dim); 
      this->groupstats.q75.fill(NA_REAL);
    }
    
    if(this->stat == "all" || this->stat == "q100") {
      this->groupstats.q100.set_size(this->groupstats.ngroups, dim); 
      this->groupstats.q100.fill(NA_REAL);
    }
    
  }
  
  // Find stats of a single group
  void summary_group_j(unsigned int j) {
    
    // Find elements in group j 
    arma::uvec membersj = arma::find(group_id == groupstats.group_id[j]); 
    
    // If this group is empty (has no members) leave the stat as NA (all stats initialized to NA in constructor)
    if(membersj.n_elem==0) return;
    
    // Otherwise, decode which stat is to be calculated, and compute it
    if(stat == "count") {
      // Note! Must set pct and cumpct after parallel processing!! 
      groupstats.count[j] = membersj.n_elem; 
      return; 
    } 
    
    if(stat == "sum") {
      groupstats.sum.row(j) = arma::sum(values.rows(membersj), 0);
      return; 
    } 
    
    if(stat == "mean") {
      groupstats.mean.row(j) = arma::mean(values.rows(membersj), 0);
      return; 
    } 
    
    if(stat == "sd") {
      if(membersj.n_elem == 1) {
        groupstats.sd.row(j).fill(0);
      } else {
        groupstats.sd.row(j) = arma::stddev(values.rows(membersj), 0, 0);
      }
      return; 
    } 
    
    if(stat == "q0") {
      groupstats.q0.row(j) = arma::min(values.rows(membersj), 0);
      return; 
    } 
    
    if(stat == "q25") {
      arma::vec probs = {0.25};
      arma::mat quants = arma::quantile(values.rows(membersj), probs, 0);
      groupstats.q25.row(j) = quants.row(0);
      return; 
    } 
    
    if(stat == "q50") {
      groupstats.q50.row(j) = arma::median(values.rows(membersj), 0);
      return; 
    } 
    
    if(stat == "q75") {
      arma::vec probs = {0.75};
      arma::mat quants = arma::quantile(values.rows(membersj), probs, 0);
      groupstats.q75.row(j) = quants.row(0);
      return; 
    } 
    
    if(stat == "q100") {
      groupstats.q100.row(j) = arma::max(values.rows(membersj), 0);
      return; 
    }
    
    // Otherwise, stat must = "all"
    groupstats.count[j] = membersj.n_elem; 
    groupstats.sum.row(j) = arma::sum(values.rows(membersj), 0);
    groupstats.mean.row(j) = groupstats.sum.row(j) / double(groupstats.count[j]); 
    if(membersj.n_elem == 1) {
      groupstats.sd.row(j).fill(0);
    } else {
      groupstats.sd.row(j) = arma::stddev(values.rows(membersj), 0, 0);
    }
    
    arma::vec probs = {0.00, 0.25, 0.50, 0.75, 1.00};
    arma::mat quants = arma::quantile(values.rows(membersj), probs, 0);
    groupstats.q0.row(j) = quants.row(0);
    groupstats.q25.row(j) = quants.row(1);
    groupstats.q50.row(j) = quants.row(2);
    groupstats.q75.row(j) = quants.row(3);
    groupstats.q100.row(j) = quants.row(4);
    
    return;
  }
  
  
  // Parallel operator - find BMU of each row of X in parallel
  void operator()(std::size_t begin, std::size_t end) {
    for(unsigned int j = begin; j < end; j++) {
      summary_group_j(j);
    }
  }
  
  // Parallel call method
  void calc_parallel() {
    if(stat == "") Rcpp::stop("Must call set_stat before computing stats");
    
    RcppParallel::parallelFor(0, groupstats.ngroups, *this);
    
    // Compute the pct and cumpct, if stat = count or all 
    if(stat == "all" || stat == "count") {
      groupstats.pct = arma::conv_to<arma::vec>::from(groupstats.count) / double(values.n_rows);
      groupstats.cumpct = arma::cumsum(groupstats.pct); 
    }
  }
  
  // Non-parallel call method
  void calc_serial() {
    if(stat == "") Rcpp::stop("Must call set_stat before computing stats");
    // Find BMU of each row of X
    for(unsigned int j=0; j<groupstats.ngroups; ++j) {
      summary_group_j(j);
    }
    
    // Compute the pct and cumpct, if stat = count or all 
    if(stat == "all" || stat == "count") {
      groupstats.pct = arma::conv_to<arma::vec>::from(groupstats.count) / double(values.n_rows);
      groupstats.cumpct = arma::cumsum(groupstats.pct); 
    }
    
  }
};





} // close namespace 

#endif


