#ifndef TOPORNET_TRNOBJ_HPP
#include "TopoRNet_types.hpp"
#endif



//RCPP_EXPOSED_CLASS(TRN)

RCPP_MODULE(trn_module){
  using namespace Rcpp; // Added (if not done globally)
  
  class_<TRN>("TRN")
    
    .constructor()
    
    // Initializer
    .method("initialize_TRN", &TRN::initialize_TRN, "Initialize a TRN object")

    // Fields involving vertex and edge lists 
    .field_readonly("nV", &TRN::nV, "Number of vertices in the graph")
    .field_readonly("CADJ_nE", &TRN::CADJ_nE, "Number of edges in CADJ graph")
    .field_readonly("CONN_nE", &TRN::CONN_nE, "Number of edges in CONN graph")
    .field_readonly("OT_nE", &TRN::OT_nE, "Number of edges in Output Topology")
    
    .field_readonly("CADJ_EL", &TRN::CADJ_EL, "CADJ edge list")
    .field_readonly("CONN_EL", &TRN::CONN_EL, "CONN edge list")
    .field_readonly("OT_EL", &TRN::OT_EL, "Output Topology edge list")
    
    .field_readonly("CADJ", &TRN::CADJ, "CADJ edge weights")
    .field_readonly("CONN", &TRN::CONN, "CONN edge weights")
    
    .field_readonly("CADJ_deg", &TRN::CADJ_deg, "Vertex degree of CADJ graph")
    .field_readonly("CONN_deg", &TRN::CONN_deg, "Vertex degree of CONN graph")
    .field_readonly("OT_deg", &TRN::OT_deg, "Vertex degree of Output Topology")
    
    .field_readonly("CADJ_wdeg", &TRN::CADJ_wdeg, "Weighted vertex degree of CADJ graph")
    .field_readonly("CONN_wdeg", &TRN::CONN_wdeg, "Weighted vertex degree of CONN graph")
    
    .field_readonly("CADJ_fwdflen", &TRN::CADJ_fwdflen, "Forward folding lengths of CADJ edges")
    .field_readonly("CADJ_bwdflen", &TRN::CADJ_bwdflen, "Backward folding lengths of Output Topology on CADJ")
    .field_readonly("CONN_fwdflen", &TRN::CONN_fwdflen, "Forward folding lengths of CONN edges")
    .field_readonly("CONN_bwdflen", &TRN::CONN_bwdflen, "Backward folding lengths of Output Topology on CONN")
    
    .field_readonly("CADJ_TP_radius", &TRN::CADJ_TP_radius, "Output Topology radius that would completely preserve CADJ topology")
    .field_readonly("CONN_TP_radius", &TRN::CONN_TP_radius, "Output Topology radius that would completely preserve CONN topology")
    
    .field_readonly("OT_geodesic_dist", &TRN::OT_geodesic_dist, "Geodesic distance between output vertices")
    .field_readonly("OT_max_dist", &TRN::OT_max_dist, "The maximum output topology distance")   
    .field_readonly("OT_nhb_sizes", &TRN::OT_nhb_sizes, "The cumulative number of neighbors of each TRN vertex at all possible geodesic distances in its output topology")
    
    //.method("update_CONN", &TRN::update_CONN, "Update CONN based on active CADJ edges")
    
    
    // CONNvis fields and methods 
    .field_readonly("CADJ_lrank", &TRN::CADJ_lrank, "CADJ local edge ranks")
    .field_readonly("CADJ_grank", &TRN::CADJ_grank, "CADJ global edge ranks")
    .method("get_CADJvis_stats", &TRN::get_CADJvis_stats, "Get the CADJvis edge statistics")
    .method("get_CADJvis_colors", &TRN::get_CADJvis_colors, "Get the CADJvis edge colors")
    .method("get_CADJvis_widths", &TRN::get_CADJvis_widths, "Get the CADJvis edge widths")
    
    .field_readonly("CONN_lrank", &TRN::CONN_lrank, "CONN local edge ranks")
    .field_readonly("CONN_grank", &TRN::CONN_grank, "CONN global edge ranks")
    .method("get_CONNvis_stats", &TRN::get_CONNvis_stats, "Get the CONNvis edge statistics")
    .method("get_CONNvis_colors", &TRN::get_CONNvis_colors, "Get the CONNvis edge colors")
    .method("get_CONNvis_widths", &TRN::get_CONNvis_widths, "Get the CONNvis edge widths")
    
    
  
    // Pruning 
    .field_readonly("CADJ_active", &TRN::CADJ_active, "Flag marking each CADJ edge as active (un-pruned, not-thresholded)")
    .method("prune_CADJ_edge_list", &TRN::prune_CADJ_edge_list, "Prune CADJ edges, given an edge list")
    .method("prune_CADJ_edge_id", &TRN::prune_CADJ_edge_id, "Prune CADJ edges, given a list of edge ids")
    .method("prune_CADJ_vertex_id", &TRN::prune_CADJ_vertex_id, "Prune CADJ edges, given a list of vertex ids")
    .method("prune_CADJ_edge_weight", &TRN::prune_CADJ_edge_weight, "Prune CADJ edges, given a minimum edge weight")
    .method("prune_CADJ_lrank", &TRN::prune_CADJ_lrank, "Prune CADJ edges, given a maximum edge local rank")
    .method("prune_CADJ_grank", &TRN::prune_CADJ_grank, "Prune CADJ edges, given a maximum edge global rank")
    .method("prune_CADJ_length", &TRN::prune_CADJ_length, "Prune CADJ edges, given a maximum edge fwdflen")
    .method("prune_CADJ_wdeg", &TRN::prune_CADJ_wdeg, "Prune CADJ edges involving vertices with a minimum wdeg")
    .method("get_CADJ", &TRN::get_CADJ, "Get CADJ adjacency matrix with only active edges")
    .method("restore_CADJ_edges", &TRN::restore_CADJ_edges, "Reset all CADJ edges to active (un-pruned)")

    .field_readonly("CONN_active", &TRN::CONN_active, "Flag marking each CONN edge as active (un-pruned, not-thresholded)")
    .method("prune_CONN_edge_list", &TRN::prune_CONN_edge_list, "Prune CONN edges, given an edge list")
    .method("prune_CONN_edge_id", &TRN::prune_CONN_edge_id, "Prune CONN edges, given a list of edge ids")
    .method("prune_CONN_vertex_id", &TRN::prune_CONN_vertex_id, "Prune CONN edges, given a list of vertex ids")
    .method("prune_CONN_edge_weight", &TRN::prune_CONN_edge_weight, "Prune CONN edges, given a minimum edge weight")
    .method("prune_CONN_lrank", &TRN::prune_CONN_lrank, "Prune CONN edges, given a maximum edge local rank")
    .method("prune_CONN_grank", &TRN::prune_CONN_grank, "Prune CONN edges, given a maximum edge global rank")
    .method("prune_CONN_length", &TRN::prune_CONN_length, "Prune CONN edges, given a maximum edge fwdflen")
    .method("prune_CONN_wdeg", &TRN::prune_CONN_wdeg, "Prune CONN edges involving vertices with a minimum wdeg")
    .method("prune_CONN_CADJ", &TRN::prune_CONN_CADJ, "Prune CONN edges that are inactive CADJ.")
    .method("get_CONN", &TRN::get_CONN, "Get CONN adjacency matrix with only active edges")
    .method("restore_CONN_edges", &TRN::restore_CONN_edges, "Reset all CONN edges to active (un-pruned)")

    .method("get_OTADJ", &TRN::get_OTADJ, "Get adjacency matrix of output topology")
     
    .method("CADJ_active_verts", &TRN::CADJ_active_verts, "List of active vertices in CADJ")
    .method("CONN_active_verts", &TRN::CONN_active_verts, "List of active vertices in CONN")
     
     
    // TPM
    .field_readonly("TopoProd", &TRN::TopoProd, "Topographic Product")
    .field_readonly("CADJ_TopoFxns", &TRN::CADJ_TopoFxns, "Topographic Functions computed on CADJ")
    .field_readonly("CONN_TopoFxns", &TRN::CONN_TopoFxns, "Topographic Functions computed on CONN")
    .method("calc_TopoMeasures", &TRN::calc_TopoMeasures, "Set the topology preservation measures")
    
    // *** DM-Prune
    //.field_readonly("DMP_CADJ_prior", &TRN::DMP_CADJ_prior, "CADJ prior for DM-Prune")
    .field_readonly("CADJ_BootSig", &TRN::CADJ_BootSig, "Bootstrap significance of CADJ edges")
    .field_readonly("CONN_BootSig", &TRN::CONN_BootSig, "Bootstrap significance of CONN edges")
    .method("set_BootSig", &TRN::set_BootSig, "Set the Bootstrap significance of CADJ & CONN edges")
  
    .field("DMP_prior", &TRN::DMP_prior, "CADJ prior for DM-Prune")
    .field_readonly("DMP_mu_G", &TRN::DMP_mu_G, "CADJ mu_G score for DM-Prune")
    .field_readonly("DMP_mu_N", &TRN::DMP_mu_N, "CADJ mu_N score for DM-Prune")
    .field_readonly("DMP_mu_L", &TRN::DMP_mu_L, "CADJ mu_L score for DM-Prune")
    .field_readonly("DMP_mu_S", &TRN::DMP_mu_S, "CADJ mu_S score for DM-Prune")
    .field_readonly("DMP_mu", &TRN::DMP_mu, "CADJ mu score for DM-Prune")
    .field_readonly("DMP_mu_rank", &TRN::DMP_mu_rank, "Rank of CADJ mu score for DM-Prune")
    
    .field_readonly("DMP_pruneStep", &TRN::DMP_pruneStep, "DM-Prune Steps")
    .field_readonly("DMP_Lambda", &TRN::DMP_Lambda, "DM-Prune Lambda")
    .method("calc_DMPrune_LambdaPath", &TRN::calc_DMPrune_LambdaPath, "Set the DM-Prune Lambda Path")
    
    .field_readonly("DMP_Lambda_error", &TRN::DMP_Lambda_error, "Error bands for Lambda path")
    .method("calc_DMPrune_LambdaPath_error", &TRN::calc_DMPrune_LambdaPath_error, "Calculate the Lambda path error bands")
  
    .field_readonly("DMP_LOO_pruneStep", &TRN::DMP_loo_order, "loo ordering")
    .field_readonly("DMP_LOO_Lambda", &TRN::DMP_loo_llk, "loo llk")
    .method("calc_DMPrune_LambdaPath_LOO", &TRN::calc_DMPrune_LambdaPath_loo, "Set loo lambda path")
  
    .method("DMPrune_CADJ_step", &TRN::DMPrune_CADJ_step, "DM-Prune CADJ edges based on minimum step")
    
    // Flags
    .field_readonly("isinit", &TRN::isinit, "Flag indicating whether $initialize_TRN has been called")
    .field_readonly("isset_TPM", &TRN::isset_TPM, "Flag indicating whether $calc_TopoMeasures has been called")
    .field_readonly("isset_DMPrune", &TRN::isset_DMPrune, "Flag indicating whether $calc_DMPrune_LambdaPath has been called")

    // IO 
    .method("save", &TRN::save, "Save a TRN object to disk")
    .method("load", &TRN::load, "Load a TRN object from disk")
    .method("as_list", &TRN::as_list, "Convert a TRN object to an R list")
    .method("load_list", &TRN::load_list, "Populate a TRN object from an R list")
    
    ;
}
