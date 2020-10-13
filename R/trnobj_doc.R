#' @name TRN
#' @title The TRN object class
#' 
#' @field parallel Flag indicating whether computations should be performed in parallel, using the RcppParallel package.
#' @field nV Number of vertices in the TRN.
#' @field CADJ_nE Number of edges in the CADJ graph.
#' @field CONN_nE Number of edges in the CONN graph.
#' @field OT_nE Number of edges in the output space topology of the TRN (i.e., the SOM lattice).
#' @field CADJ_EL List of CADJ edges. 
#' This is a 2-column matrix whose rows define the vertices connected by each CADJ edge. 
#' @field CONN_EL List of CONN edges. 
#' This is a 2-column matrix whose rows define the vertices connected by each CONN edge. 
#' @field OT_EL List of edges in the output space topology of the TRN. 
#' This is a 2-column matrix whose rows define the vertices connected by each output space edge.
#' @field CADJ Weights of each CADJ edge defined in \code{CADJ_EL}.
#' This is a vector of \code{length = nrow(CADJ_EL)}. 
#' @field CONN Weights of each CONN edge defined in \code{CONN_EL}.
#' This is a vector of \code{length = nrow(CONN_EL)}. 
#' @field CADJ_deg The degree (number of incident edges) of each TRN vertex, according to the CADJ graph. 
#' @field CONN_deg The degree (number of incident edges) of each TRN vertex, according to the CONN graph. 
#' @field OT_deg The degree (number of incident edges) of each TRN vertex, according to the output space topology. 
#' @field CADJ_wdeg The weighted degree of each TRN vertex (sum of edge weights emanating from each vertex), according to the CADJ graph. 
#' @field CONN_wdeg The weighted degree of each TRN vertex (sum of edge weights emanating from each vertex), according to the CONN graph. 
#' @field CADJ_fwdflen The forward folding length of each CADJ edge listed in \code{CADJ_EL}. 
#' This is the geodesic distance of each edge, measured on the TRN's output topology. 
#' @field CONN_fwdflen The forward folding length of each CONN edge listed in \code{CONN_EL}. 
#' This is the geodesic distance of each edge, measured on the TRN's output topology. 
#' @field CADJ_bwdflen The backward folding length (geodesic distance) of each edge in the TRN's output topology, measured on the CADJ graph. 
#' \code{length = nrow(OT_EL)}. 
#' @field CONN_bwdflen The backward folding length (geodesic distance) of each edge in the TRN's output topology, measured on the CONN graph. 
#' \code{length = nrow(OT_EL)}. 
#' @field CADJ_TP_radius The neighborhood size on the TRN's output topology which would preserve all of the topological adjacencies 
#' found in the CADJ graph, if all of these adjacencies were packed into a neighborhood of this radius.  
#' @field CONN_TP_radius The neighborhood size on the TRN's output topology which would preserve all of the topological adjacencies 
#' found in the CONN graph, if all of these adjacencies were packed into a neighborhood of this radius.  
#' @field OT_geodesic_dist The geodesic distance between TRN vertices, measured according to its output topology. 
#' @field OT_max_dist The maximum geodesic distance between any two TRN vertices, measured according to its output topology. 
#' @field OT_nhb_sizes The cumulative number of neighbors of each TRN vertex at all possible geodesic distances in its output topology,  
#' stored in a \code{nV} x \code{OT_max_dist} matrix. 
#' @field CADJ_lrank The CADJvis local rank of each CADJ edge in \code{CADJ_EL}. 
#' @field CADJ_grank The CADJvis global rank of each CADJ edge in \code{CADJ_EL}. 
#' @field CONN_lrank The CONNvis local rank of each CONN edge in \code{CONN_EL}. 
#' @field CONN_grank The CONNvis global rank of each CONN edge in \code{CONN_EL}. 
#' @field CADJ_active Indicator flag identifying whether each edge in \code{CADJ_EL} is active (=1) or not (=0). 
#' "active" = un-pruned, i.e. not removed by thresholding (so all edges are active until any of the \code{prune_CADJ_*} methods are called). 
#' @field CONN_active Indicator flag identifying whether each edge in \code{CONN_EL} is active (=1) or not (=0). 
#' "active" = un-pruned, i.e. not removed by thresholding (so all edges are active until any of the \code{prune_CONN_*} methods are called). 
#' @field TopoProd The Topographic Product of the TRN. See \code{?set_TopoPresMeasures} for details. 
#' @field CADJ_TopoFxns The suite of Topographic Functions of the CADJ graph. See \code{?set_TopoPresMeasures} for details. 
#' @field CONN_TopoFxns The suite of Topographic Functions of the CONN graph. See \code{?set_TopoPresMeasures} for details. 
#' @field CADJ_BootSig The significance values for each edge in \code{CADJ_EL}, as derived from a bootstrapped re-sampling procedure. 
#' @field CONN_BootSig The significance values for each edge in \code{CONN_EL}, as derived from a bootstrapped re-sampling procedure. 
#' @field DMP_prior The DM-Prune priors for each edge in \code{CADJ_EL}. 
#' @field DMP_mu_G The DM-Prune global \eqn{\mu} score for each edge in \code{CADJ_EL}. 
#' @field DMP_mu_N The DM-Prune local \eqn{\mu} score for each edge in \code{CADJ_EL}. 
#' @field DMP_mu_L The DM-Prune length \eqn{\mu} score for each edge in \code{CADJ_EL}. 
#' @field DMP_mu_S The DM-Prune bootstrap significance \eqn{\mu} score for each edge in \code{CADJ_EL}. 
#' @field DMP_mu The DM-Prune composite \eqn{\mu} score for each edge in \code{CADJ_EL}. 
#' @field DMP_mu_rank The rank of each edge in \code{CADJ_EL}, according to their DM-Prune composite \eqn{\mu} scores. 
#' The ranking mechanism is ascending and dense, meaning the edge with the lowest \eqn{\mu} has \code{rank = 1}, 
#' and any edges with identical \eqn{mu} scores share the same rank.   
#' @field DMP_pruneStep The unique prune steps at which the DM-Prune \eqn{\Lambda}-path is calculated. 
#' These steps are the unique values in \code{DMP_mu_rank}.  
#' @field DMP_Lambda The DM-Prune \eqn{\Lambda} path, computed at each step in \code{DMP_pruneStep}. 
#' 
#' @field isinit Flag indicating whether \code{initialize_TRN} has been called.
#' @field isset_TPM Flag indicating whether \code{calc_TopoMeasures} has been called. 
#' @field isset_DMPrune Flag indicating whether \code{calc_DMPrune_LambdaPath} has been called.
#'
#' @section Methods:
#' Each class method has its own documentation, accessible via \code{?TopoRNet::<method_name>}. 
#' For completeness, the list is repeated here in entirety.  Additional functionality 
#' for visualizing a TRN object is available through the \code{vis_*} functions. 
#' See their documentation for more information. 
#' 
#' \describe{
#' \item{\code{TRN$new}}{Instantiate a TRN object.}
#' \item{\code{set_parallel}}{Set the parallel computation flag.}
#' \item{\code{initialize_TRN}}{Initialize a TRN object with CADJ and output topology adjacency matrices.}
# \item{\code{update_CONN}}{Update CONN based on active CADJ edges.}
#' \item{\code{get_CADJvis_stats}}{Get the CADJvis summary statistics, grouped by \code{lrank}, \code{grank}, or \code{length}.}
#' \item{\code{get_CADJvis_colors}}{Get the CADJvis colors of each edge in \code{CADJ_EL}, used for producing the CADJvis.}
#' \item{\code{get_CADJvis_widths}}{Get the CADJvis line widths of each edge in \code{CADJ_EL}, used for producing the CADJvis.}
#' \item{\code{get_CONNvis_stats}}{Get the CONNvis summary statistics, grouped by \code{lrank}, \code{grank}, or \code{length}.}
#' \item{\code{get_CONNvis_colors}}{Get the CONNvis colors of each edge in \code{CONN_EL}, used for producing the CONNvis.}
#' \item{\code{get_CONNvis_widths}}{Get the CONNvis line widths of each edge in \code{CONN_EL}, used for producing the CONNvis.}
#' \item{\code{prune_CADJ_edge_list}}{Prune the CADJ graph given a list of edges to remove.}
#' \item{\code{prune_CADJ_edge_id}}{Prune the CADJ graph given a list of edge IDs to remove.}
#' \item{\code{prune_CADJ_vertex_id}}{Prune the CADJ graph by removing edges connecting vertices in a list of given vertex IDs.}
#' \item{\code{prune_CADJ_edge_weight}}{Prune the CADJ graph by removing edges whose weights are < a given minimum edge weight.}
#' \item{\code{prune_CADJ_lrank}}{Prune the CADJ graph by removing edges whose \code{CADJ_lrank} is > a given maximum lrank.}
#' \item{\code{prune_CADJ_grank}}{Prune the CADJ graph by removing edges whose \code{CADJ_grank} is > a given maximum grank.}
#' \item{\code{prune_CADJ_length}}{Prune the CADJ graph by removing edges whose \code{CADJ_fwdflen} is > a given maximim forward folding length.}
#' \item{\code{prune_CADJ_wdeg}}{Prune the CADJ graph by removing edges incident to vertices whose \code{CADJ_wdeg} is < a given minimum weighted degree.}
#' \item{\code{get_CADJ}}{Get the CADJ adjacency matrix of all \code{CADJ_active} edges (i.e., the adjacency matrix returned excludes any pruned edges).}
#' \item{\code{restore_CADJ_edges}}{Set all edges in \code{CADJ_EL} as active, which erases any graph pruning that has been done via the \code{prune_CADJ_*} methods.}
#' \item{\code{prune_CONN_edge_list}}{Prune the CONN graph given a list of edges to remove.}
#' \item{\code{prune_CONN_edge_id}}{Prune the CONN graph given a list of edge IDs to remove.}
#' \item{\code{prune_CONN_vertex_id}}{Prune the CONN graph by removing edges connecting vertices in a list of given vertex IDs.}
#' \item{\code{prune_CONN_edge_weight}}{Prune the CONN graph by removing edges whose weights are < a given minimum edge weight.}
#' \item{\code{prune_CONN_lrank}}{Prune the CONN graph by removing edges whose \code{CONN_lrank} is > a given maximum lrank.}
#' \item{\code{prune_CONN_grank}}{Prune the CONN graph by removing edges whose \code{CONN_grank} is > a given maximum grank.}
#' \item{\code{prune_CONN_length}}{Prune the CONN graph by removing edges whose \code{CONN_fwdflen} is > a given maximim forward folding length.}
#' \item{\code{prune_CONN_wdeg}}{Prune the CONN graph by removing edges incident to vertices whose \code{CONN_wdeg} is < a given minimum weighted degree.}
#' \item{\code{prune_CONN_CADJ}}{Prune the CONN graph by removing CADJ inactive edges.}
#' \item{\code{get_CONN}}{Get the CONN adjacency matrix of all \code{CONN_active} edges (i.e., the adjacency matrix returned excludes any pruned edges).}
#' \item{\code{restore_CONN_edges}}{Set all edges in \code{CONN_EL} as active, which erases any graph pruning that has been done via the \code{prune_CONN_*} methods.}
#' \item{\code{get_OTADJ}}{Get the adjacency matrix of the output space topology of the TRN. This is the same adjacency that was given (and stored) during a call to \code{initialize_TRN}.}
#' \item{\code{CADJ_active_verts}}{Return a list of vertex IDs that are connected by any \code{CADJ_active} edge.}
#' \item{\code{CONN_active_verts}}{Return a list of vertex IDs that are connected by any \code{CONN_active} edge.}
#' \item{\code{calc_TopoMeasures}}{Sets the Topology Preserving Measures of both the CADJ and CONN graphs (Topographic Product and several Topographic Functions).}
#' \item{\code{set_BootSig}}{Set the Bootstrap significance values of each CADJ and CONN edge.}
#' \item{\code{calc_DMPrune_LambdaPath}}{Calculate the DM-Prune Lambda path.}
#' \item{\code{DMPrune_CADJ_step}}{Prune the CADJ graph by removing edges in \code{CADJ_EL} whose \code{DMP_mu_rank} is < a given minimum prune step.}
#' \item{\code{save}}{Save a TRN object to disk}
#' \item{\code{load}}{Load a previously saved TRN object from disk}
#' \item{\code{as_list}}{Convert and return all fields of a TRN object to an R list.}
#' \item{\code{load_list}}{Populate an instance of a TRN object from an R list.}
#' }
#' 
NULL


# ***** Constructor *****

#' @name new
#' @title Create an empty TRN object
#' @return an empty TRN object templated for initialization (via \link{initialize_TRN}).
#' @usage TRN$new()
NULL 

#' @name set_parallel
#' @title Setter function for parallel computation.  
#' @description Sets an internal flag controlling whether computations are performed in parallel. 
#' @param parallel either TRUE or FALSE, as desired
#' @details 
#' Parallel computation is supported via RcppParallel, 
#' see \link[RcppParallel]{setThreadOptions} for details to control threading. 
#' By default, \code{parallel = TRUE} at TRNobj instantiation. 
#' @return None, \code{parallel} field is set internally 
#' @usage TRNobj$set_parallel(parallel)
NULL

#' @name initialize_TRN
#' @title Initialize a TRN object 
#' @description Sets the CADJ, CONN, and output topology graphs. 
#' @param CADJ The CADJ adjacency matrix of the TRN. 
#' @param OTADJ The adjacency matrix of the TRN's output topology (e.g., the SOM neuron lattice adjacency). 
#' Only binary output topology adjacency matrices are supported at this time. 
#' (\eqn{OTADJ_{ij} = 1} indicate an edge between vertices \code{i} and \code{j}). 
#' @details 
#' The following fields are computed and stored internally: 
#' \itemize{
#' \item \code{nV}, set = nrow(CADJ). 
#' \item \code{CADJ_EL} and \code{CONN_EL}, the CADJ and CONN edge lists. 
#' \item \code{CADJ} and \code{CONN}, the CADJ and CONN edge weights. 
#' \item \code{CADJ_nE} and \code{CONN_nE}, the number of edges in the CADJ and CONN graphs. 
#' \item \code{CADJ_deg} and \code{CONN_deg}, the degree of each TRN vertex in the CADJ and CONN graphs. 
#' \item \code{CADJ_wdeg} and \code{CONN_wdeg}, the weighted degree of each TRN vertex in the CADJ and CONN graphs. 
#' \item \code{OT_EL} and \code{OT_nE} and \code{OT_deg} the edge list, number of edges, and vertex degrees in the TRN's output topology. 
#' \item \code{OT_max_dist}, the maximum geodesic distance of any two vertices in the TRN's output topology. 
#' \item \code{CADJ_fwdflen} and \code{CONN_fwdflen}, the forward folding lengths of each CADJ and CONN edge. 
#' \item \code{CADJ_bwdflen} and \code{CONN_bwdflen}, the backward folding lengths of each output topology edge on the CADJ and CONN graphs. 
#' \item \code{CADJ_TP_radius} and \code{CONN_TP_radius}, the Topology Preserving Radius of the CADJ and CONN graphs. 
#' \item \code{CADJ_lrank} and \code{CONN_lrank} The local ranks of each CADJ and CONN edge. 
#' \item \code{CADJ_grank} and \code{CONN_grank} The global ranks of each CADJ and CONN edge. 
#' }
#' @details Both \code{CADJ} and \code{OTADJ} should be square, and have the same dimensions. 
#' In addition to the above list of internally computed quantities, \code{initialize_TRN} also computes the 
#' CADJvis and CONNvis summary statistics (retrievable via \link{get_CADJvis_stats} and \link{get_CONNvis_stats}, respectively), and 
#' sets all \code{CADJ_active} and \code{CONN_active} flags = 1.  
#' @return None
#' @usage TRNobj$initialize_TRN(CADJ, OTADJ)
#' \insertRef{MartinetzSchulten1994}{TopoRNet}
NULL 


# @name update_CONN
# @title Update CONN 
# @description Updates the CONN Topology and CONNvis statistics based on \code{CADJ_active} edges 
# @return None
# @usage TRNobj$update_CONN
#NULL 


# ***** CONNvis *****

#' @name get_CADJvis_stats
#' @title Get the CADJvis summary statistics. 
#' @description The CADJvis summary statistics are computed for each unique grouping found in \code{CADJ_lrank}, \code{CADJ_grank}, and \code{CADJ_fwdflen}. 
#' @param group a string identifying which summary statistics to return. Can be one of 'lrank', 'grank', or 'length'. 
#' 
#' @return A data frame with summary statistics for each group in its rows, and columns: 
#' \itemize{
#' \item \code{lrank} (or \code{grank} or \code{length} as applicable), identifying the grouping for each set of summary stats
#' \item \code{count} the number of CADJ edges in the group 
#' \item \code{pct} the proportion of the total number of CADJ edges in each group 
#' \item \code{cumpct} the cumulative proportion of CADH edges in each group (as ordered by the grouping variable)
#' \item \code{mean} the average of CADJ edge weights in each group 
#' \item \code{sd} the standard deviation of CADJ edge weights in each group 
#' \item \code{q0} the 0.00 quantile (minimum) of CADJ edge weights in each group 
#' \item \code{q25} the 0.25 quantile (Q1) of CADJ edge weights in each group 
#' \item \code{q50} the 0.50 quantile (median) of CADJ edge weights in each group 
#' \item \code{q75} the 0.75 quantile (Q3) of CADJ edge weights in each group
#' \item \code{q100} the 1.00 quantile (maximum) of CADJ edge weights in each group
#' }
#' @usage TRNobj$get_CADJvis_stats(group)
#' @references 
#' \insertRef{TasdemirMerenyi2009}{TopoRNet}
NULL 

#' @name get_CADJvis_colors
#' @title Get the CADJvis edge colors
#' @description The edge colorings (red, blue, green, yellow, grayscale) are assigned in increasing order of an edge's \code{CADJ_lrank}. 
#' @return A vector containing the CADJvis (hex) color of each edge in \code{CADJ_EL}.  
#' @usage TRNobj$get_CADJvis_colors()
#' @references 
#' \insertRef{TasdemirMerenyi2009}{TopoRNet}
NULL 

#' @name get_CADJvis_widths
#' @title Get the CADJvis edge widths
#' @description The edge widths are assigned in reverse of each edge's \code{CADJ_grank}. That is, if a CADJ graph has unique global ranks 
#' of 1,2,3,4,5, then edges with those global ranks will have widths = 5,4,3,2,1.  
#' @return A vector containing the CADJvis widths of each edge in \code{CADJ_EL}.  
#' @usage TRNobj$get_CADJvis_widths()
#' @references 
#' \insertRef{TasdemirMerenyi2009}{TopoRNet}
NULL 

#' @name get_CONNvis_stats
#' @title Get the CONNvis summary statistics. 
#' @description The CONNvis summary statistics are computed for each unique grouping found in \code{CONN_lrank}, \code{CONN_grank}, and \code{CONN_fwdflen}. 
#' @param group a string identifying which summary statistics to return. Can be one of 'lrank', 'grank', or 'length'. 
#' 
#' @return A data frame with summary statistics for each group in its rows, and columns: 
#' \itemize{
#' \item \code{lrank} (or \code{grank} or \code{length} as applicable), identifying the grouping for each set of summary stats
#' \item \code{count} the number of CONN edges in the group 
#' \item \code{pct} the proportion of the total number of CONN edges in each group 
#' \item \code{cumpct} the cumulative proportion of CADH edges in each group (as ordered by the grouping variable)
#' \item \code{mean} the average of CONN edge weights in each group 
#' \item \code{sd} the standard deviation of CONN edge weights in each group 
#' \item \code{q0} the 0.00 quantile (minimum) of CONN edge weights in each group 
#' \item \code{q25} the 0.25 quantile (Q1) of CONN edge weights in each group 
#' \item \code{q50} the 0.50 quantile (median) of CONN edge weights in each group 
#' \item \code{q75} the 0.75 quantile (Q3) of CONN edge weights in each group
#' \item \code{q100} the 1.00 quantile (maximum) of CONN edge weights in each group
#' }
#' @usage TRNobj$get_CONNvis_stats(group)
#' @references 
#' \insertRef{TasdemirMerenyi2009}{TopoRNet}
NULL 

#' @name get_CONNvis_colors
#' @title Get the CONNvis edge colors
#' @description The edge colorings (red, blue, green, yellow, grayscale) are assigned in increasing order of an edge's \code{CONN_lrank}. 
#' @return A vector containing the CONNvis (hex) color of each edge in \code{CONN_EL}.  
#' @usage TRNobj$get_CONNvis_colors()
#' @references 
#' \insertRef{TasdemirMerenyi2009}{TopoRNet}
NULL 

#' @name get_CONNvis_widths
#' @title Get the CONNvis edge widths
#' @description The edge widths are assigned in reverse of each edge's \code{CONN_grank}. That is, if a CONN graph has unique global ranks 
#' of 1,2,3,4,5, then edges with those global ranks will have widths = 5,4,3,2,1.  
#' @return A vector containing the CONNvis widths of each edge in \code{CONN_EL}.  
#' @usage TRNobj$get_CONNvis_widths()
#' @references 
#' \insertRef{TasdemirMerenyi2009}{TopoRNet}
NULL 



# ***** Pruning *****

#' @name prune_CADJ_edge_list
#' @title Prune CADJ by an edge list 
#' @description Prune the CADJ graph given a list of edges to remove. 
#' @param edge_list a 2-column matrix rows define the vertices connected by each CADJ edge to prune 
#' @details The vertex indices in \code{edge_list} should be 1-based. 
#' @return None
#' @usage TRNobj$prune_CADJ_edge_list(edge_list)
NULL 


#' @name prune_CADJ_edge_id
#' @title Prune CADJ by edge ids
#' @description Prune the CADJ graph given a list of edge IDs to remove.
#' @param edge_id a vector containing the rows indices of \code{CADJ_EL} which contain the edges desired to prune. 
#' @return None
#' @usage TRNobj$prune_CADJ_edge_id(edge_id)
NULL 

#' @name prune_CADJ_vertex_id
#' @title Prune CADJ by a set of vertex ids
#' @description Prune the CADJ graph by removing edges connecting vertices in a list of given vertex IDs.
#' @param vertex_id a vector containing the vertices to prune. 
#' All connections in \code{CADJ_EL} from or to the vertices in \code{vertex_id} will be pruned. 
#' @return None
#' @usage TRNobj$prune_CADJ_vertex_id(vertex_id)
NULL 

#' @name prune_CADJ_edge_weight
#' @title Prune CADJ by a minimum edge weight
#' @description Prune the CADJ graph by removing edges whose weights are \strong{strictly less than} a given minimum edge weight.
#' @param min_weight minimum weight allowed in the pruned graph. 
#' Any edges in \code{CADJ_EL} whose \code{CADJ} value is < \code{min_weight} will be pruned. 
#' @return None
#' @usage TRNobj$prune_CADJ_edge_weight(min_weight)
NULL 

#' @name prune_CADJ_lrank
#' @title Prune CADJ by a maximum local rank
#' @description Prune the CADJ graph by removing edges whose \code{CADJ_lrank} is \strong{strictly greater than} a given maximum lrank.
#' @param max_lrank maximum \code{CADJ_lrank} allowed in the pruned graph. 
#' Any edges whose \code{CADJ_lrank} is > \code{max_lrank} will be pruned. 
#' @return None
#' @usage TRNobj$prune_CADJ_lrank(max_lrank)
NULL 

#' @name prune_CADJ_grank
#' @title Prune CADJ by a maximum global rank 
#' @description Prune the CADJ graph by removing edges whose \code{CADJ_grank} is \strong{strictly greater than} a given maximum grank.
#' @param max_grank maximum \code{CADJ_grank} allowed in the pruned graph. 
#' Any edges whose \code{CADJ_grank} is > \code{max_grank} will be pruned. 
#' @return None
#' @usage TRNobj$prune_CADJ_grank(max_grank)
NULL 

#' @name prune_CADJ_length
#' @title Prune CADJ by a maximum output topology length
#' @description Prune the CADJ graph by removing edges whose \code{CADJ_fwdflen} is \strong{stricly greater than} a given maximim forward folding length.
#' @param max_fwdflen maximum \code{CADJ_fwdflen} allowed in the pruned graph. 
#' Any edges whose \code{CADJ_fwdflen} is > \code{max_fwdflen} will be pruned. 
#' @return None
#' @usage TRNobj$prune_CADJ_length(max_fwdflen)
NULL 

#' @name prune_CADJ_wdeg
#' @title Prune CADJ by a minimum vertex weight 
#' @description Prune the CADJ graph by removing edges incident to vertices whose \code{CADJ_wdeg} is \strong{stricly less than} a given minimum weighted degree.
#' @param min_weight minimum vertex weight allowed in the pruned graph. 
#' Any vertices whose \code{CADJ_wdeg} is < \code{min_weight} will be pruned. 
#' @return None
#' @usage TRNobj$prune_CADJ_wdeg(min_weight)
NULL 

#' @name get_CADJ 
#' @title Return the CADJ adjacency matrix
#' @description This method returns the CADJ adjacency matrix of the TRN (as set during \code{initialize_TRN}) whose \code{CADJ_active} flag = 1. 
#' @return an adjacency matrix (nrows = ncols = \code{nV}) whose (i,j) entries are the \code{CADJ} edge weights. 
#' @usage TRNobj$get_CADJ()
NULL 

#' @name restore_CADJ_edges
#' @title Restore all CADJ edges
#' @description When pruning of CADJ edges occurs (via any of the \code{prune_CADJ_*} methods), the \code{CADJ_active} flag corresponding 
#' to the pruned edges is changed from 1 to 0. This method restores all \code{CADJ_active} flags to 1 (reversing any existing pruning). 
#' @return None
#' @usage TRNobj$restore_CADJ_edges()
NULL 


#' @name prune_CONN_edge_list
#' @title Prune CONN by an edge list 
#' @description Prune the CONN graph given a list of edges to remove. 
#' @param edge_list a 2-column matrix rows define the vertices connected by each CONN edge to prune 
#' @details The vertex indices in \code{edge_list} should be 1-based. 
#' @return None
#' @usage TRNobj$prune_CONN_edge_list(edge_list)
NULL 


#' @name prune_CONN_edge_id
#' @title Prune CONN by edge ids
#' @description Prune the CONN graph given a list of edge IDs to remove.
#' @param edge_id a vector containing the rows indices of \code{CONN_EL} which contain the edges desired to prune. 
#' @return None
#' @usage TRNobj$prune_CONN_edge_id(edge_id)
NULL 

#' @name prune_CONN_vertex_id
#' @title Prune CONN by a set of vertex ids
#' @description Prune the CONN graph by removing edges connecting vertices in a list of given vertex IDs.
#' @param vertex_id a vector containing the vertices to prune. 
#' All connections in \code{CONN_EL} from or to the vertices in \code{vertex_id} will be pruned. 
#' @return None
#' @usage TRNobj$prune_CONN_vertex_id(vertex_id)
NULL 

#' @name prune_CONN_edge_weight
#' @title Prune CONN by a minimum edge weight
#' @description Prune the CONN graph by removing edges whose weights are \strong{strictly less than} a given minimum edge weight.
#' @param min_weight minimum weight allowed in the pruned graph. 
#' Any edges in \code{CONN_EL} whose \code{CONN} value is < \code{min_weight} will be pruned. 
#' @return None
#' @usage TRNobj$prune_CONN_edge_weight(min_weight)
NULL 

#' @name prune_CONN_lrank
#' @title Prune CONN by a maximum local rank
#' @description Prune the CONN graph by removing edges whose \code{CONN_lrank} is \strong{strictly greater than} a given maximum lrank.
#' @param max_lrank maximum \code{CONN_lrank} allowed in the pruned graph. 
#' Any edges whose \code{CONN_lrank} is > \code{max_lrank} will be pruned. 
#' @return None
#' @usage TRNobj$prune_CONN_lrank(max_lrank)
NULL 

#' @name prune_CONN_grank
#' @title Prune CONN by a maximum global rank 
#' @description Prune the CONN graph by removing edges whose \code{CONN_grank} is \strong{strictly greater than} a given maximum grank.
#' @param max_grank maximum \code{CONN_grank} allowed in the pruned graph. 
#' Any edges whose \code{CONN_grank} is > \code{max_grank} will be pruned. 
#' @return None
#' @usage TRNobj$prune_CONN_grank(max_grank)
NULL 

#' @name prune_CONN_length
#' @title Prune CONN by a maximum output topology length
#' @description Prune the CONN graph by removing edges whose \code{CONN_fwdflen} is \strong{stricly greater than} a given maximim forward folding length.
#' @param max_fwdflen maximum \code{CONN_fwdflen} allowed in the pruned graph. 
#' Any edges whose \code{CONN_fwdflen} is > \code{max_fwdflen} will be pruned. 
#' @return None
#' @usage TRNobj$prune_CONN_length(max_fwdflen)
NULL 

#' @name prune_CONN_wdeg
#' @title Prune CONN by a minimum vertex weight 
#' @description Prune the CONN graph by removing edges incident to vertices whose \code{CONN_wdeg} is \strong{stricly less than} a given minimum weighted degree.
#' @param min_weight minimum vertex weight allowed in the pruned graph. 
#' Any vertices whose \code{CONN_wdeg} is < \code{min_weight} will be pruned. 
#' @return None
#' @usage TRNobj$prune_CONN_wdeg(min_weight)
NULL 

#' @name prune_CONN_CADJ
#' @title Prune CONN by CADJ 
#' @description Prune the CONN graph by removing CADJ inactive edges. 
#' @return None
#' @usage TRNobj$prune_CONN_CADJ()
NULL 

#' @name get_CONN 
#' @title Return the CONN adjacency matrix
#' @description This method returns the CONN adjacency matrix of the TRN (as set during \code{initialize_TRN}) whose \code{CONN_active} flag = 1. 
#' @return an adjacency matrix (nrows = ncols = \code{nV}) whose (i,j) entries are the \code{CONN} edge weights. 
#' @usage TRNobj$get_CONN()
NULL 

#' @name restore_CONN_edges
#' @title Restore all CONN edges
#' @description When pruning of CONN edges occurs (via any of the \code{prune_CONN_*} methods), the \code{CONN_active} flag corresponding 
#' to the pruned edges is changed from 1 to 0. This method restores all \code{CONN_active} flags to 1 (reversing any existing pruning). 
#' @return None
#' @usage TRNobj$restore_CONN_edges()
NULL 


#' @name get_OTADJ 
#' @title Return the TRN's output space adjacency matrix
#' @description This method returns the adjacency matrix of the output space topology of the TRN 
#' (as set during \code{initialize_TRN}), 
#' @return an adjacency matrix (nrows = ncols = \code{nV})  
#' @usage TRNobj$get_OTADJ()
NULL 

#' @name CADJ_active_verts
#' @title Return the CADJ active vertices
#' @description Active CADJ vertices participate in active CADJ edges; 
#' Thus if there are vertices in the graph with no active CADJ edges connecting them 
#' they are deemed inactive. 
#' @return a list of vertex IDs that are connected by any \code{CADJ_active} edge. 
#' @usage TRNobj$CADJ_active_verts()
NULL 

#' @name CONN_active_verts
#' @title Return the CONN active vertices
#' @description Active CONN vertices participate in active CONN edges; 
#' Thus if there are vertices in the graph with no active CONN edges connecting them 
#' they are deemed inactive. 
#' @return a list of vertex IDs that are connected by any \code{CONN_active} edge. 
#' @usage TRNobj$CONN_active_verts()
NULL 

#' @name calc_TopoMeasures
#' @title Calculate Topology Preserving Measures of the TRN 
#' @description Sets the Topology Preserving Measures of both the CADJ and CONN graphs (Topographic Product and several Topographic Functions). 
#' See \code{?Topographic_Product} and \code{?Topographic_Functions} for more details. 
#' 
#' @usage TRNobj$calc_TopoMeasures
#' 
#' \insertRef{BauerPawelzikGeisel1992}{TopoRNet}
#' \insertRef{VillmannDerHerrmannMartinetz1997}{TopoRNet}
#' \insertRef{ZhangMerenyi2006}{TopoRNet}
NULL 


# ***** DM Prune ***** 

#' @name set_BootSig
#' @title Set the CADJ significances 
#' @description Set the Bootstrap significance values of each CADJ and CONN edge.
#' @param BSADJ an adjacency matrix (nrows = ncols = \code{nV}) whose (i,j) entries contain significance values 
#' of each CADJ edge. 
#' 
#' @usage TRNobj$set_BootSig(BSADJ)
NULL 

#' @name calc_DMPrune_LambdaPath
#' @title Calculate the DM-Prune Lambda Path 
#' @description The DM-Prune Lambda path, at each prune step \eqn{t}, is the (normalized) Dirichlet-Multinomial likelihood of CADJ edges 
#' after pruning edges whose \code{DMP_mu_rank} is \strong{strictly less than} \eqn{t}.  
#' @param priorADJ an adjacency matrix (nrows = ncols = \code{nV}) whose (i,j) entries contain the Dirichlet prior CADJ values. 
#' 
#' @usage TRNobj$calc_DMPrune_LambdaPath(priorADJ)
NULL 


#' @name DMPrune_CADJ_step
#' @title Prune CADJ at a DM-Prune step 
#' @description  Prune the CADJ graph by removing edges in \code{CADJ_EL} whose \code{DMP_mu_rank} is \strong{strictly less than} a given minimum prune step.
#' @param min_step 
#' 
#' @usage TRNobj$DMPrune_CADJ_step(min_step)
NULL


# ***** Saving / Loading ***** 

#' @name save 
#' @title Save a TRN object 
#' @description All fields in a TRN object can be saved to disk with this function, which allows them to be re-loaded 
#' into a new R environment at a later time.  
#' @param trnfile a string indicating the file path and name in which to save the TRN object.  
#' This must end in extension ".trn", otherwise an error is returned. 
#' @details The TRN object is saved to disk as an R list, with each field occupying a corresponding field of the list.  
#' The file is saved in .rds format (can check its details with \link[base]{infoRDS}).  
#' 
#' Saved TRNs can be re-loaded with \link{load}
#' 
#' @return None, the TRN object is saved to disk 
#' @usage TRNobj$save(trnfile)
NULL 

#' @name load
#' @title Load an existing TRN object 
#' @description TRN objects previously written to disk via the \link{save} method can be re-loaded into 
#' a new R environment with this function. All fields of the internal C++ class will be populated, and 
#' all methods can be called on the loaded TRN object.  
#' @param trnfile a string indicating the file path and name of the saved TRN object.  
#' @details Because the .trn file is in .rds format it can, technically, be loaded directly into an R environment 
#' as a list via \link[base]{readRDS}. This can be useful for spot checking the contents of a saved TRN object, but 
#' does not allow use of any of its methods (or visualizations).  The \code{load} method allows for 
#' proper restoration of a previously saved TRN 
#' 
#' @return None, the TRN object is loaded 
#' @usage TRNobj$load(TRNfile)
NULL 

#' @name as_list
#' @title Convert an TRN object to a list 
#' @description This method extracts all fields of an TRN object and places their stored values 
#' into the fields of an R list, with list field names matching the TRN object field names. 
#' @return A list
#' @usage TRNobj$as_list()
NULL 

#' @name load_list 
#' @title Populate a TRN object from a list 
#' @description This method populates all fields of a TRN object from the fields of an R list object. 
#' The list must have field names which exactly match the TRN field names. 
#' @param TRNList a TRN object converted to a list, e.g., with \code{as_list}
#' @return None
#' @usage TRNobj$load_list(TRNList)
NULL 












