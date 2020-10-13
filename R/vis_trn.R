#' Visualize CADJ Statistics
#' 
#' @description The CADJvis statistics will be plotted by each local rank, global rank, or output topology distance, as requested. 
#' 
#' @param TRN a TRN object 
#' @param group_by which grouping to visualize, default = 'lrank'. Can also be 'grank' or 'length'
#' 
#' @export 
vis_CADJvis_stats = function(TRN, group_by = "lrank") {
  if(!(class(TRN)=="Rcpp_TRN")) stop("Input TRN must be a valid TRN object constructed from TRN$new()")
  
  ## Decode group by 
  if(!(group_by == "lrank" || group_by == "grank" || group_by == "length")) stop("group_by must be one of 'lrank','grank','length'")
  
  if(group_by=="lrank") {
    plotdf = TRN$get_CADJvis_stats("lrank")
    xlabstr = "Local Rank"
  } else if(group_by == "grank") {
    plotdf = TRN$get_CADJvis_stats("grank")
    xlabstr = "Global Rank"
  } else if(group_by == "length") {
    plotdf = TRN$get_CADJvis_stats("length")
    xlabstr = "Output Topology Length"
  }
  
  ## Add in summary statistics 
  plotdf$IQR = plotdf$q75 - plotdf$q25
  plotdf$whisker_lo = pmax(plotdf$q0, plotdf$q25-1.5*plotdf$IQR)
  plotdf$whisker_hi = pmin(plotdf$q100, plotdf$q75+1.5*plotdf$IQR)
  
  
  ## Plot defaults 
  titlestr = "CADJ Edge Weight Summary"
  boxwidth=0.9
  centers = plotdf[,1]
  
  ## Pre-defined colors
  tableaugray = "#88807b"
  tableaublue = "#1F77B4"
  tableaured = "#d62728"
  tableaulightteal = "#c8ffff"
  
  ## Build plot 
  plot(x=c(boxwidth/2,centers,max(centers)+boxwidth/2), y=c(0,plotdf$whisker_hi,0), cex=0, axes = F, xlab = xlabstr, ylab = NA)
  # Whiskers
  segments(x0 = centers, x1 = centers, y0 = plotdf$whisker_lo, y1 = plotdf$whisker_hi, col = tableaugray)
  # Box 
  rect(xleft = centers - boxwidth/2, xright = centers + boxwidth/2, ybottom = plotdf$q25, ytop = plotdf$q75, border = tableaublue, col = tableaulightteal, lwd = 1.2)
  # Median 
  segments(x0 = centers - boxwidth*0.45, x1 = centers + boxwidth*0.45, y0 = plotdf$q50, y1 = plotdf$q50, col = tableaured, lwd=1.5)
  axis(side = 1, at = centers)
  axis(side = 2, at = pretty(plotdf$hi))
  title(main = titlestr, line = -1)
  
}

#' Visualize CONN Statistics
#' 
#' @description The CONNvis statistics will be plotted by each local rank, global rank, or output topology distance, as requested. 
#' 
#' @param TRN a TRN object 
#' @param group_by which grouping to visualize, default = 'lrank'. Can also be 'grank' or 'length'
#' 
#' @export 
vis_CONNvis_stats = function(TRN, group_by = "lrank") {
  if(!(class(TRN)=="Rcpp_TRN")) stop("Input TRN must be a valid TRN object constructed from TRN$new()")
  
  ## Decode group by 
  if(!(group_by == "lrank" || group_by == "grank" || group_by == "length")) stop("group_by must be one of 'lrank','grank','length'")
  
  if(group_by=="lrank") {
    plotdf = TRN$get_CONNvis_stats("lrank")
    xlabstr = "Local Rank"
  } else if(group_by == "grank") {
    plotdf = TRN$get_CONNvis_stats("grank")
    xlabstr = "Global Rank"
  } else if(group_by == "length") {
    plotdf = TRN$get_CONNvis_stats("length")
    xlabstr = "Output Topology Length"
  }
  
  ## Add in summary statistics 
  plotdf$IQR = plotdf$q75 - plotdf$q25
  plotdf$whisker_lo = pmax(plotdf$q0, plotdf$q25-1.5*plotdf$IQR)
  plotdf$whisker_hi = pmin(plotdf$q100, plotdf$q75+1.5*plotdf$IQR)
  
  
  ## Plot defaults 
  titlestr = "CONN Edge Weight Summary"
  boxwidth=0.9
  centers = plotdf[,1]
  
  ## Pre-defined colors
  tableaugray = "#88807b"
  tableaublue = "#1F77B4"
  tableaured = "#d62728"
  tableaulightteal = "#c8ffff"
  
  ## Build plot 
  plot(x=c(boxwidth/2,centers,max(centers)+boxwidth/2), y=c(0,plotdf$whisker_hi,0), cex=0, axes = F, xlab = xlabstr, ylab = NA)
  # Whiskers
  segments(x0 = centers, x1 = centers, y0 = plotdf$whisker_lo, y1 = plotdf$whisker_hi, col = tableaugray)
  # Box 
  rect(xleft = centers - boxwidth/2, xright = centers + boxwidth/2, ybottom = plotdf$q25, ytop = plotdf$q75, border = tableaublue, col = tableaulightteal, lwd = 1.2)
  # Median 
  segments(x0 = centers - boxwidth*0.45, x1 = centers + boxwidth*0.45, y0 = plotdf$q50, y1 = plotdf$q50, col = tableaured, lwd=1.5)
  axis(side = 1, at = centers)
  axis(side = 2, at = pretty(plotdf$hi))
  title(main = titlestr, line = -1)
  
}




# # Return a vector of CONNvis edge colors, based on input vector of edge local rankings 
# .CONNvis_colors = function(lrank) {
#   ## Assign CONNvis colors
#   connvis_colors = c("red","blue","green","yellow")
#   if(max(lrank) > 4) {
#     connvis_colors = c(connvis_colors, gray.colors(max(0, max(lrank)-4)))  
#   } 
#   return(connvis_colors[lrank])
# }




#' CONNvis Visualization 
#' 
#' @param TRN a TRN object 
#' @param add whether to create a new plotting device (=FALSE, default), or add to an existing one (=TRUE)
#' @param vertex.xy a matrix of (x,y) coordinates (nrows = TRN$nV, ncols = 2) defining the plot coordinates 
#' of the vertices of the TRN. It is assumed that the \code{i-th} row of \code{vertex.xy} gives the coordinates for vertex \code{i} of the graph. 
#' @param vertex.pch the pch symbol plotted vertices, default = 16. 
#' @param vertex.cex the cex of plotted vertices, default = 1. Set = 0 to suppress vertex plotting. 
#' @param vertex.col the color of plotted vertices, default = 'black'
#' @param edge.lwd_range the min/max range of the plotted line widths. 
#' Default = NULL means line widths inherit from the \code{TRN$CONN_grank} of each edge. 
#' If supplied as something other than NULL, it must be a length=2 vector giving the (lower,upper) bounds of plotted edge widths. 
#' @param vertex.active whether to restrict plotted vertices to those which have an edge connecting them, default = TRUE.
#' @param vertex.subset a vector of vertex indices to restrict the plotting to. Default = NULL, meaning all edges and vertices are plotted. 
#' If given, any edges connecting vertices in this list (along with the vertices themselves) will not be plotted. 
#' 
#' @details Only active CONN edges are plotted (i.e., those with \code{TRN$CONN_active = 1}), 
#' so that any previously pruned CONN edges are not shown.  If all edges are desired, called \code{TRN$restore_CONN_edges} 
#' prior to plotting.  
#'
#' @export
vis_CONNvis = function(TRN, add = F, vertex.xy, 
                       vertex.pch = 16, vertex.cex=1, vertex.col = "black", 
                       edge.lwd_range = NULL, 
                       vertex.active = T, vertex.subset = NULL) {
  if(!(class(TRN)=="Rcpp_TRN")) stop("Input TRN must be a valid TRN object constructed from TRN$new()")
  
  vertex.xy = as.matrix(vertex.xy)
  stopifnot(ncol(vertex.xy)==2)
  stopifnot(nrow(vertex.xy)==TRN$nV)
  edge.lwd_range = sort(edge.lwd_range, decreasing = F)
  
  ## Setup plotting dataframe
  # Add vertices (x,y) locations and whether each edge is active 
  plotdf = data.frame(TRN$CONN_EL, vertex.xy[TRN$CONN_EL[,1],], vertex.xy[TRN$CONN_EL[,2],], TRN$CONN_active, TRN$CONN_grank, TRN$CONN_lrank)
  names(plotdf) = c('v1','v2','x1','y1','x2','y2', 'active', 'grank', 'lrank')
  # Add the colors 
  plotdf$color = TRN$get_CONNvis_colors()
  # Add the line widths, scale if requested 
  plotdf$width = TRN$get_CONNvis_widths()
  if(!is.null(edge.lwd_range)) {
    if(length(edge.lwd_range) != 2) stop("length(edge.lwd_range) must = 2 giving lower / upper bounds for edge widths")
    edge.lwd_range = sort(edge.lwd_range)
    plotdf$width = (plotdf$width - min(plotdf$width)) / diff(range(plotdf$width)) * diff(edge.lwd_range) + edge.lwd_range[1]
  }
  
  ## Sort the data frame by descending global & local rank
  ## The first sort ensures thin lines are overplotted, the second ensures red edges are always on top 
  plotdf = dplyr::arrange(plotdf, dplyr::desc(grank), dplyr::desc(lrank))
  
  ## Toss out any inactive edges 
  plotdf = subset(plotdf, active==1)
  
  ## Subset the vertices 
  ## If given as null, then vertex.subset is reset to all vertices
  ## Otherwise, it is assumed to contain vertex indices to restrict plotting to
  if(is.null(vertex.subset)) {
    vertex.subset = 1:TRN$nV
  } else {
    if(!all(vertex.subset %in% 1:TRN$nV)) stop("vertex.subset contains out of range vertex indices")
  }
  ## Further remove any inactive vertices, if requested 
  if(vertex.active)
    vertex.subset = intersect(vertex.subset, TRN$CONN_active_verts())
  
  plotdf = subset(plotdf, v1 %in% vertex.subset & v2 %in% vertex.subset)
  
  ## Temporarily change line end types to "butt" 
  oldlend = par()$lend 
  par(lend = "butt")
  
  ## Open a plot window, if requested 
  if(!add)
    plot(x = vertex.xy[vertex.subset,1,drop=F], y = vertex.xy[vertex.subset,2,drop=F], cex=0, col='black', xlab = 'x', ylab = 'y', frame.plot = F)
  
  ## Plot the edges 
  with(plotdf, segments(x0 = x1, y0 = y1, x1 = x2, y1 = y2, col = color, lty = 1, lwd = width)) 
  
  if(vertex.cex > 0)
    points(x = vertex.xy[vertex.subset,1,drop=F], y = vertex.xy[vertex.subset,2,drop=F], pch = vertex.pch, cex = vertex.cex)
  
  par(lend = oldlend)
}



#' CADJvis Visualization 
#' 
#' @description In CADJvis, plotted edges only extend to the midpoint between source and sink vertices 
#' in order to highlight the asymmetries in the CADJ graph.  The colors and line widths are computed 
#' relative to CADJ (not CONN).  
#' 
#' @param TRN a TRN object 
#' @param add whether to create a new plotting device (=FALSE, default), or add to an existing one (=TRUE)
#' @param vertex.xy a matrix of (x,y) coordinates (nrows = TRN$nV, ncols = 2) defining the plot coordinates 
#' of the vertices of the TRN. It is assumed that the \code{i-th} row of \code{vertex.xy} gives the coordinates for vertex \code{i} of the graph. 
#' @param vertex.pch the pch symbol plotted vertices, default = 16. 
#' @param vertex.cex the cex of plotted vertices, default = 1. Set = 0 to suppress vertex plotting. 
#' @param vertex.col the color of plotted vertices, default = 'black'
#' @param edge.lwd_range the min/max range of the plotted line widths. 
#' Default = NULL means line widths inherit from the \code{TRN$CADJ_grank} of each edge. 
#' If supplied as something other than NULL, it must be a length=2 vector giving the (lower,upper) bounds of plotted edge widths. 
#' @param vertex.active whether to restrict plotted vertices to those which have an edge connecting them, default = TRUE.
#' @param vertex.subset a vector of vertex indices to restrict the plotting to. Default = NULL, meaning all edges and vertices are plotted. 
#' If given, any edges connecting vertices in this list (along with the vertices themselves) will not be plotted. 
#' 
#' @details Only active CADJ edges are plotted (i.e., those with \code{TRN$CADJ_active = 1}), 
#' so that any previously pruned CADJ edges are not shown.  If all edges are desired, called \code{TRN$restore_CADJ_edges} 
#' prior to plotting.  
#'
#' @export
vis_CADJvis = function(TRN, add = F, vertex.xy, 
                       vertex.pch = 16, vertex.cex=1, vertex.col = "black", 
                       edge.lwd_range = NULL, 
                       vertex.active = T, vertex.subset = NULL) {
  if(!(class(TRN)=="Rcpp_TRN")) stop("Input TRN must be a valid TRN object constructed from TRN$new()")
  
  vertex.xy = as.matrix(vertex.xy)
  stopifnot(ncol(vertex.xy)==2)
  stopifnot(nrow(vertex.xy)==TRN$nV)
  edge.lwd_range = sort(edge.lwd_range, decreasing = F)
  
  ## Setup plotting dataframe
  # Add vertices (x,y) locations and whether each edge is active 
  plotdf = data.frame(TRN$CADJ_EL, vertex.xy[TRN$CADJ_EL[,1],], vertex.xy[TRN$CADJ_EL[,2],], TRN$CADJ_active, TRN$CADJ_grank, TRN$CADJ_lrank)
  names(plotdf) = c('v1','v2','x1','y1','x2','y2', 'active', 'grank', 'lrank')
  # Modify the ending portion of the line segments to be the midpoint 
  plotdf = dplyr::mutate(plotdf, x2 = (x1 + x2)/2, y2 = (y1 + y2)/2)
  #with(plotdf, x2 = (x1 + x2)/2, y2 = (y1 + y2)/2)
  # Add the colors 
  plotdf$color = TRN$get_CADJvis_colors()
  # Add the line widths, scale if requested 
  plotdf$width = TRN$get_CADJvis_widths()
  if(!is.null(edge.lwd_range)) {
    if(length(edge.lwd_range) != 2) stop("length(edge.lwd_range) must = 2 giving lower / upper bounds for edge widths")
    edge.lwd_range = sort(edge.lwd_range)
    plotdf$width = (plotdf$width - min(plotdf$width)) / diff(range(plotdf$width)) * diff(edge.lwd_range) + edge.lwd_range[1]
  }
  
  ## Sort the data frame by descending global & local rank
  ## The first sort ensures thin lines are overplotted, the second ensures red edges are always on top 
  plotdf = dplyr::arrange(plotdf, dplyr::desc(grank), dplyr::desc(lrank))
  
  ## Toss out any inactive edges 
  plotdf = subset(plotdf, active==1)
  
  ## Subset the vertices 
  ## If given as null, then vertex.subset is reset to all vertices
  ## Otherwise, it is assumed to contain vertex indices to restrict plotting to
  if(is.null(vertex.subset)) {
    vertex.subset = 1:TRN$nV
  } else {
    if(!all(vertex.subset %in% 1:TRN$nV)) stop("vertex.subset contains out of range vertex indices")
  }
  ## Further remove any inactive vertices, if requested 
  if(vertex.active)
    vertex.subset = intersect(vertex.subset, TRN$CONN_active_verts())
  
  plotdf = subset(plotdf, v1 %in% vertex.subset & v2 %in% vertex.subset)
  
  ## Temporarily change line end types to "butt" 
  oldlend = par()$lend 
  par(lend = "butt")
  
  ## Open a plot window, if requested 
  if(!add)
    plot(x = vertex.xy[vertex.subset,1,drop=F], y = vertex.xy[vertex.subset,2,drop=F], cex=0, col='black', xlab = 'x', ylab = 'y', frame.plot = F)
  
  ## Plot the edges 
  with(plotdf, segments(x0 = x1, y0 = y1, x1 = x2, y1 = y2, col = color, lty = 1, lwd = width)) 
  
  if(vertex.cex > 0)
    points(x = vertex.xy[vertex.subset,1,drop=F], y = vertex.xy[vertex.subset,2,drop=F], pch = vertex.pch, cex = vertex.cex)
  
  par(lend = oldlend)
}



#' TopoView Visualization 
#' 
#' @param ADJ a (possibly weighted) adjacency matrix of a Topology Representing Network.   
#' @param add whether to create a new plotting device (=FALSE, default), or add to an existing one (=TRUE)
#' @param vertex.xy a matrix of (x,y) coordinates (nrows = TRN$nV, ncols = 2) defining the plot coordinates of the vertices of the TRN. 
#' It is assumed that the \code{i-th} row of \code{vertex.xy} gives the coordinates for vertex \code{i} of the graph. 
#' @param vertex.pch the pch symbol plotted vertices, default = 16. 
#' @param vertex.cex the cex of plotted vertices, default = 1. Set = 0 to suppress vertex plotting. 
#' @param vertex.col the color of plotted vertices. 
#' @param edge.color line color of plotted edges, default = "darkorange".
#' @param edge.lwd_range the min/max range of the plotted line widths, Default = c(1, 5). 
#' This parameter is only valid if input \code{ADJ} is weighted, in which case the line widths represent 
#' the edge weights via a linear scaling from \code{(min weight, max weight)} to \code{edge.lwd_range}. 
#' If \code{ADJ} is unweighted, the plotted edges have width = \code{max(edge.lwd_range)}. 
#' @param vertex.active whether to restrict plotted vertices to those which have an edge connecting them, default = TRUE.
#' @param vertex.subset a vector of vertex indices to restrict the plotting to. Default = NULL, meaning all edges and vertices are plotted. 
#' If given, any edges connected vertices in this list (along with the vertices themselves) will not be plotted. 
#' 
#' @details It is assumed that any 0 value in the input \eqn{ADJ_{ij}} means no edge connects vertices \eqn{i} and \eqn{j}. 
#'
#' @export
vis_TopoView = function(ADJ, add = F, vertex.xy, 
                         vertex.pch = 16, vertex.cex=1, vertex.col = "black", 
                         edge.col = "darkorange", edge.lwd_range = c(1,5), 
                         vertex.active = T, vertex.subset = NULL) {
  
  vertex.xy = as.matrix(vertex.xy)
  stopifnot(ncol(vertex.xy)==2)
  stopifnot(nrow(vertex.xy)==nrow(ADJ))
  stopifnot(nrow(ADJ)==ncol(ADJ))
  stopifnot(length(edge.lwd_range)==2)
  edge.lwd_range = sort(edge.lwd_range, decreasing = F)
  
  
  ## Setup plotting dataframe 
  EL = which(ADJ > 0, arr.ind = T)
  plotdf = data.frame(EL, vertex.xy[EL[,1],], vertex.xy[EL[,2],], ADJ[EL])
  names(plotdf) = c('v1', 'v2', 'x1', 'y1', 'x2', 'y2', 'w')
  
  plotdf$color = edge.col
  
  ## Assign line widths 
  w.range = diff(range(plotdf$w))
  if(w.range > 0) {
    plotdf$width = (plotdf$w - min(plotdf$w)) / diff(range(plotdf$w)) * diff(edge.lwd_range) + edge.lwd_range[1]  
  } else {
    plotdf$width = max(edge.lwd_range)
  }
  
  
  
  ## Subset the vertices 
  ## If given as null, then vertex.subset is reset to all vertices
  ## Otherwise, it is assumed to contain vertex indices to restrict plotting to
  if(is.null(vertex.subset)) {
    vertex.subset = 1:nrow(ADJ)
  } else {
    if(!all(vertex.subset %in% 1:nrow(ADJ))) stop("vertex.subset contains out of range vertex indices")
  }
  plotdf = subset(plotdf, v1 %in% vertex.subset & v2 %in% vertex.subset)
  ## Further remove any inactive vertices, if requested 
  if(vertex.active) {
    active_verts = which(apply(ADJ,1,sum) > 0)
    vertex.subset = intersect(vertex.subset, active_verts)
  }
    
  ## Temporarily change line end types to "butt" 
  oldlend = par()$lend 
  par(lend = "butt")
  
  ## Open a plot window, if requested 
  if(!add)
    plot(x = vertex.xy[vertex.subset,1,drop=F], y = vertex.xy[vertex.subset,2,drop=F], cex=0, col='black', xlab = 'x', ylab = 'y', frame.plot = F)
  
  ## Plot the edges 
  with(plotdf, segments(x0 = x1, y0 = y1, x1 = x2, y1 = y2, col = color, lty = 1, lwd = width))  
  
  if(vertex.cex > 0)
    points(x = vertex.xy[vertex.subset,1,drop=F], y = vertex.xy[vertex.subset,2,drop=F], pch = vertex.pch, cex = vertex.cex)
  
  par(lend = oldlend)
}


#' Visualize the Topographic Functions for CADJ 
#' 
#' @param TRN a TRN object
#' 
#' @details Three subplots are returned, corresponding to the Topographic Function, the Discrete Topographic Function, and 
#' the Weighted Discrete Topographic Function, respectively. 
#' See \code{?Topographic_Functions} for further details. 
#' 
#' @export
vis_CADJ_TopoFunctions = function(TRN) {
  if(!(class(TRN)=="Rcpp_TRN")) stop("Input TRN must be a valid TRN object constructed from TRN$new()")
  if(!TRN$isset_TPM) stop("Must call $calc_TopoMeasures before calling vis_TopoFunctions")
  
  lightblue = "#73C2FB"
  blue = "#1F77B4"
  lightgreen = "#98fb98"
  green = "#2ca02c"
  lightorange = "#ffd27f"
  orange = "#ff7f0e"
  
  oldmfrow = par()$mfrow
  
  par(mfrow=c(1,3))
  barplot(height = TRN$CADJ_TopoFxns$TF, names.arg = TRN$CADJ_TopoFxns$k, main = 'TF', col = lightblue, border = blue)
  barplot(height = TRN$CADJ_TopoFxns$DTF, names.arg = TRN$CADJ_TopoFxns$k, main = 'DTF', col = lightgreen, border = green)
  barplot(height = TRN$CADJ_TopoFxns$WDTF, names.arg = TRN$CADJ_TopoFxns$k, main = 'WDTF', col = lightorange, border = orange)
  par(mfrow=oldmfrow)
  
}


#' Visualize the Topographic Functions for CONN 
#' 
#' @param TRN a TRN object
#' 
#' @details Three subplots are returned, corresponding to the Topographic Function, the Discrete Topographic Function, and 
#' the Weighted Discrete Topographic Function, respectively. 
#' See \code{?Topographic_Functions} for further details. 
#' 
#' @export
vis_CONN_TopoFunctions = function(TRN) {
  if(!(class(TRN)=="Rcpp_TRN")) stop("Input TRN must be a valid TRN object constructed from TRN$new()")
  if(!TRN$isset_TPM) stop("Must call $calc_TopoMeasures before calling vis_TopoFunctions")
  
  lightblue = "#73C2FB"
  blue = "#1F77B4"
  lightgreen = "#98fb98"
  green = "#2ca02c"
  lightorange = "#ffd27f"
  orange = "#ff7f0e"
  
  oldmfrow = par()$mfrow
  
  par(mfrow=c(1,3))
  barplot(height = TRN$CONN_TopoFxns$TF, names.arg = TRN$CONN_TopoFxns$k, main = 'TF', col = lightblue, border = blue)
  barplot(height = TRN$CONN_TopoFxns$DTF, names.arg = TRN$CONN_TopoFxns$k, main = 'DTF', col = lightgreen, border = green)
  barplot(height = TRN$CONN_TopoFxns$WDTF, names.arg = TRN$CONN_TopoFxns$k, main = 'WDTF', col = lightorange, border = orange)
  par(mfrow=oldmfrow)
  
}




