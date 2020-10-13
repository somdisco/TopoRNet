#' Visualize the DM-Prune Lambda Path 
#' 
#' @param TRN a TRN object 
#' @param min_step the minimum pruning step to show 
#' @param max_step the maximum pruning step to show 
#' @param scale a scale factor controlling the size of plotted lines, 
#' points, and text labels. Higher values increase the size. 
#' @param plot.cpts whether to compute and plot the variance changepoints of the 
#' second-order differenced Lambda path. Default = TRUE. 
#' 
#' @return None, a plot is produced. 
#' @export
vis_DMPrune_LambdaPath = function(TRN, min_step = NULL, max_step = NULL, scale = 1, plot.cpts = T) {
  ## Make sure the Lambda path exists 
  if(!(class(TRN)=="Rcpp_TRN")) stop("Input TRN must be a valid TRN object constructed from TRN$new()")
  if(!TRN$isset_DMPrune) stop("Must call $calc_DMPrune_LambdaPath first")
  
  ## Plot colors 
  red = "#d62728" 
  lightred = "#ffb09c"
  blue = "#1F77B4"
  lightblue = "#73C2FB"
  green = "#2ca02c" 
  lightgreen = "#98fb98"
  gray = "#88807b"
  lightgray = "#c7c6c1"
  
  library(ggplot2)
  
  
  cat(sprintf("!! Lambda Path Visualization !!\n"))
  cat(sprintf("----------------------------------------------------------\n"))
  
  
  ## Setup plot dataframe for Lambdas
  plotdf = data.frame(step = TRN$DMP_pruneStep, Lambda = TRN$DMP_Lambda)
  plotdf$dLambda = c(NA, diff(plotdf$Lambda))
  plotdf$d2Lambda = c(NA, NA, diff(diff(plotdf$Lambda)))
  
  ## Work out the special points
  L0 = TRN$DMP_Lambda[1]
  Lmax = Inf; Lintersect = Inf; 
  # if(max(TRN$DMP_Lambda > L0)) {
  #   Lmax = max(TRN$DMP_Lambda)
  #   stepmax = which.max(TRN$DMP_Lambda)
  #   cat(sprintf("max(Lambda) = %0.2g found at step = %d\n", Lmax, stepmax))
  #   
  #   stepintersect = stepmax + which.min(abs(TRN$DMP_Lambda[(stepmax+1):length(TRN$DMP_Lambda)] - L0))
  #   Lintersect = TRN$DMP_Lambda[stepintersect]
  #   cat(sprintf("Lambda path intersects its starting value at step = %d\n", stepintersect))
  # }
  stepmax = which.max(plotdf$Lambda)
  Lmax = plotdf$Lambda[stepmax]
  cat(sprintf("max(Lambda) = %0.2g found at step = %d\n", Lmax, stepmax))
  if(Lmax > L0) {
    stepintersect = stepmax + which.min(abs(plotdf$Lambda[(stepmax+1):nrow(plotdf)] - L0))
    Lintersect = plotdf$Lambda[stepintersect]
    cat(sprintf("Lambda path intersects its starting value at step = %d\n", stepintersect))
  }
  
  if(is.null(min_step)) min_step = 0
  if(is.null(max_step)) max_step = max(plotdf$step)
  
  ## COMPUTE CHANGEPOINTS HERE! 
  if(plot.cpts) {
    cpts_idx = changepoint::cpts(changepoint::cpt.var(data = plotdf$d2Lambda[3:nrow(plotdf)], method = "PELT", penalty = "MBIC"))
    cpts_idx = cpts_idx + 2 # +2 because we threw away the 1st two NA values in d2Lambda when computing cpts 
    cptsdf = data.frame(step = plotdf$step[cpts_idx], Lambda = plotdf$Lambda[cpts_idx], dLambda = plotdf$dLambda[cpts_idx], d2Lambda = plotdf$d2Lambda[cpts_idx])
    
    cat(sprintf("Change points identified at steps: \n"))
    for(i in 1:nrow(cptsdf)) {
      cat(sprintf("%d\t",cptsdf$step[i]))
      if(i %% 5 == 0 || i == nrow(cptsdf)) cat(sprintf("\n"))
    }
  }
  cat(sprintf("----------------------------------------------------------\n"))
  
  ## Restrict data frames to requested plotting range 
  plotdf = subset(plotdf, step >= min_step & step <= max_step)
  cptsdf = subset(cptsdf, step >= min_step & step <= max_step)
  
  
  ## *** Build Lambda plot 
  pL = ggplot(data = plotdf) + 
    geom_hline(yintercept = L0, size = 0.35*scale, linetype = "dashed", color = lightred) + 
    geom_line(aes(x = step, y = Lambda), size = 1*scale, color = lightblue)
  
  # Add point at max, if max > starting L0 
  if(is.finite(Lmax) && (stepmax >= min_step && stepmax <= max_step)) {
    tmp_plotdf = data.frame(step = stepmax, Lambda = Lmax)
    pL = pL + ggrepel::geom_text_repel(data = tmp_plotdf, aes(x=step, y=Lambda, label=step), color = blue, 
                                       size = 3*scale, fontface = 2, hjust=1, box.padding = 1, segment.size = 0.2, segment.color = gray)
    pL = pL + geom_point(data = tmp_plotdf, aes(x = step, y = Lambda), size = 2*scale, color = blue)
    rm(tmp_plotdf)
  }
  
  # Add point where Lambda intersects its starting values, after its max, if it exists 
  if(is.finite(Lintersect) && (stepintersect >= min_step && stepintersect <= max_step)) {
    tmp_plotdf = data.frame(step = stepintersect, Lambda = Lintersect)
    pL = pL + ggrepel::geom_text_repel(data = tmp_plotdf, aes(x=step, y=Lambda, label=step), color = red, 
                                       size = 3*scale, fontface = 2, hjust=1, box.padding = 1, segment.size = 0.2, segment.color = gray)
    pL = pL + geom_point(data = tmp_plotdf, aes(x = step, y = Lambda), size = 2*scale, color = red)
    rm(tmp_plotdf)
  }
  
  # Add changepoints 
  if(plot.cpts) {
    pL = pL + ggrepel::geom_text_repel(data = cptsdf, aes(x=step, y=Lambda, label=step), color = green, direction = 'both', 
                                       size = 3*scale, fontface = 2, hjust=1, box.padding = 1, segment.size = 0.2, segment.color = gray)
    pL = pL + geom_point(data = cptsdf, aes(x = step, y = Lambda), size = 2*scale, color = green)
  }
  
  
  pL = pL + xlab("step") + ylab(NULL) + ggtitle(latex2exp::TeX("DM-Prune $\\Lambda(t)$")) + theme_linedraw()
  
  
  ## *** Build dLambda plot 
  pdL = ggplot(data = plotdf) + 
    geom_line(aes(x = step, y = dLambda), size = 0.9*scale, color = lightblue)  
  
 
  # Add changepoints 
  if(plot.cpts) {
    pdL = pdL + ggrepel::geom_text_repel(data = cptsdf, aes(x=step, y=dLambda, label=step), color = green, direction = 'both', 
                                         size = 3*scale, fontface = 2, hjust=1, box.padding = 1, segment.size = 0.2, segment.color = gray)
    pdL = pdL + geom_point(data = cptsdf, aes(x = step, y = dLambda), size = 2*scale, color = green)
  }
  
  pdL = pdL + xlab("step") + ylab(NULL) + ggtitle(latex2exp::TeX("$d\\Lambda(t)$")) + theme_linedraw()
  #ggpubr::theme_pubclean()
  
  
  ## *** Build d2Lambda plot 
  pd2L = ggplot(data = plotdf) + 
    geom_line(aes(x = step, y = d2Lambda), size = 0.9*scale, color = lightblue)
  
  # Add changepoints 
  if(plot.cpts) {
    pd2L = pd2L + ggrepel::geom_text_repel(data = cptsdf, aes(x=step, y=d2Lambda, label=step), color = green, direction = 'both', 
                                           size = 3*scale, fontface = 2, hjust=1, box.padding = 1, segment.size = 0.2, segment.color = gray)
    pd2L = pd2L + geom_point(data = cptsdf, aes(x = step, y = d2Lambda), size = 2*scale, color = green)
  }
  
  pd2L = pd2L + xlab("step") + ylab(NULL) + ggtitle(latex2exp::TeX("$d^2\\Lambda(t)$")) + theme_linedraw()
  
  
  ## Draw 
  gridExtra::grid.arrange(pL, pdL, pd2L, layout_matrix = matrix(c(1,1,2,3),ncol=2))
}



#' Visualize the DM-Prune prior 
#' 
#' @description Plots the CADJ prior values vs. the CADJ values for comparison. 
#' 
#' @param TRN a TRN object 
#' 
#' @return None, a plot is produced. 
#' @export
vis_DMPrune_prior = function(TRN) {
  if(!(class(TRN)=="Rcpp_TRN")) stop("Input TRN must be a valid TRN object constructed from TRN$new()")
  blue = "#1F77B4"
  plot(x = TRN$CADJ, y = TRN$DMP_prior, pch=16, col = blue, xlab = 'CADJ', ylab = 'Prior', main = 'DM-Prune Prior Comparison')
}
