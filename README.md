# TopoRNet

`TopoRNet` is an R package providing a C++ template for manipulation of and computation with a Topology Representing Network (or TRN, which is also introduced in this document).  As TRNs represent mappings from an **input space** to an **output space** they possess an additional layer of complexity when compared to standard mathematical graph representations.  While several R packages for graph analysis currently exist (the [`igraph`](https://igraph.org/r/) package being the most popular) none (inherently) contain the machinery required to properly handle this additional complexity.  `TopoRNet` addresses this shortcoming in the R community by providing the following functionality which, to our knowledge, no other package does: 

+ A templated C++ class to facilitate storage of TRNs and derivative products.
+ Calculation of the **Topographic Product** and several **Topographic Functions.** Collectively, these Topology Preserving Measures (TPMs) provide a metric for assessing the *degree* to which topology is preserved, when it is represented by the TRN.
+ The **CONNvis** visualization (and calculation of its required statistics) which aides human and automated cluster inference from TRNs. 
+ Intelligent sparsification of a TRN (which improves the quality of cluster extraction) via its CONNvis statistics and the recently introduced `DM-Prune` methodology.  
+ C++ implementations of the above (based on [Rcpp](https://cran.r-project.org/web/packages/Rcpp/index.html) and [RcppArmadillo](https://cran.r-project.org/web/packages/RcppArmadillo/index.html))
+ Optional parallel computation of the above (as applicable, via [RcppParallel](https://cran.r-project.org/web/packages/RcppParallel/index.html))

# Installation 

```r
devtools::install_github("somdisco/TopoRNet")
```

# Documentation 

See the [TopoRNet homepage](https://somdisco.github.io/TopoRNet/output/index.html) for more information.
