#' flowCore: Basic structures for flow cytometry data
#' 
#' Provides S4 data structures and basic infrastructure and functions to deal
#' with flow cytometry data.
#' 
#' 
#' Define important flow cytometry data classes:
#' \code{\link[flowCore:flowFrame-class]{flowFrame}},
#' \code{\link[flowCore:flowSet-class]{flowSet}} and their accessors.
#' 
#' Provide important transformation, filter, gating, workflow, and summary
#' functions for flow cytometry data analysis.
#' 
#' Most of flow cytometry related Bioconductor packages (such as flowStats,
#' flowFP, flowQ, flowViz, flowMerge, flowClust) are heavily dependent on this
#' package.
#' 
#' \tabular{ll}{ Package: \tab flowCore \cr Type: \tab Package \cr Version:
#' \tab 1.11.20 \cr Date: \tab 2009-09-16 \cr License: \tab Artistic-2.0\cr }
#' 
#' @useDynLib flowCore, .registration=TRUE 
#' @import methods
#' @import Rcpp
#' @importFrom matrixStats colMins colMaxs
#' @importClassesFrom Biobase AnnotatedDataFrame
#' @importFrom grDevices devAskNewPage
#' @importFrom graphics abline
#' @importFrom stats density na.omit
#' @importFrom utils read.csv
#' @importMethodsFrom Biobase description "description<-" exprs "exprs<-"
#' featureNames pData "pData<-" phenoData "phenoData<-" sampleNames 
#' "sampleNames<-" varLabels "varLabels<-" varMetadata "varMetadata<-"
#' @importMethodsFrom graph addEdge addNode adj edgeData "edgeData<-"
#' "edgeDataDefaults<-" nodes removeNode
#' @importFrom Biobase listLen read.AnnotatedDataFrame
#' @importFrom BiocGenerics colnames eval get mget nrow ncol sort
#' @importFrom stats4 plot summary
#' @importFrom graph "edgeRenderInfo<-" "nodeRenderInfo<-"
#' @importFrom graphics hist text smoothScatter
#' @importFrom rrcov CovMcd
#' @importFrom MASS cov.rob
#' @importFrom stats kmeans mad median qnorm quantile runif sd uniroot var
#' @importFrom utils data find head tail write.table
#' @importFrom corpcor pseudoinverse
#' 
#' 
#' 
#' @name flowCore-package
#' @aliases flowCore flowCore-package
#' @docType package
#' @author 
#' Maintainer: Florian Hahne <fhahne@@fhcrc.org>
#' 
#' Authors: B. Ellis, P. Haaland, F. Hahne, N. Le Meur, N. Gopalakrishnan
#' @keywords package
NULL
