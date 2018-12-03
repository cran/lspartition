#'Nonparametric Estimation and Inference using Partitioning-Based Least Squares Regression
#'@description This package provides tools for statistical analysis using B-splines, wavelets, and
#'             piecewise polynomials (generalized regressogram) as described in
#'             \href{https://sites.google.com/site/nppackages/lspartition/Cattaneo-Farrell-Feng_2018_Partitioning.pdf?attredirects=0}{Cattaneo, Farrell and Feng (2018a)}.
#'             \code{\link{lsprobust}} for least squares point estimation with robust bias-corrected pointwise and
#'             uniform inference procedures. \code{\link{lspkselect}} for data-driven procedures
#'             for selecting the IMSE-optimal number of partitioning knots. \code{\link{lsprobust.plot}}
#'             for regression plots with robust confidence intervals and confidence bands.
#'             \code{\link{lsplincom}} for estimation and inference for linear combination of regression
#'             functions of different groups.
#'
#'             The companion software article, \href{https://sites.google.com/site/nppackages/lspartition/Cattaneo-Farrell-Feng_2018_lspartition.pdf?attredirects=0}{Cattaneo, Farrell and Feng (2018b)},
#'             provides further implementation details and empirical illustration.
#'
#'@importFrom stats lm sd complete.cases quantile qnorm poly rbinom rnorm
#'@importFrom matrixStats rowProds colProds colMaxs colMins rowMins rowMaxs
#'@importFrom combinat xsimplex
#'@importFrom splines splineDesign
#'@importFrom pracma bernoulli
#'@importFrom mgcv tensor.prod.model.matrix
#@importFrom graphics lines plot polygon
#'@importFrom MASS ginv
#@importFrom dplyr distinct
#'@import ggplot2
#'@docType package
#'@name lspartition-package
#'@author Matias D. Cattaneo, University of Michigan, Ann Arbor, MI. \email{ cattaneo@umich.edu}.
#'
#'        Max H. Farrell, University of Chicago, Chicago, IL. \email{max.farrell@chicagobooth.edu}.
#'
#'        Yingjie Feng, University of Michigan, Ann Arbor, MI. \email{yjfeng@umich.edu}.
#'
#'@references
#' Cattaneo, M. D., M. H. Farrell, and Y. Feng (2018a): \href{https://sites.google.com/site/nppackages/lspartition/Cattaneo-Farrell-Feng_2018_Partitioning.pdf?attredirects=0}{Large Sample Properties of Partitioning-Based Series Estimators}. Working paper.
#'
#' Cattaneo, M. D., M. H. Farrell, and Y. Feng (2018b): \href{https://sites.google.com/site/nppackages/lspartition/Cattaneo-Farrell-Feng_2018_lspartition.pdf?attredirects=0}{lspartition: Partitioning-Based Least Squares Regression}. Working paper.
#'
NULL
