% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/lspartition-package.R
\docType{package}
\name{lspartition-package}
\alias{lspartition-package}
\title{Nonparametric Estimation and Inference using Partitioning-Based Least Squares Regression}
\description{
This package provides tools for statistical analysis using B-splines, wavelets, and
            piecewise polynomials as described in
            \href{https://arxiv.org/abs/1804.04916}{Cattaneo, Farrell and Feng (2019a)}:
            \code{\link{lsprobust}} for least squares point estimation with robust bias-corrected pointwise and
            uniform inference procedures; \code{\link{lspkselect}} for data-driven procedures
            for selecting the IMSE-optimal number of partitioning knots; \code{\link{lsprobust.plot}}
            for regression plots with robust confidence intervals and confidence bands;
            \code{\link{lsplincom}} for estimation and inference for linear combination of regression
            functions of different groups.

            The companion software article, \href{https://arxiv.org/abs/1906.00202}{Cattaneo, Farrell and Feng (2019b)},
            provides further implementation details and empirical illustrations.
}
\references{
Cattaneo, M. D., M. H. Farrell, and Y. Feng (2019a): \href{https://arxiv.org/abs/1804.04916}{Large Sample Properties of Partitioning-Based Series Estimators}. Annals of Statistics, forthcoming. arXiv:1804.04916.

Cattaneo, M. D., M. H. Farrell, and Y. Feng (2019b): \href{https://arxiv.org/abs/1906.00202}{lspartition: Partitioning-Based Least Squares Regression}. R Journal, forthcoming. arXiv:1906.00202.
}
\author{
Matias D. Cattaneo, Princeton University, Princeton, NJ. \email{cattaneo@princeton.edu}.

       Max H. Farrell, University of Chicago, Chicago, IL. \email{max.farrell@chicagobooth.edu}.

       Yingjie Feng (maintainer), Princeton University, Princeton, NJ. \email{yingjief@princeton.edu}.
}
