#'Tuning Parameter Selection Procedures for Partitioning-Based Regression Estimation and Inference
#'
#'@description \code{lspkselect} implements data-driven procedures to select the Integrated Mean Squared Error (IMSE) optimal number of partitioning knots for partitioning-based least squares regression estimators. Three series methods are supported: B-splines, compactly supported wavelets, and piecewise polynomials.
#'             See \href{https://sites.google.com/site/nppackages/lspartition/Cattaneo-Farrell_2013_JoE.pdf?attredirects=0}{Cattaneo and Farrell (2013)} and \href{https://arxiv.org/abs/1804.04916}{Cattaneo, Farrell and Feng (2019a)} for complete details.
#'
#'             Companion commands: \code{\link{lsprobust}} for partitioning-based least squares regression estimation and inference; \code{\link{lsprobust.plot}} for plotting results; \code{\link{lsplincom}} for multiple sample estimation and inference.
#'
#'             A detailed introduction to this command is given in \href{https://arxiv.org/abs/1906.00202}{Cattaneo, Farrell and Feng (2019b)}.
#'
#'             For more details, and related Stata and R packages useful for empirical analysis,
#'             visit \url{https://sites.google.com/site/nppackages/}.
#'
#'@param y Outcome variable.
#'@param x Independent variable. A matrix or data frame.
#'@param m Order of basis used in the main regression. Default is \code{m=2}.
#'@param m.bc Order of basis used to estimate leading bias. Default is \code{m.bc=m+1}.
#'@param smooth Smoothness of B-splines for point estimation. When \code{smooth=s}, B-splines have \code{s}-order
#'              continuous derivatives. Default is \code{smooth=m-2}.
#'@param bsmooth Smoothness of B-splines for bias correction. Default is \code{bsmooth=m.bc-2}.
#'@param deriv Derivative order of the regression function to be estimated. A vector object of the same
#'             length as \code{ncol(x)}. Default is \code{deriv=c(0,...,0)}.
#'@param method Type of basis used for expansion. Options are \code{"bs"} for B-splines,
#'              \code{"wav"} for compactly supported wavelets (Cohen, Daubechies and Vial, 1993),
#'              and \code{"pp"} for piecewise polynomials. Default is \code{method="bs"}.
#'@param bc Bias correction method. Options are \code{"bc1"} for higher-order-basis bias correction,
#'          \code{"bc2"} for least squares bias correction, and \code{"bc3"} for plug-in bias correction.
#'          Defaults are \code{"bc3"} for splines and piecewise polynomials and \code{"bc2"}
#'          for wavelets.
#'@param ktype Knot placement. Options are \code{"uni"} for evenly spaced knots over the
#'             support of \code{x} and \code{"qua"} for quantile-spaced knots. Default is \code{ktype="uni"}.
#'@param kselect Method for selecting the number of inner knots used by \code{lspkselect}. Options
#'               are \code{"imse-rot"} for a rule-of-thumb (ROT) implementation of IMSE-optimal number of knots,
#'               \code{"imse-dpi"} for second generation direct plug-in (DPI) implementation of IMSE-optimal number
#'               of knots, and \code{"all"} for both. Default is \code{kselect="imse-dpi"}.
#'@param proj If \code{TRUE}, projection of leading approximation error onto the lower-order approximating space
#'            is included for bias correction (splines and piecewise polynomial only). Default is \code{proj=TRUE}.
#'@param vce Procedure to compute the heteroskedasticity-consistent (HCk) variance-covariance matrix estimator with plug-in residuals. Options are
#'           \itemize{
#'           \item \code{"hc0"} for unweighted residuals (HC0).
#'           \item \code{"hc1"} for HC1 weights.
#'           \item \code{"hc2"} for HC2 weights. Default.
#'           \item \code{"hc3"} for HC3 weights.
#'           }
#'@param subset Optional rule specifying a subset of observations to be used.
#'@param rotnorm If \code{TRUE}, ROT selection is adjusted using normal densities.
#'@return \item{\code{ks}}{A matrix may contain \code{k.rot} (IMSE-optimal number of knots for
#'                         the main regression through ROT implementation), \code{k.bias.rot}
#'                         (IMSE-optimal number of knots for bias correction through ROT
#'                         implementation), \code{k.dpi} (IMSE-optimal number of knots for the
#'                         main regression through DPI implementation), \code{k.bias.dpi} (IMSE-optimal
#'                         number of knots for bias correction through DPI implementation)}
#'        \item{\code{opt}}{A list containing options passed to the function.}
#'@author
#' Matias D. Cattaneo, Princeton University, Princeton, NJ. \email{cattaneo@princeton.edu}.
#'
#' Max H. Farrell, University of Chicago, Chicago, IL. \email{max.farrell@chicagobooth.edu}.
#'
#' Yingjie Feng (maintainer), Princeton University, Princeton, NJ. \email{yingjief@princeton.edu}.
#'
#'@references
#'
#' Cattaneo, M. D., and M. H. Farrell (2013): \href{https://sites.google.com/site/nppackages/lspartition/Cattaneo-Farrell_2013_JoE.pdf?attredirects=0}{Optimal convergence rates, Bahadur representation, and asymptotic normality of partitioning estimators}. Journal of Econometrics 174(2): 127-143.
#'
#' Cattaneo, M. D., M. H. Farrell, and Y. Feng (2019a): \href{https://arxiv.org/abs/1804.04916}{Large Sample Properties of Partitioning-Based Series Estimators}. Annals of Statistics, forthcoming. arXiv:1804.04916.
#'
#' Cattaneo, M. D., M. H. Farrell, and Y. Feng (2019b): \href{https://arxiv.org/abs/1906.00202}{lspartition: Partitioning-Based Least Squares Regression}. R Journal, forthcoming. arXiv:1906.00202.
#'
#' Cohen, A., I. Daubechies, and P.Vial (1993): Wavelets on the Interval and Fast Wavelet Transforms. Applied and Computational Harmonic Analysis 1(1): 54-81.
#'
#'@seealso \code{\link{lsprobust}}, \code{\link{lsprobust.plot}}, \code{\link{lsplincom}}
#'
#'@examples
#'x   <- data.frame(runif(500), runif(500))
#'y   <- sin(4*x[,1])+cos(x[,2])+rnorm(500)
#'est <- lspkselect(y, x)
#'summary(est)
#'
#'@export
#'
# version 0.4 Aug2019
lspkselect = function(y, x, m=NULL, m.bc=NULL, smooth=NULL, bsmooth=NULL, deriv=NULL, method="bs", ktype="uni",
                      kselect="imse-dpi", proj = TRUE, bc="bc3", vce="hc2", subset=NULL, rotnorm=TRUE) {

  if (!is.null(m)) m <- m - 1
  q <- m.bc
  if (!is.null(q)) q <- m.bc - 1

  x <- as.matrix(x)
  if (is.data.frame(y)) y <- y[,1]
  # substract subset
  if (!is.null(subset)) {
    x <- x[subset, , drop = F]
    y <- y[subset]
  }

  na.ok <- complete.cases(x) & complete.cases(y)

  x <- x[na.ok, , drop = F]
  y <- y[na.ok]

  method   <- tolower(method)
  ktype    <- tolower(ktype)
  kselect  <- tolower(kselect)
  vce      <- tolower(vce)

  method.type <- "B-spline"
  if (method == "localpoly" | method == "pp") method.type <- "Piecewise Polynomial"
  if (method == "wavelet" | method == "wav")  method.type <- "Wavelet"

  knot.type <- ""
  if (ktype == "uniform" | ktype =="uni")    knot.type <- "Uniform"
  if (ktype == "quantile" | ktype =="qua")   knot.type <- "Quantile"
  if (method == "wav" | method == "wavelet") knot.type <- "Uniform"
  #############################################
  #############################################

  if (method == "bspline")   method <- "bs"
  if (method == "localpoly") method <- "pp"
  if (method == "wavelet")   method <- "wav"

  if (ktype == "uniform")  ktype <- "uni"
  if (ktype == "quantile") ktype <- "qua"

  N <- nrow(x); d <- ncol(x)

  # order of basis
  if (is.null(deriv) | method == "wav")  deriv <- rep(0, d)
  if (method == "bs" | method == "pp") {
     if (is.null(m))        m <- max(deriv) + 1
     if (is.null(q))        q <- m + 1
  } else if (method == "wav") {
     if (is.null(m)) m <- max(deriv) + 1
     # q <- m + 1
     if (bc == "bc3") bc <- "bc2"
  }

  if (method == "bs") {
    if (is.null(smooth)) {
      smooth <- m - 1
    } else {
      if (smooth < m-1)   proj <- TRUE
      if (m < smooth + 1) m <- smooth + 1
    }
    if (is.null(bsmooth)) {
      bsmooth <- q - 1
    } else {
      if (bsmooth < q-1)   proj <- TRUE
      if (q < bsmooth + 1) q <- bsmooth + 1
    }
  }

  deriv.str <- paste0("(", toString(deriv), ")")

  if (kselect == "imse-rot") {
    ks <- matrix(NA, 1, 2)
    colnames(ks) <- c("k.rot", "k.bias.rot")
    ks[1,1] <- lspkselect.imse.rot(y, x, m, method, rotnorm, ktype, deriv)
    if (bc == "bc3") {
      deriv.bias <- genIndex(m, d, method)
      k.bias.seq <- apply(deriv.bias, 1, function(z)
                          lspkselect.imse.rot(y, x, q, method, rotnorm, ktype, deriv=z))
      ks[1,2] <- ceiling(mean(k.bias.seq))
    }
  } else if (kselect == "imse-dpi") {
    ks <- matrix(NA, 1, 2)
    colnames(ks) <- c("k.dpi", "k.bias.dpi")
    ks[1,1] <- lspkselect.imse.dpi(y, x, m, method, ktype, vce, deriv, proj=proj, smooth=smooth, rotnorm)
    if (bc == "bc3") {
      deriv.bias <- genIndex(m, d, method)
      k.bias.seq <- apply(deriv.bias, 1, function(z)
                          lspkselect.imse.dpi(y, x, q, method, ktype, vce, deriv=z,
                                              proj=proj, smooth=bsmooth, rotnorm=rotnorm))
      ks[1,2] <- ceiling(mean(k.bias.seq))
    }
  } else if (kselect == "all") {
    ks <- matrix(NA, 1, 4)
    colnames(ks) <- c("k.rot", "k.bias.rot", "k.dpi", "k.bias.dpi")
    ks[1,1] <- lspkselect.imse.rot(y, x, m, method, rotnorm, ktype, deriv)
    ks[1,3] <- lspkselect.imse.dpi(y, x, m, method, ktype, vce, deriv, proj=proj, smooth=smooth, rotnorm=rotnorm)
    if  (bc == "bc3") {
        deriv.bias <- genIndex(m, d, method)
        k.bias.seq <- apply(deriv.bias, 1, function(z)
                            lspkselect.imse.rot(y, x, q, method, rotnorm, ktype, deriv=z))
        ks[1,2] <- ceiling(mean(k.bias.seq))

        k.bias.seq <- apply(deriv.bias, 1, function(z)
                            lspkselect.imse.dpi(y, x, q, method, ktype, vce, deriv=z,
                                                proj=proj, smooth=bsmooth, rotnorm=rotnorm))
        ks[1,4] <- ceiling(mean(k.bias.seq))
    }
  }

  out        <- list(ks = ks, opt = list(n=N, m=m+1, m.bc=q+1, deriv=deriv.str, ktype=knot.type,
                                         d=d, method=method.type, kselect=kselect))
  out$call   <- match.call()
  class(out) <- "lspkselect"
  return(out)
}

#'@describeIn lspkselect \code{print} method for class "\code{lspkselect}".
#'@param ... further arguments
#'@export
print.lspkselect <- function(x, ...) {
  cat("Call: lspkselect\n\n")

  cat(paste("Sample size (n)                            =    ", x$opt$n,        "\n", sep=""))
  cat(paste("Basis function (method)                    =    ", x$opt$method,   "\n", sep=""))
  cat(paste("Order of basis point estimation (m)        =    ", x$opt$m,        "\n", sep=""))
  cat(paste("Order of derivative (deriv)                =    ", x$opt$deriv,    "\n", sep=""))
  cat(paste("Order of basis bias correction (m.bc)      =    ", x$opt$m.bc,     "\n", sep=""))
  cat(paste("Knot placement (ktype)                     =    ", x$opt$ktype,    "\n", sep=""))
  cat(paste("Knot method  (kselect)                     =    ", x$opt$kselect,  "\n", sep=""))
  cat("\n")
}

#'@describeIn  lspkselect \code{summary} method for class "\code{lspkselect}".
#'@param object class \code{lspkselect} objects.
#'@export
summary.lspkselect <- function(object,...) {
  x <- object
  args <- list(...)
  if (is.null(args[['sep']]))   { sep <- 5 } else { sep <- args[['sep']] }

  cat("Call: lspkselect\n\n")

  cat(paste("Sample size (n)                            =    ", x$opt$n,        "\n", sep=""))
  cat(paste("Basis function (method)                    =    ", x$opt$method,   "\n", sep=""))
  cat(paste("Order of basis point estimation (m)        =    ", x$opt$m,        "\n", sep=""))
  cat(paste("Order of derivative (deriv)                =    ", x$opt$deriv,    "\n", sep=""))
  cat(paste("Order of basis bias correction (m.bc)      =    ", x$opt$m.bc,     "\n", sep=""))
  cat(paste("Knot placement (ktype)                     =    ", x$opt$ktype,    "\n", sep=""))
  cat(paste("Knot method  (kselect)                     =    ", x$opt$kselect,  "\n", sep=""))
  cat("\n")

  if (x$opt$kselect=="all") {
    col1.names = c("", "IMSE-ROT","", "IMSE-DPI")
    col2.names = rep(c("k", "k.bc"), 2)
  } else if (x$opt$kselect=="imse-rot") {
    col1.names = c("", "IMSE-ROT")
    col2.names = c("k", "k.bc")
  } else {
    col1.names = c("", "IMSE-DPI")
    col2.names = c("k", "k.bc")
  }

  ### print output
  if (x$opt$kselect=="imse-dpi" | x$opt$kselect=="imse-rot") {
    cat(paste(rep("=", 15 + 8), collapse="")); cat("\n")
  } else {
    cat(paste(rep("=", 15 + 8*ncol(x$ks)), collapse="")); cat("\n")
  }

  cat(format(col1.names, width=8, justify="right"))
  cat("\n")

  cat(format(col2.names, width=8, justify="right"))
  cat("\n")

  if (x$opt$kselect=="imse-dpi" | x$opt$kselect=="imse-rot") {
    cat(paste(rep("=", 15 + 8), collapse="")); cat("\n")
  } else {
    cat(paste(rep("=", 15 + 8*ncol(x$ks)), collapse="")); cat("\n")
  }

  if (x$opt$kselect=="imse-dpi" | x$opt$kselect=="imse-rot") {
    cat(format(sprintf("%3.0f", x$ks[1:2])  , width=8, justify="right"))
    cat("\n")
    cat(paste(rep("=", 15 + 8), collapse="")); cat("\n")
  } else {
    cat(format(sprintf("%3.0f", x$ks[1:4])  , width=8, justify="right"))
    cat("\n")
    cat(paste(rep("=", 15 + 8*ncol(x$ks)), collapse="")); cat("\n")
  }
}
