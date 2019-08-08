#'Partitioning-Based Least Squares Regression with Robust Inference.
#'
#'@description \code{lsprobust} implements partitioning-based least squares point estimators for the regression function and its derivatives. It also provides robust bias-corrected (pointwise and uniform) inference, including simulation-based confidence bands. Three series methods are supported: B-splines, compact supported wavelets, and piecewise polynomials.
#'             See \href{https://sites.google.com/site/nppackages/lspartition/Cattaneo-Farrell_2013_JoE.pdf?attredirects=0}{Cattaneo and Farrell (2013)} and \href{https://arxiv.org/abs/1804.04916}{Cattaneo, Farrell and Feng (2019a)} for complete details.
#'
#'             Companion commands: \code{\link{lspkselect}} for data-driven IMSE-optimal selection of the number of knots on rectangular partitions; \code{\link{lsprobust.plot}} for plotting results; \code{\link{lsplincom}} for multiple sample estimation and inference.
#'
#'             A detailed introduction to this command is given in \href{https://arxiv.org/abs/1906.00202}{Cattaneo, Farrell and Feng (2019b)}.
#'
#'             For more details, and related Stata and R packages useful for empirical analysis,
#'             visit \url{https://sites.google.com/site/nppackages/}.
#'
#'@param y Outcome variable.
#'@param x Independent variable. A matrix or data frame.
#'@param eval Evaluation points. A matrix or data frame.
#'@param neval Number of quantile-spaced evaluating points.
#'@param method Type of basis used for expansion. Options are \code{"bs"} for B-splines,
#'              \code{"wav"} for compactly supported wavelets (Cohen, Daubechies and Vial, 1993),
#'              and \code{"pp"} for piecewise polynomials. Default is \code{method="bs"}.
#'@param m Order of basis used in the main regression. Default is \code{m=2}. For B-splines,
#'         if \code{smooth} is specified but \code{m} is unspecified, default is \code{m=smooth+2}.
#'@param m.bc Order of basis used to estimate leading bias. Default is \code{m.bc=m+1}. For B-splines,
#'            if \code{bsmooth} is specified but \code{m.bc} is unspecified, default is \code{m.bc=bsmooth+2}.
#'@param deriv Derivative order of the regression function to be estimated. A vector object of the same
#'             length as \code{ncol(x)}. Default is \code{deriv=c(0,...,0)}.
#'@param smooth Smoothness of B-splines for point estimation. When \code{smooth=s}, B-splines have \code{s}-order
#'              continuous derivatives. Default is \code{smooth=m-2}.
#'@param bsmooth Smoothness of B-splines for bias correction. Default is \code{bsmooth=m.bc-2}.
#'@param ktype Knot placement. Options are \code{"uni"} for evenly-spaced knots over the
#'             support of \code{x} and \code{"qua"} for quantile-spaced knots. Default is \code{ktype="uni"}.
#'@param knot A list of numeric vectors giving the knot positions (including boundary knots) for each dimension
#'            which are used in the main regression. The length of the list is equal to \code{ncol(x)}.
#'            If not specified, it uses the number of knots either specified by users
#'            or computed by the companion command \code{lspkselect} to generate the
#'            corresponding knots according to the rule specified by \code{ktype}.
#'@param nknot A numeric vector of the same length as \code{ncol(x)}. Each element corresponds to
#'             the number of \emph{inner} partitioning knots for each dimension used in the main regression.
#'             If not specified, \code{nknot} is computed by the companion command \code{lspkselect}.
#'@param same If \code{TRUE}, the same knots are used for bias correction as that for the
#'            main regression. Default is \code{same=TRUE}.
#'@param bknot A list of numeric vectors giving knot positions used for bias correction. If not
#'             specified and \code{same=FALSE}, it uses the number of knots either specified by
#'             users or computed by the companion command \code{lspkselect} to generate
#'             knots according to the rule specified by \code{ktype}.
#'@param bc Bias correction method. Options are \code{"bc1"} for higher-order-basis bias correction,
#'          \code{"bc2"} for least squares bias correction, and \code{"bc3"} for plug-in bias correction.
#'          Default are \code{"bc3"} for splines and piecewise polynomials and \code{"bc2"}
#'          for wavelets.
#'@param proj If \code{TRUE}, projection of leading approximation error onto the lower-order approximation space
#'            is included for bias correction (splines and piecewise polynomials only). Default is
#'            \code{proj=TRUE}.
#'@param bnknot A numeric vector of the same length as \code{ncol(x)}. Each element corresponds
#'              to the number of \emph{inner} partitioning knots for each dimension used for bias
#'              correction. If not specified, \code{bnknot} is computed by the companion command \code{lspkselect}.
#'@param J A numeric vector containing resolution levels of father wavelets for each dimension.
#'@param kselect Method for selecting the number of \emph{inner} knots used by \code{lspkselect}. Options
#'               are \code{"imse-rot"} for ROT implementation of IMSE-optimal number of knots and
#'               \code{"imse-dpi"} for second generation of DPI implementation of IMSE-optimal number
#'               of knots. Default is \code{kselect="imse-dpi"}.
#'@param vce Procedure to compute the heteroskedasticity-consistent (HCk) variance-covariance matrix estimator with plug-in residuals. Options are
#'           \itemize{
#'           \item \code{"hc0"} for unweighted residuals (HC0).
#'           \item \code{"hc1"} for HC1 weights.
#'           \item \code{"hc2"} for HC2 weights. Default.
#'           \item \code{"hc3"} for HC3 weights.
#'           }
#'@param level Confidence level used for confidence intervals; default is \code{level=95}.
#'@param uni.method Method used to implement uniform inference. Options are \code{"pl"} for
#'                  a simulation-based plug-in procedure, \code{"wb"} for a wild bootstrap
#'                  procedure. If unspecified, neither procedure is
#'                  implemented. Default is \code{uni.method=NULL}.
#'@param uni.grid  A matrix containing all grid points used to implement uniform inference. Each row
#'                 correponds to the coordinates of one grid point.
#'@param uni.ngrid A numeric vector of the same length as \code{ncol(x)}. Each element corresponds to
#'                 the number of grid points for each dimension used to implement uniform inference.
#'                 Default is \code{uni.ngrid=50}.
#'@param uni.out If \code{TRUE}, the quantities used to implement uniform inference is outputted. Default is
#'               \code{uni.out=FALSE}.
#'@param band If \code{TRUE}, the critical value for constructing confidence band is calculated. Default
#'            is \code{band=FALSE}. If \code{band=TRUE} with \code{uni.method} unspecified,
#'            default is \code{uni.method="pl"}.
#'@param B    Number of simulated samples used to obtain the critical value for confidence bands.
#'            Default is \code{B=1000}.
#'@param subset Optional rule specifying a subset of observations to be used.
#'@param rotnorm If \code{TRUE}, ROT selection is adjusted using normal densities.
#'@return \item{\code{Estimate}}{ A matrix containing eval (grid points), N (effective sample sizes),
#'                             tau.cl (point estimates with a basis of order \code{m}), tau.bc (bias corrected point
#'                             estimates with a basis of order \code{m.bc}), se.cl (standard error corresponding
#'                             to tau.cl), and se.rb (robust standard error).}
#'        \item{\code{k.num}}{ A matrix containing the number of inner partitioning knots used in the main
#'                             regression and bias correction for each covariate.}
#'        \item{\code{knot}}{ A list of knots for point estimation.}
#'        \item{\code{bknot}}{ A list of knots for bias correction.}
#'        \item{\code{sup.cval}}{ Critical value for constructing confidence band.}
#'        \item{\code{uni.output}}{ A list containing quantities used to implement uniform inference.}
#'        \item{\code{opt}}{ A list containing options passed to the function.}
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
#'@seealso \code{\link{lspkselect}}, \code{\link{lsprobust.plot}}, \code{\link{lsplincom}}
#'
#'@examples
#'x   <- data.frame(runif(500), runif(500))
#'y   <- sin(4*x[,1])+cos(x[,2])+rnorm(500)
#'est <- lsprobust(y, x)
#'summary(est)
#'
#'@export

# Version 0.4 Aug2019
lsprobust = function(y, x, eval=NULL, neval=NULL, method="bs", m=NULL, m.bc=NULL, deriv=NULL, smooth=NULL,
                     bsmooth=NULL, ktype="uni", knot=NULL, nknot=NULL, same = TRUE, bknot=NULL, bnknot=NULL,
                     J=NULL, bc="bc3", proj=TRUE, kselect="imse-dpi", vce="hc2", level=95, uni.method=NULL,
                     uni.grid=NULL, uni.ngrid=50, uni.out=FALSE, band=FALSE, B=1000, subset=NULL, rotnorm=TRUE) {

  if (!is.null(m)) m <- m - 1
  q <- m.bc
  if (!is.null(q)) q <- q - 1

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

  x.max <- colMaxs(x); x.min <- colMins(x)
  d <- ncol(x)

  if (!is.null(eval) & d == 1) eval <- as.matrix(eval)

  if (d == 1 & is.vector(knot, mode="numeric")) {
    knot <- list(knot)
  }

  if (d == 1 & is.vector(bknot, mode="numeric")) {
    bknot <- list(bknot)
  }

  # drop x outside boundary
  if (is.list(knot)) {
    knot <- lapply(knot, sort)
    knot.min <- sapply(knot, function(z) z[1])
    knot.max <- sapply(knot, function(z) z[length(z)])
    if (any(c(x.min < knot.min, x.max > knot.max))) {
      warning("x outside boundary knots dropped")
      x.ind <- rowSums((sweep(x, 2, knot.min) >= 0) + (sweep(x, 2, knot.max) <= 0)) == 2*d
      x <- x[x.ind,, drop=F]
      y <- y[x.ind]
    }
  }

  if (is.list(bknot)) {
    bknot <- lapply(bknot, sort)
    bknot.min <- sapply(bknot, function(x) x[1])
    bknot.max <- sapply(bknot, function(x) x[length(x)])
    if (any(c(x.min < bknot.min, x.max > bknot.max))) {
      warning("x outside boundary knots dropped")
      x.ind <- rowSums((sweep(x, 2, bknot.min) >= 0) + (sweep(x, 2, bknot.max) <= 0)) == 2*d
      x <- x[x.ind,, drop=F]
      y <- y[x.ind]
    }
  }
  N <- nrow(x); x.max <- colMaxs(x); x.min <- colMins(x)

  method   <- tolower(method)
  ktype    <- tolower(ktype)
  kselect  <- tolower(kselect)
  vce      <- tolower(vce)

  ######################################################   CHECK ERRORS
  exit <- 0
  if (!is.null(eval)) {
    if (ncol(eval) != d) {
      print("evaluating points incorrectly")
      exit <- 1
    }
  }

  if (method!="bs" & method!="bspline" & method!="pp" &
      method!="localpoly" & method!="wav" & method!="wavelet" & method!="" ){
    print("method incorrectly specified")
    exit <- 1
  }

  if  (kselect!="imse-dpi" & kselect!="imse-rot" & kselect!=""){
    print("kselect incorrectly specified")
    exit <- 1
  }

  if  (ktype!="uni" & ktype!="uniform" & ktype!="qua" & ktype!="quantile" & ktype!=""){
    print("ktype incorrectly specified")
    exit <- 1
  }

  if (vce!="" & vce!="hc1" & vce!="hc2" & vce!="hc3" & vce!="hc0"){
    print("vce incorrectly specified")
    exit <- 1
  }

  if(!is.null(m)) {
    if (m < 0) {
      print("m incorrectly specified")
      exit <- 1
    }
  }

  if ((!is.null(m) | !is.null(q)) & (method == "wav" | method == "wavelet")) {
     if (m > 3 | q > 3) {
       print("wavelet of order greater than 4 not allowed")
       exit <- 1
     }
  }

  if (d > 1 & !is.null(knot) & !is.list(knot)) {
    print("knot is not a list")
    exist <- 1
  }

  if (d > 1 & !is.null(bknot) & !is.list(bknot)) {
    print("bknot is not a list")
    exist <- 1
  }

  if (!is.null(knot) & (method == "wav" | method == "wavelet")) {
    print("knot option not allowed for wavelets; try J or nknot")
    exit <- 1
  }

  if ((method == "wav" | method == "wavelet") & !is.null(deriv)) {
    if (!all(deriv == rep(0, d))) {
      print("deriv estimates not allowed for wavelets")
      exit <- 1
    }
  }

  if(!is.null(m) & !is.null(smooth)) {
    if (m < smooth+1) {
      print("smoothness incorrectly specified")
      exit <- 1
    }
  }

  if(!is.null(q) & !is.null(bsmooth)) {
    if (q < bsmooth+1) {
      print("smoothness incorrectly specified")
      exit <- 1
    }
  }

  if (!is.null(smooth) & method != "bs" & method != "bspline") {
    print("smoothness incorrectly specified")
    exit <- 1
  }

  if (!is.null(bsmooth) & method != "bs" & method != "bspline") {
    print("smoothness incorrectly specified")
    exit <- 1
  }

  if(!is.null(q)) {
    if (q < 0) {
      print("m.bc incorrectly specified")
      exit <- 1
    }
  }

  if(!is.null(deriv)) {
    if (!all(deriv >= 0)) {
      print("deriv should be positive integers")
      exit <- 1
    }
  }

  if (!is.null(m) & !is.null(q)) {
    if (m >= q) {
       print("m.bc should be greater than m")
       exit <- 1
    }
  }

  if (!is.null(deriv) & !is.null(m)) {
     if (!all(deriv <= m)) {
        print("deriv cannot be greater than m")
        exit <- 1
     }
  }

  if (level>100 | level<=0) {
    print("level should be set between 0 and 100")
    exit <- 1
  }

  if (!bc %in% c("bc1", "bc2", "bc3", "none")) {
    print("correction method incorrectly specified")
    exit <- 1
  }

  if (!is.null(m) & !is.null(J)) {
    if (any(2^J < 2*(m+2))) {
       print("resolution level incorrectly specified")
       exit <- 1
    }
  }

  if (!is.null(nknot)) {
    if (any(nknot < 0)) {
      print("number of knots should be positive")
      exit <- 1
    }
  }

  if (!is.null(bnknot)) {
    if (any(bnknot < 0)) {
      print("number of knots should be positive")
      exit <- 1
    }
  }

  if (!is.null(uni.method)) {
    if (uni.method != "pl" & uni.method != "wb") {
      print("uniform method incorrectly specified")
      exit <- 1
    }
  }

  if (exit>0) stop()
  #####################################################
  # store option info
  if (method == "bspline")   method <- "bs"
  if (method == "localpoly") method <- "pp"
  if (method == "wavelet")   method <- "wav"

  method.type <- "B-spline"
  if (method == "pp")  method.type <- "Piecewise Polynomial"
  if (method == "wav") method.type <- "Wavelet"

  if (ktype == "uniform")    ktype  <- "uni"
  if (ktype == "quantile")   ktype  <- "qua"

  knot.type <- ""
  if (!is.null(knot))  knot.type <- "User-specified"
  if (ktype == "uni")  knot.type <- "Uniform"
  if (ktype == "qua")  knot.type <- "Quantile"
  if (method == "wav") knot.type <- "Uniform"

  kselect.type <- kselect
  if (!is.null(knot) | !is.null(nknot) | !is.null(J)) kselect.type <- "User-specified"

  smooth.p <- NA
  smooth.q <- NA

  #### some default options###############
  if (is.null(deriv))  deriv  <- rep(0, d)
  if (method != "wav") {
    if (is.null(m)) m <- max(deriv) + 1
    if (is.null(q)) q <- m + 1
  } else {
    deriv <- rep(0, d)
    if (is.null(m)) m <- 1
    if (is.null(q)) q <- m + 1
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
    smooth.p <- smooth
    smooth.q <- bsmooth
  }

  if (m == 0 & method == "wav") {
    method <- "bs"; ktype <- "uni"; smooth <- -1; bsmooth <- q - 1
    if (!is.null(J)) nknot <- 2^J-1
  }

  ############################################## Prepare knots!!!
  if (method != "wav") {
    if (is.null(knot)) {
      if (is.null(nknot)) {
         lsk <- lspkselect(y, x, m+1, q+1, smooth, bsmooth, deriv, method, ktype, kselect=kselect, proj=proj,
                           bc=bc, vce=vce, rotnorm=rotnorm)
         nknot <- rep(lsk$ks[,1], d)
      }
      if (ktype == "uni") {
           knot <- genKnot.u(x.min, x.max, d, nknot)
      } else if (ktype == "qua") {
           knot <- genKnot.q(x, d, nknot)
      }
    }

    if (bc != "bc3" | (same & is.null(bknot) & is.null(bnknot))) {
       bknot  <- knot
    }
    if (!same & is.null(bknot) & is.null(bnknot)) bnknot <- rep(lsk$ks[,2], d)
    if (is.null(bknot)) {
      if (ktype == "uni") {
        bknot <- genKnot.u(x.min, x.max, d, bnknot)
      } else if (ktype == "qua") {
        bknot <- genKnot.q(x, d, bnknot)
      }
    }
  } else {
    if (!is.null(nknot)) J <- pmax(ceiling(log2(nknot+1)), ceiling(log2(2*(m+2))))
    if (is.null(J)) {
      J <- lspkselect(y, x, m+1, q+1, smooth, bsmooth, deriv, method, ktype,
                      kselect=kselect, proj=proj, bc=bc, vce=vce, rotnorm=rotnorm)$ks[,1]
      J <- pmax(ceiling(log2(rep(J, d))), ceiling(log2(2*(m+2))))
    }
    knot  <- genKnot.u(x.min, x.max, d, rep(2^J-1, d))
    bknot <- knot
  }

  if (method != "wav") {
    k.p <- sapply(knot,  function(z) length(z)-2)
    k.q <- sapply(bknot, function(z) length(z)-2)
    k.num <- cbind(k.p, k.q)
  } else {
    k.num <- cbind(2^J-1, 2^J-1)
  }
  colnames(k.num) <- c("k", "k.bc")
  k.m.str <- paste0("(", toString(k.num[,1]), ")")
  k.b.str <- paste0("(", toString(k.num[,2]), ")")
  #####################################################

  # define evaluation points
  if (is.null(eval)) {
    if (is.null(neval)) {
      if (d == 1) {
        eval <- lapply(knot, function(z) cumsum(c(z[1], rep(diff(z)/11, each=11)))[-seq(1, by=11, length.out=length(z))])
        eval <- as.matrix(expand.grid(eval))
      } else {
        #eval <- unique(x)
        qseq <- seq(0, 1, 1/(20+1))
        eval <- apply(x, 2, function(z) quantile(z, qseq[2:(length(qseq)-1)], names=F))
      }
    } else {
      #eval <- seq(x.min,x.max,length.out=neval)
      qseq <- seq(0, 1, 1/(neval+1))
      eval <- apply(x, 2, function(z) quantile(z, qseq[2:(length(qseq)-1)], names=F))
      if (neval == 1) eval <- t(eval)
    }
  }

  eval  <- as.matrix(eval)
  neval <- nrow(eval)

  ##############################################
  deriv.str <- paste0("(", toString(deriv), ")")

  if (band & is.null(uni.method)) uni.method <- "pl"
  uni.method.str <- uni.method
  if (is.null(uni.method.str)) {
    uni.method.str <- NA
  } else if (uni.method.str == "pl") {
    uni.method.str <- "Plug-in"
  } else {
    uni.method.str <- "Bootstrap"
  }

  ################################
  # effective sample size
  if (method == "pp") {
    eN <- gen.eN(eval, x, d, knot)
  } else {
    eN <- N
  }

  #########################
  #Estimation and Inference
  #########################
  # classical estimate
  P <- genDesign(x, method, m, J, smooth, knot, rep(0, d))
  invG.p <- qrXXinv(P)
  basis.p <- genDesign(eval, method, m, J, smooth, knot, deriv)

  # define leverage
  hii.p <- rowSums((P %*% invG.p) * P)
  hii.p[hii.p >= 0.9] <- 0.9; hii.p[hii.p < 0] <- 0
  d.p <- ncol(P)

  # point estimate and inference
  beta.p <- invG.p %*% crossprod(P, y)
  tau.cl <- basis.p %*% beta.p
  predict.p <- P %*% beta.p
  res.p <- lsprobust.res(y, predict.p, hii.p, vce, d.p)
  se.cl <- sqrt(rowSums((basis.p %*% invG.p %*% lsprobust.vce(P, res.p) %*% invG.p) * basis.p))

  # bias correction
  if (bc !="none") {
    Q <- genDesign(x, method, q, J, bsmooth, bknot, rep(0, d))
    invG.q <- qrXXinv(Q)
    if (bc == "bc3") {
      basis.q  <- genBias(eval, method, m, q, knot, bknot, bsmooth, deriv)
      bias.all <- genBias(x, method, m, q, knot, bknot, bsmooth, rep(0, d))
      hii.q <- hii.p + rowSums((bias.all %*% invG.q) * Q)
      if ((method == "bs" | method=="pp") & proj == TRUE) {
        proj.bias <- basis.p %*% invG.p %*% crossprod(P, bias.all)
        basis.q   <- basis.q - proj.bias
        hii.q <- hii.q - rowSums((P %*% invG.p %*% crossprod(P, bias.all) %*% invG.q) * Q)
      }
    } else if (bc == "bc2") {
      basis.q <- genDesign(eval, method, q, J, bsmooth, bknot, deriv) -
                 genDesign(eval, method, m, J, smooth,  knot, deriv) %*% invG.p %*% crossprod(P, Q)
      hii.q <- rowSums((Q %*% invG.q) * Q)
      hii.q <- hii.p + hii.q - rowSums((P %*% invG.p %*% crossprod(P, Q) %*% invG.q) * Q)
    } else if (bc == "bc1") {
      basis.q <- genDesign(eval, method, q, J, bsmooth, bknot, deriv)
      hii.q <- rowSums((Q %*% invG.q) * Q)
    }
    hii.q[hii.q >= 0.9] <- 0.9; hii.q[hii.q < 0] <- 0
    d.q <- sum(hii.q)

    beta.q <- invG.q %*% crossprod(Q, y)
    if (bc == "bc1") tau.bc <- basis.q %*% beta.q
    else             tau.bc <- tau.cl + basis.q %*% beta.q

    if (bc != "bc3") {
      predict.q <- Q %*% beta.q
      if (bc == "bc2") predict.q <- predict.p + predict.q - P %*% invG.p %*% crossprod(P, Q) %*% beta.q
    } else {
      predict.q <- predict.p + genBias(x, method, m, q, knot, bknot, bsmooth, rep(0, d)) %*% beta.q
      if ((method == "bs" | method == "pp") & proj == TRUE) {
        predict.q <- predict.q - P %*% invG.p %*% crossprod(P, bias.all) %*% beta.q
      }
    }
    res.q <- lsprobust.res(y, predict.q, hii.q, vce, d.q)
    predict.p <- NULL; predict.q <- NULL

    if (bc == "bc1") {
      se.bc <- sqrt(rowSums((basis.q %*% invG.q %*% lsprobust.vce(Q, res.q) %*% invG.q) * basis.q))
    } else {
      se.bc <- sqrt(rowSums((basis.q %*% invG.q %*% lsprobust.vce(Q, res.q) %*% invG.q) * basis.q) +
               2*rowSums((basis.p %*% invG.p %*% lsprobust.cov(P, Q, res.q) %*% invG.q) * basis.q) +
               rowSums((basis.p %*% invG.p %*% lsprobust.vce(P, res.q) %*% invG.p) * basis.p))
    }
    basis.p <- NULL; basis.q <- NULL
  }


  #####################
  ##uniform inference##
  #####################
  uni.output   <- NA
  sup.quantile <- NA
  pl.num <- NA; wb.num.1 <- NA; wb.num.2 <- NA
  if (!is.null(uni.method)) {
    if (is.null(uni.grid)) {
      uni.grid <- list()
      for (j in 1:d) {
        uni.grid[[j]] <- seq(x.min[j], x.max[j], length.out = uni.ngrid+2)[2:(uni.ngrid+1)]
      }
      uni.grid <- as.matrix(expand.grid(uni.grid))
    }

    if (bc == "bc3") {
      design.p <- genDesign(uni.grid, method, m, J, smooth, knot, deriv)
      design.q <- genBias(uni.grid, method, m, q, knot, bknot, bsmooth, deriv)
      if ((method == "bs" | method == "pp") & proj == TRUE) {
        design.q      <- design.q - design.p %*% invG.p %*% crossprod(P, bias.all)
      }
      denom <- sqrt(rowSums((design.q %*% invG.q %*% lsprobust.vce(Q, res.q) %*% invG.q) * design.q) +
                    2*rowSums((design.p %*% invG.p %*% lsprobust.cov(P, Q, res.q) %*% invG.q) * design.q) +
                    rowSums((design.p %*% invG.p %*% lsprobust.vce(P, res.q) %*% invG.p) * design.p))
    } else if (bc == "bc2") {
      design.p <- genDesign(uni.grid, method, m, J, smooth, knot, deriv)
      design.q <- genDesign(uni.grid, method, q, J, bsmooth, bknot, deriv) -
                  genDesign(uni.grid, method, m, J, smooth, knot, deriv) %*% invG.p %*% crossprod(P, Q)
      denom    <- sqrt(rowSums((design.q %*% invG.q %*% lsprobust.vce(Q, res.q) %*% invG.q) * design.q) +
                       2*rowSums((design.p %*% invG.p %*% lsprobust.cov(P, Q, res.q) %*% invG.q) * design.q) +
                       rowSums((design.p %*% invG.p %*% lsprobust.vce(P, res.q) %*% invG.p) * design.p))
    } else if (bc == "bc1") {
      design.q <- genDesign(uni.grid, method, q, J, bsmooth, bknot, deriv)
      denom    <- sqrt(rowSums((design.q %*% invG.q %*% lsprobust.vce(Q, res.q) %*% invG.q) * design.q))
    } else {
      design.p <- genDesign(uni.grid, method, m, J, smooth, knot, deriv)
      denom    <- sqrt(rowSums((design.p %*% invG.p %*% lsprobust.vce(P, res.p) %*% invG.p) * design.p))
    }

    if (uni.method == "pl") {
      if (bc == "bc1") {
        Sigma      <- lsprobust.vce(Q, res.q)
        Sigma.root <- lssqrtm(Sigma)
        pl.num     <- design.q %*% invG.q %*% Sigma.root
      } else if (bc == "none") {
        Sigma      <- lsprobust.vce(P, res.p)
        Sigma.root <- lssqrtm(Sigma)
        pl.num     <- design.p %*% invG.p %*% Sigma.root
      } else {
        Sigma      <- lsprobust.vce(cbind(P, Q), res.q)
        Sigma.root <- lssqrtm(Sigma)
        pl.num     <- cbind(design.p %*% invG.p, design.q %*% invG.q) %*% Sigma.root
      }
    } else {
      if (bc == "bc1") {
        wb.num.1 <- design.q %*% invG.q
        wb.num.2 <- Q
      } else if (bc == "none") {
        wb.num.1 <- design.p %*% invG.p
        wb.num.2 <- P
      } else {
        wb.num.1 <- cbind(design.p %*% invG.p, design.q %*% invG.q)
        wb.num.2 <- Q; Q <- NULL
        wb.num.2 <- cbind(P, wb.num.2); P <- NULL
      }
    }
  }

  if (uni.out ==  TRUE) {
    if (bc == "none") {
      uni.output <- list(t.num.pl=pl.num, t.num.wb1=wb.num.1, t.num.wb2=wb.num.2,
                         t.denom=denom, res=res.p)
    } else {
      uni.output <- list(t.num.pl=pl.num, t.num.wb1=wb.num.1, t.num.wb2=wb.num.2,
                         t.denom=denom, res=res.q)
    }
  }

  if (band == TRUE) {
    if (uni.method == "pl") {
      sup.quantile <- lsprobust.sup.pl(pl.num, denom, ncol(Sigma.root), B, level)
    } else {
      if (bc == "none") {
        sup.quantile <- lsprobust.sup.wb(wb.num.1, wb.num.2, denom, res.p, N, B, level)
      } else {
        sup.quantile <- lsprobust.sup.wb(wb.num.1, wb.num.2, denom, res.q, N, B, level)
      }
    }
  }

  ################
  ####output######
  ################
  namlist  <- paste("X", 1:d, sep = "")
  if (bc == "none") {
    Estimate <- cbind(eval, eN, tau.cl, se.cl)
    colnames(Estimate) <- c(namlist, "N", "tau.cl", "se.cl")
  } else {
    Estimate <- cbind(eval, eN, tau.cl, tau.bc, se.cl, se.bc)
    colnames(Estimate) <- c(namlist, "N", "tau.cl", "tau.bc", "se.cl", "se.rb")
  }

  out <- list(Estimate=Estimate, k.num=k.num, sup.cval=sup.quantile, uni.output=uni.output,
              knot=knot, bknot=bknot, opt=list(m=m+1, m.bc=q+1, deriv=deriv.str,
              method=method.type, ktype=knot.type, n=N, d=d, neval=neval, J=J,
              k=k.m.str, k.bc=k.b.str, kselect=kselect.type, bc=bc,
              smooth=smooth.p, bsmooth=smooth.q, uni.method=uni.method.str))
  out$call <- match.call()
  class(out) <- "lsprobust"
  return(out)
}

#'@describeIn lsprobust \code{print} method for class "\code{lsprobust}"
#'@param ... further arguments
#'@export
print.lsprobust <- function(x, ...){
  cat("Call: lsprobust\n\n")

  cat(paste("Sample size (n)                            =    ", x$opt$n,         "\n", sep=""))
  cat(paste("Num. covariates (d)                        =    ", x$opt$d,         "\n", sep=""))
  cat(paste("Basis function (method)                    =    ", x$opt$method,    "\n", sep=""))
  cat(paste("Order of basis point estimation (m)        =    ", x$opt$m,         "\n", sep=""))
  cat(paste("Order of derivative (deriv)                =    ", x$opt$deriv,     "\n", sep=""))
  cat(paste("Order of basis bias correction (m.bc)      =    ", x$opt$m.bc,      "\n", sep=""))
  cat(paste("Smoothness point estimation (smooth)       =    ", x$opt$smooth,    "\n", sep=""))
  cat(paste("Smoothness bias correction (bsmooth)       =    ", x$opt$bsmooth,   "\n", sep=""))
  cat(paste("Knot placement (ktype)                     =    ", x$opt$ktype,     "\n", sep=""))
  cat(paste("Knots method (kselect)                     =    ", x$opt$kselect,   "\n", sep=""))
  cat(paste("Uniform inference method (uni.method)      =    ", x$opt$uni.method,"\n", sep=""))
  cat(paste("Num. knots point estimation (nknot)        =    ", x$opt$k,         "\n", sep=""))
  cat(paste("Num. knots bias correction (bnknot)        =    ", x$opt$k.bc,      "\n", sep=""))
  cat("\n")

  # cat("Use summary(...) to show estimates.\n")

}

#'@describeIn lsprobust \code{summary} method for class "\code{lsprobust}"
#'@param object class \code{lsprobust} objects.
#'@export
summary.lsprobust <- function(object, ...) {
  x    <- object
  args <- list(...)
  if (is.null(args[['alpha']])) { alpha <- 0.05 } else { alpha <- args[['alpha']] }
  if (is.null(args[['sep']]))   { sep <- 5 } else { sep <- args[['sep']] }

  cat("Call: lprobust\n\n")

  cat(paste("Sample size (n)                             =    ", x$opt$n,          "\n", sep=""))
  cat(paste("Num. covariates (d)                         =    ", x$opt$d,          "\n", sep=""))
  cat(paste("Basis function (method)                     =    ", x$opt$method,     "\n", sep=""))
  cat(paste("Order of basis point estimation (m)         =    ", x$opt$m,          "\n", sep=""))
  cat(paste("Order of derivative (deriv)                 =    ", x$opt$deriv,      "\n", sep=""))
  cat(paste("Order of basis bias correction (m.bc)       =    ", x$opt$m.bc,       "\n", sep=""))
  cat(paste("Smoothness point estimation (smooth)        =    ", x$opt$smooth,     "\n", sep=""))
  cat(paste("Smoothness bias correction (bsmooth)        =    ", x$opt$bsmooth,    "\n", sep=""))
  cat(paste("Knot placement (ktype)                      =    ", x$opt$ktype,      "\n", sep=""))
  cat(paste("Knots method (kselect)                      =    ", x$opt$kselect,    "\n", sep=""))
  cat(paste("Uniform inference method (uni.method)       =    ", x$opt$uni.method, "\n", sep=""))
  cat(paste("Num. knots point estimation (nknot)         =    ", x$opt$k,          "\n", sep=""))
  cat(paste("Num. knots bias correction (bnknot)         =    ", x$opt$k.bc,       "\n", sep=""))
  cat("\n")

  ### compute CI
  z    <- qnorm(1 - alpha / 2)
  d <- x$opt$d; bc <- x$opt$bc

  if (bc == "none") {
    CI_l <- x$Estimate[, "tau.cl"] - x$Estimate[, "se.cl"] * z
    CI_r <- x$Estimate[, "tau.cl"] + x$Estimate[, "se.cl"] * z
  } else {
    CI_l <- x$Estimate[, "tau.bc"] - x$Estimate[, "se.rb"] * z
    CI_r <- x$Estimate[, "tau.bc"] + x$Estimate[, "se.rb"] * z
  }

  ### print output
  cat(paste(rep("=", 4+8*d + 8 + 10 + 10 + 25), collapse="")); cat("\n")

  cat(format(" "            , width= 4))
  cat(format("Eval"         , width= 8*d, justify="centre"))
  cat(format(" "            , width= 8))
  cat(format("Point"        , width= 10,  justify="right"))
  cat(format("Std."         , width= 10,  justify="right"))
  if (bc == "none") {
    cat(format("Classical"  , width=25,   justify="centre"))
  } else {
    cat(format("Robust B.C.", width=25,   justify="centre"))
  }
  cat("\n")

  cat(format(" "        , width=4))
  cat(format("X1"       , width=8, justify="centre"))
  if (d > 1) {
    for (j in 2:d) {
      cat(format(paste("X", j, sep = ""), width=8, justify="centre"))
    }
  }
  if (x$opt$method  == "Piecewise Polynomial") {
    cat(format("Eff.n"  , width=8 ,  justify="right"))
  } else {
    cat(format("n"      , width=8 ,  justify="right"))
  }
  cat(format("Est."     , width=10,  justify="right"))
  cat(format("Error"    , width=10,  justify="right"))
  cat(format(paste("[ ", floor((1-alpha)*100), "%", " C.I. ]", sep=""), width=25, justify="centre"))
  cat("\n")

  cat(paste(rep("=", 4+8*d + 8 + 10 + 10 + 25), collapse="")); cat("\n")

  for (j in 1:nrow(x$Estimate)) {
    cat(format(toString(j), width=4))
    for (i in 1:d) {
      cat(format(sprintf("%3.3f", x$Estimate[j, paste("X", i, sep="")]), width=8, justify="right"))
    }
    cat(format(sprintf("%3.0f", x$Estimate[j, "N"])  , width=8 , justify="right"))
    cat(format(sprintf("%3.3f", x$Estimate[j, "tau.cl"]) , width=10, justify="right"))
    cat(format(paste(sprintf("%3.3f", x$Estimate[j, "se.cl"]), sep=""), width=10, justify="right"))
    cat(format(paste("[", sprintf("%3.3f", CI_l[j]), " , ", sep="")  , width=14, justify="right"))
    cat(format(paste(sprintf("%3.3f", CI_r[j]), "]", sep=""), width=11, justify="left"))
    cat("\n")
    if (is.numeric(sep)) if (sep > 0) if (j %% sep == 0) {
      cat(paste(rep("-", 4+8*d + 8 + 10 + 10 + 25), collapse="")); cat("\n")
    }
  }

  cat(paste(rep("=", 4+8*d + 8 + 10 + 10 + 25), collapse="")); cat("\n")
}
