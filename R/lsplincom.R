#'Linear Combination of Estimators for \pkg{lspartition} Package
#'
#'@description \code{lsplincom} implements user-specified linear combinations across different data sub-groups for regression functions estimation, and computes corresponding (pointwise and uniform) robust bias-corrected inference measures. Estimation and inference is implemented using the \pkg{lspartition} package.
#'             See \href{https://sites.google.com/site/nppackages/lspartition/Cattaneo-Farrell_2013_JoE.pdf?attredirects=0}{Cattaneo and Farrell (2013)} and \href{https://arxiv.org/abs/1804.04916}{Cattaneo, Farrell and Feng (2019a)} for complete details.
#'
#'             A detailed introduction to this command is given in \href{https://arxiv.org/abs/1906.00202}{Cattaneo, Farrell and Feng (2019b)}.
#'
#'             For more details, and related Stata and R packages useful for empirical analysis,
#'             visit \url{https://sites.google.com/site/nppackages/}.
#'
#'@param y Outcome variable.
#'@param x Independent variable. A matrix or data frame.
#'@param G Group indicator. It may take on multiple discrete values.
#'@param R A numeric vector giving the linear combination of interest. Each element is the coefficient
#'         of the conditional mean estimator of one group, and they are ordered ascendingly along the value
#'         of \code{G}.
#'@param eval Evaluation points. A matrix or data frame.
#'@param neval Number of quantile-spaced evaluating points.
#'@param level Confidence level used for confidence intervals; default is \code{level=95}.
#'@param band If \code{TRUE}, the critical value for constructing confidence band is calculated. Default
#'            is \code{band=FALSE}.
#'@param cb.method Method used to calculate the critical value for confidence bands.
#'                 Options are \code{"pl"} for a simulation-based plug-in procedure, and
#'                 \code{"wb"} for a wild bootstrap procedure. If \code{band=TRUE} with
#'                 \code{cb.method} unspecified, default is \code{cb.method="pl"}.
#'@param cb.grid  A matrix containing all grid points used to construct confidence bands. Each row
#'                 correponds to the coordinates of one grid point.
#'@param cb.ngrid A numeric vector of the same length as \code{ncol(x)}. Each element corresponds to
#'                 the number of grid points for each dimension used to implement uniform inference.
#'                 Default is \code{uni.ngrid=50}.
#'@param B    Number of simulated samples used to obtain the critical value for confidence bands.
#'            Default is \code{B=1000}.
#'@param subset Optional rule specifying a subset of observations to be used.
#'@param knot A list of numeric vectors giving the knot positions (including boundary knots) for each dimension
#'            which are used in the main regression. The length of the list is equal to \code{ncol(x)}.
#'            If not specified, it uses the number of knots either specified by users
#'            or computed by the companion command \code{lspkselect} to generate the
#'            corresponding knots according to the rule specified by \code{ktype}. See help for \code{\link{lsprobust}}.
#'@param ... Arguments to be passed to the function. See \code{lsprobust}.
#'@return \item{\code{Estimate}}{ A matrix containing eval (grid points), N (effective sample sizes),
#'                             tau.cl (point estimates with a basis of order \code{m}), tau.bc (bias corrected point
#'                             estimates with a basis of order \code{m.bc}), se.cl (standard error corresponding
#'                             to tau.cl), and se.rb (robust standard error).}
#'        \item{\code{sup.cval}}{ Critical value for constructing confidence bands.}
#'        \item{\code{opt}}{ A list containing options passed to the function.}
#'@author
#' Matias D. Cattaneo, Princeton University, Princeton, NJ. \email{cattaneo@princeton.edu}.
#'
#' Max H. Farrell, University of Chicago, Chicago, IL. \email{max.farrell@chicagobooth.edu}.
#'
#' Yingjie Feng (maintainer), Princeton University, Princeton, NJ. \email{yjfeng@umich.edu}.
#'
#'@references
#'
#' Cattaneo, M. D., M. H. Farrell, and Y. Feng (2019a): \href{https://arxiv.org/abs/1804.04916}{Large Sample Properties of Partitioning-Based Series Estimators}. Annals of Statistics, forthcoming. arXiv:1804.04916.
#'
#' Cattaneo, M. D., M. H. Farrell, and Y. Feng (2019b): \href{https://arxiv.org/abs/1906.00202}{lspartition: Partitioning-Based Least Squares Regression}. Working paper.
#'
#'@seealso \code{\link{lsprobust}}, \code{\link{lspkselect}}, \code{\link{lsprobust.plot}},
#'
#'@examples
#'x   <- runif(500)
#'y   <- sin(4*x)+rnorm(500)
#'z   <- c(rep(0, 250), rep(1, 250))
#'est <- lsplincom(y, x, z, c(-1, 1))
#'summary(est)
#'
#'@export


# version 0.3 May2019
lsplincom <- function(y, x, G, R, eval=NULL, neval=NULL, level=95, band=FALSE, cb.method=NULL,
                      cb.grid=NULL, cb.ngrid=50, B=1000, subset=NULL, knot=NULL, ...) {

  x <- as.matrix(x); R <- as.matrix(R)
  if (is.data.frame(y)) y <- y[,1]
  if (is.data.frame(G)) G <- G[,1]
  #substract subset
  if (!is.null(subset)) {
    x <- x[subset, , drop = F]
    y <- y[subset]
    G <- G[subset]
  }

  na.ok <- complete.cases(x) & complete.cases(y) & complete.cases(G)

  x <- x[na.ok, , drop = F]
  y <- y[na.ok]
  G <- G[na.ok]

  N <- nrow(x); d <- ncol(x); G.val <- sort(unique(G)); g <- length(G.val)
  x.max <- sapply(G.val, function(j) colMaxs(as.matrix(x[G==j,])))
  x.min <- sapply(G.val, function(j) colMins(as.matrix(x[G==j,])))
  x.max <- rowMins(matrix(x.max, d, g))
  x.min <- rowMaxs(matrix(x.min, d, g))

  if (d == 1 & is.vector(knot, mode="numeric")) {
    knot <- list(knot)
  }

  if (is.null(eval)) {
    if (!is.null(knot)) {
      if (is.null(neval)) {
        eval <- lapply(knot, function(z) cumsum(c(z[1], rep(diff(z)/11, each=11)))[-seq(1, by=11, length.out=length(z))])
        eval <- as.matrix(expand.grid(eval))
      } else {
        qseq <- seq(0, 1, 1/(neval+1))
        eval <- apply(x, 2, function(z) quantile(z, qseq[2:(length(qseq)-1)], names=F))
        if (neval == 1) eval <- t(eval)
      }
    } else {
      ind <- rowSums((sweep(x, 2, x.min)>=0) + (sweep(x, 2, x.max)<=0)) == 2*d
      if (is.null(neval)) {
        qseq <- seq(0, 1, 1/(20+1))
        eval <- apply(x[ind,, drop=F], 2, function(z) quantile(z, qseq[2:(length(qseq)-1)]))
      } else {
        qseq <- seq(0, 1, 1/(neval+1))
        eval <- apply(x[ind,, drop=F], 2, function(z) quantile(z, qseq[2:(length(qseq)-1)]))
        if (neval == 1) eval <- t(eval)
      }
    }
  }

  eval  <- as.matrix(eval)
  neval <- nrow(eval)

  ng <- c()
  tau.cl <- matrix(NA, nrow=neval, ncol=g)
  tau.bc <- matrix(NA, nrow=neval, ncol=g)
  se.cl  <- matrix(NA, nrow=neval, ncol=g)
  se.bc  <- matrix(NA, nrow=neval, ncol=g)
  pl.num.list <- list(); wb1.num.list <- list(); wb2.num.list <- list()
  denom.list <- list(); res.list <- list()
  uni.out <- FALSE; sup.quantile <- NA

  if (band == TRUE) {
    uni.out <- TRUE
    if (is.null(cb.method)) cb.method <- "pl"
    if (is.null(cb.grid)) {
        cb.grid <- list()
        for (j in 1:d) {
          cb.grid[[j]] <- seq(x.min[j], x.max[j], length.out = cb.ngrid+2)[2:(cb.ngrid+1)]
        }
        cb.grid <- as.matrix(expand.grid(cb.grid))
    }
  }

  for (j in 1:g) {
    y.j        <- y[G==G.val[j]]; x.j <- x[G==G.val[j], , drop = F]
    fit        <- lsprobust(y.j, x.j, eval, neval, level=level, uni.out=uni.out, knot=knot,
                            uni.method=cb.method, uni.grid=cb.grid, band=FALSE, ...)
    ng[j]      <- fit$opt$n
    tau.cl[,j] <- fit$Estimate[, "tau.cl"]
    tau.bc[,j] <- fit$Estimate[, "tau.bc"]
    se.cl[,j]  <- fit$Estimate[, "se.cl"]
    se.bc[,j]  <- fit$Estimate[, "se.rb"]
    if (band == TRUE) {
      if (cb.method == "pl") {
        pl.num.list[[j]] <- fit$uni.output$t.num.pl
        denom.list[[j]]  <- fit$uni.output$t.denom
      } else {
        wb1.num.list[[j]] <- fit$uni.output$t.num.wb1
        wb2.num.list[[j]] <- fit$uni.output$t.num.wb2
        denom.list[[j]]   <- fit$uni.output$t.denom
        res.list[[j]]     <- fit$uni.output$res
      }
    }
  }

  m    <- fit$opt$m
  m.bc <- fit$opt$m.bc
  deriv.str <- fit$opt$deriv
  method.type <- fit$opt$method
  knot.type <- fit$opt$ktype
  J <- fit$opt$J
  kselect.type <- fit$opt$kselect
  smooth.p <- fit$opt$smooth
  smooth.q <- fit$opt$bsmooth

  eval <- fit$Estimate[, 1:d]
  tau.cl <- tau.cl %*% R; tau.bc <- tau.bc %*% R
  se.cl <- sqrt(rowSums(sweep(se.cl^2, 2, R^2, FUN = "*")))
  se.bc <- sqrt(rowSums(sweep(se.bc^2, 2, R^2, FUN = "*")))

  Estimate <- cbind(eval, tau.cl, tau.bc, se.cl, se.bc)
  namlist  <- paste("X", 1:d, sep = "")
  colnames(Estimate) <- c(namlist, "tau.cl", "tau.bc", "se.cl", "se.rb")

  cb.method.str <- cb.method
  if (is.null(cb.method.str)) {
    cb.method.str <- NA
  } else if (cb.method.str == "pl") {
    cb.method.str <- "Plug-in"
  } else if (cb.method.str == "wb") {
    cb.method.str <- "Bootstrap"
  }

  if (band == TRUE) {
    temp.sup <- rep(NA, B)
    if (cb.method == "pl") {
      cols <- sapply(pl.num.list, FUN = ncol)
      rows <- nrow(pl.num.list[[1]])
      for (i in 1:B) {
        zeta    <- matrix(rnorm(sum(cols), 0, 1), ncol = 1)
        index.r <- cumsum(cols); index.l <- c(0, index.r[-length(index.r)]) + 1
        num     <- sapply(1:g, function(z) pl.num.list[[z]] %*% zeta[(index.l[z]: index.r[z]),] * R[z,])
        denom   <- sapply(1:g, function(z) denom.list[[z]]^2 * R[z,]^2)
        num     <- rowSums(matrix(unlist(num), nrow=rows))
        denom   <- sqrt(rowSums(matrix(unlist(denom), nrow=rows)))
        temp.sup[i] <- max(abs(num / denom))
      }
    } else {
      cols <- sapply(wb2.num.list, FUN = nrow)
      rows <- nrow(wb1.num.list[[1]])
      for (i in 1:B) {
        w <- rbinom(sum(cols), 1, 0.5) * 2 - 1
        index.r <- cumsum(cols); index.l <- c(0, index.r[-length(index.r)]) + 1
        num     <- sapply(1:g, function(z) wb1.num.list[[z]] %*%
                               crossprod(wb2.num.list[[z]], (res.list[[z]] *
                               w[index.l[z]: index.r[z]])) * R[z,])
        denom   <- sapply(1:g, function(z) denom.list[[z]]^2 * R[z, ]^2)
        num     <- rowSums(matrix(unlist(num), nrow=rows))
        denom   <- sqrt(rowSums(matrix(unlist(denom), nrow=rows)))
        temp.sup[i] <- max(abs(num / denom))
      }
    }
    sup.quantile <- quantile(temp.sup, level/100)
  }

  out <- list(Estimate=Estimate, sup.cval=sup.quantile,
              opt=list(m=m, m.bc=m.bc, deriv=deriv.str, method=method.type, ktype=knot.type,
                       n=N, d=d, g=g, ng=ng, neval=neval, J=J, kselect=kselect.type, smooth=smooth.p,
                       bsmooth=smooth.q, cb.method=cb.method.str))
  out$call <- match.call()
  class(out) <- "lsplincom"
  return(out)
}

#'@describeIn lsplincom \code{print} method for class "\code{lsplincom}".
#'@export
print.lsplincom <- function(x, ...){
  cat("Call: lsplincom\n\n")

  cat(paste("Sample size (n)                            =    ", x$opt$n,         "\n", sep=""))
  cat(paste("Num. covariates (d)                        =    ", x$opt$d,         "\n", sep=""))
  cat(paste("Num. groups (G)                            =    ", x$opt$g,         "\n", sep=""))
  cat(paste("Basis function (method)                    =    ", x$opt$method,    "\n", sep=""))
  cat(paste("Order of basis point estimation (m)        =    ", x$opt$m,         "\n", sep=""))
  cat(paste("Order of derivative (deriv)                =    ", x$opt$deriv,     "\n", sep=""))
  cat(paste("Order of basis bias correction (m.bc)      =    ", x$opt$m.bc,      "\n", sep=""))
  cat(paste("Smoothness point estimation (smooth)       =    ", x$opt$smooth,    "\n", sep=""))
  cat(paste("Smoothness bias correction (bsmooth)       =    ", x$opt$bsmooth,   "\n", sep=""))
  cat(paste("Knot placement (ktype)                     =    ", x$opt$ktype,     "\n", sep=""))
  cat(paste("Knots method (kselect)                     =    ", x$opt$kselect,   "\n", sep=""))
  cat(paste("Confidence band method (cb.method)         =    ", x$opt$cb.method, "\n", sep=""))
  #cat(paste("Num. knots for the main regression        =    ", x$k.num[1,1],    "\n", sep=""))
  #cat(paste("Num. knots for bias correction            =    ", x$k.num[1,2],    "\n", sep=""))
  cat("\n")

  # cat("Use summary(...) to show estimates.\n")
}

#'@describeIn lsplincom \code{summary} method for class "\code{lsplincom}"
#'@param object class \code{lsplincom} objects.
#'@export
summary.lsplincom <- function(object,...) {
  x    <- object
  args <- list(...)
  if (is.null(args[['alpha']])) { alpha <- 0.05 } else { alpha <- args[['alpha']] }
  if (is.null(args[['sep']]))   { sep <- 5 } else { sep <- args[['sep']] }

  cat("Call: lprobust\n\n")

  cat(paste("Sample size (n)                            =    ", x$opt$n,          "\n", sep=""))
  cat(paste("Num. covariates (d)                        =    ", x$opt$d,          "\n", sep=""))
  cat(paste("Num. groups (G)                            =    ", x$opt$g,          "\n", sep=""))
  cat(paste("Basis function (method)                    =    ", x$opt$method,     "\n", sep=""))
  cat(paste("Order of basis point estimation (m)        =    ", x$opt$m,          "\n", sep=""))
  cat(paste("Order of derivative (deriv)                =    ", x$opt$deriv,      "\n", sep=""))
  cat(paste("Order of basis bias correction (m.bc)      =    ", x$opt$m.bc,       "\n", sep=""))
  cat(paste("Smoothness point estimation (smooth)       =    ", x$opt$smooth,     "\n", sep=""))
  cat(paste("Smoothness bias correction (bsmooth)       =    ", x$opt$bsmooth,    "\n", sep=""))
  cat(paste("Knot placement (ktype)                     =    ", x$opt$ktype,      "\n", sep=""))
  cat(paste("Knots method (kselect)                     =    ", x$opt$kselect,    "\n", sep=""))
  cat(paste("Confidence band method (cb.method)         =    ", x$opt$cb.method,  "\n", sep=""))
  #cat(paste("Num of knots for the main regression      =    ", x$k.num[1,1],     "\n", sep=""))
  #cat(paste("Num of knots for bias correction          =    ", x$k.num[1,2],     "\n", sep=""))
  cat("\n")

  ### compute CI
  z    <- qnorm(1 - alpha / 2)
  CI_l <- x$Estimate[, "tau.bc"] - x$Estimate[, "se.rb"] * z;
  CI_r <- x$Estimate[, "tau.bc"] + x$Estimate[, "se.rb"] * z;

  d <- x$opt$d
  ### print output
  cat(paste(rep("=", 4+8*d + 10 + 10 + 25), collapse="")); cat("\n")

  cat(format(" ", width= 4))
  cat(format("Eval", width= 8*d, justify="centre"))
  cat(format("Point", width= 10, justify="right"))
  cat(format("Std." , width= 10, justify="right"))
  cat(format("Robust B.C.", width=25, justify="centre"))
  cat("\n")

  cat(format(" ", width=4))
  cat(format("X1"              , width=8, justify="centre"))
  if (d > 1) {
    for (j in 2:d) {
      cat(format(paste("X", j, sep = ""), width=8, justify="centre"))
    }
  }
  cat(format("Est."            , width=10,  justify="right"))
  cat(format("Error"           , width=10,  justify="right"))
  cat(format(paste("[ ", floor((1-alpha)*100), "%", " C.I. ]", sep=""), width=25, justify="centre"))
  cat("\n")

  cat(paste(rep("=", 4+8*d + 10 + 10 + 25), collapse="")); cat("\n")

  for (j in 1:nrow(x$Estimate)) {
    cat(format(toString(j), width=4))
    for (i in 1:d) {
      cat(format(sprintf("%3.3f", x$Estimate[j, paste("X", i, sep="")]), width=8, justify="right"))
    }
    cat(format(sprintf("%3.3f", x$Estimate[j, "tau.cl"]) , width=10, justify="right"))
    cat(format(paste(sprintf("%3.3f", x$Estimate[j, "se.cl"]), sep=""), width=10, justify="right"))
    cat(format(paste("[", sprintf("%3.3f", CI_l[j]), " , ", sep="")  , width=14, justify="right"))
    cat(format(paste(sprintf("%3.3f", CI_r[j]), "]", sep=""), width=11, justify="left"))
    cat("\n")
    if (is.numeric(sep)) if (sep > 0) if (j %% sep == 0) {
      cat(paste(rep("-", 4+8*d + 10 + 10 + 25), collapse="")); cat("\n")
    }
  }

  cat(paste(rep("=", 4+8*d + 10 + 10 + 25), collapse="")); cat("\n")
}
