#'Graphic Presentation of Results for \pkg{lspartition} Package
#'
#'@description \code{lsprobust.plot} plots estimated regression functions and confidence regions using the \pkg{lspartition} package.
#'             See \href{https://sites.google.com/site/nppackages/lspartition/Cattaneo-Farrell_2013_JoE.pdf?attredirects=0}{Cattaneo and Farrell (2013)} and \href{https://arxiv.org/abs/1804.04916}{Cattaneo, Farrell and Feng (2019a)} for complete details.
#'
#'             Companion command: \code{\link{lsprobust}} for partitioning-based least squares regression
#'             estimation and inference; \code{\link{lsprobust.plot}} for plotting results; \code{\link{lsplincom}} for multiple sample estimation and inference.
#'
#'             A detailed introduction to this command is given in \href{https://arxiv.org/abs/1906.00202}{Cattaneo, Farrell and Feng (2019b)}.
#'
#'             For more details, and related Stata and R packages useful for empirical analysis,
#'             visit \url{https://sites.google.com/site/nppackages/}.
#'
#'@param ... Objects returned by \code{\link{lsprobust}}.
#'@param alpha Numeric scalar between 0 and 1, the significance level for plotting
#'             confidence regions. If more than one is provided, they will be applied
#'             to data series accordingly.
#'@param type String, one of \code{"line"} (default), \code{"points"}, \code{"binscatter"},
#'            \code{"none"} or \code{"both"}, how the point estimates are plotted. If more
#'            than one is provided, they will be applied to data series accordingly.
#'@param CS String, type of confidence sets. Options are \code{"ci"} for pointwise confidence
#'          intervals, \code{"cb"} for uniform confidence bands, and \code{"all"} for both.
#'@param CStype String, one of \code{"region"} (shaded region, default), \code{"line"}
#'              (dashed lines), \code{"ebar"} (error bars), \code{"all"} (all of the previous)
#'              or \code{"none"} (no confidence region), how the confidence region should
#'              be plotted. If more than one is provided, they will be applied to data series accordingly.
#'              If \code{CS = "all"}, pointwise confidence intervals are forced to be represented by error bars,
#'              and uniform bands are represented by both lines and regions.
#'@param title String, title of the plot.
#'@param xlabel Strings, labels for x-axis.
#'@param ylabel Strings, labels for y-axis.
#'@param lty Line type for point estimates, only effective if \code{type} is \code{"line"} or
#'           \code{"both"}. \code{1} for solid line, \code{2} for dashed line, \code{3}
#'           for dotted line. For other options, see the instructions for \code{\link{ggplot2}}
#'           or \code{\link{par}}. If more than one is provided, they will be applied to data
#'           series accordingly.
#'@param lwd Line width for point estimates, only effective if \code{type} is \code{"line"}
#'           or \code{"both"}. Should be strictly positive. For other options, see the
#'           instructions for \code{\link{ggplot2}} or \code{\link{par}}. If more than one
#'           is provided, they will be applied to data series accordingly.
#'@param lcol Line color for point estimates, only effective if \code{type} is \code{"line"} or
#'            \code{"both"}. \code{1} for black, \code{2} for red, \code{3} for green,
#'            \code{4} for blue. For other options, see the instructions for \code{\link{ggplot2}}
#'            or \code{\link{par}}. If more than one is provided, they will be applied to
#'            data series accordingly.
#'@param pty Scatter plot type for point estimates, only effective if \code{type} is
#'           \code{"points"} or \code{"both"}. For options, see the instructions for
#'           \code{\link{ggplot2}} or \code{\link{par}}. If more than one is provided,
#'           they will be applied to data series accordingly.
#'@param pwd Scatter plot size for point estimates, only effective if \code{type} is
#'           \code{"points"} or \code{"both"}. Should be strictly positive. If more than
#'           one is provided, they will be applied to data series accordingly.
#'@param pcol Scatter plot color for point estimates, only effective if \code{type} is
#'            \code{"points"} or \code{"both"}. \code{1} for black, \code{2} for red,
#'            \code{3} for green, \code{4} for blue. For other options, see the instructions
#'            for \code{\link{ggplot2}} or \code{\link{par}}. If more than one is provided,
#'            they will be applied to data series accordingly.
#'@param CSshade Numeric, opaqueness of the confidence region, should be between 0
#'               (transparent) and 1. Default is 0.2. If more than one is provided,
#'               they will be applied to data series accordingly.
#'@param CScol Color for confidence region. \code{1} for black, \code{2} for red, \code{3}
#'             for green, \code{4} for blue. For other options, see the instructions for
#'             \code{\link{ggplot2}} or \code{\link{par}}. If more than one is provided,
#'             they will be applied to data series accordingly.
#'@param legendTitle String, title of legend.
#'@param legendGroups String vector, group names used in legend.
#'
#'@details Companion command: \code{\link{lsprobust}} for partition-based least-squares regression
#'         estimation.
#'
#'@return A standard \code{\link{ggplot2}} object is returned, hence can be used for further
#'        customization.
#'
#'@author
#' Matias D. Cattaneo, Princeton University, Princeton, NJ. \email{cattaneo@princeton.edu}.
#'
#' Max H. Farrell, University of Chicago, Chicago, IL. \email{max.farrell@chicagobooth.edu}.
#'
#' Yingjie Feng (maintainer), Princeton University, Princeton, NJ. \email{yingjief@princeton.edu}.
#'
#'@seealso \code{\link{lsprobust}}, \code{\link{lspkselect}}, \code{\link{lsplincom}}, \code{\link{ggplot2}}.
#'
#'@references
#'
#' Cattaneo, M. D., M. H. Farrell, and Y. Feng (2019a): \href{https://arxiv.org/abs/1804.04916}{Large Sample Properties of Partitioning-Based Series Estimators}. Annals of Statistics, forthcoming. arXiv:1804.04916.
#'
#' Cattaneo, M. D., M. H. Farrell, and Y. Feng (2019b): \href{https://arxiv.org/abs/1906.00202}{lspartition: Partitioning-Based Least Squares Regression}. R Journal, forthcoming. arXiv:1906.00202.
#'
#'@examples
#'x   <- runif(500)
#'y   <- sin(4*x)+rnorm(500)
#'est <- lsprobust(y, x)
#'lsprobust.plot(est)
#'
#'@export

# version 0.4 Aug2019
lsprobust.plot <- function(..., alpha=NULL, type=NULL, CS="ci", CStype=NULL,
                          title="", xlabel="", ylabel="",
                          lty=NULL, lwd=NULL, lcol=NULL, pty=NULL, pwd=NULL, pcol=NULL,
                          CSshade=NULL, CScol=NULL, legendTitle=NULL, legendGroups=NULL) {

  ########################################
  # check how many series are passed in
  ########################################

  x <- list(...)
  nfig <- length(x)
  if (nfig == 0) stop("Nothing to plot.\n")

  ########################################
  # error handling
  ########################################
  # alpha
  if (length(alpha) == 0) {
    alpha <- rep(0.05, nfig)
  } else if (!all(alpha>0 & alpha<1)) {
    stop("Significance level incorrectly specified.\n")
  } else {
    alpha <- rep(alpha, length.out=nfig)
  }

  # plot type
  if (length(type) == 0) {
    type <- rep("line", nfig)
  } else {
    if (!all(type%in%c("line", "points", "both", "none", "binscatter"))) {
       stop("Plotting type incorrectly specified.\n")
    }
    type <- rep(type, length.out=nfig)
  }

  # CS type
  if (length(CStype) == 0) {
    CStype <- rep("region", nfig)
  } else {
    if (!all(CStype%in%c("region", "line", "ebar", "all", "none"))) {
      stop("Confidence set type incorrectly specified.\n")
    }
    if (CS == "cb" & (!all(CStype%in%c("region", "line", "all", "none")))) {
      stop("Confidence band type incorrectly specified.\n")
    }
    CStype <- rep(CStype, length.out=nfig)
  }
  CBtype <- rep("none", length.out=nfig)
  if (CS == "all") {
    CStype <- rep("ebar", length.out=nfig)
    CBtype <- rep("all", length.out=nfig)
    #if (nfig > 1)  CBtype <- rep("line", length.out=nfig)
  }

  # line style, line width, line color
  if (length(lty) == 0) {
    lty <- rep(1, nfig)
  } else {
    lty <- rep(lty, length.out=nfig)
  }
  if (length(lwd) == 0) {
    lwd <- rep(0.5, nfig)
  } else {
    lwd <- rep(lwd, length.out=nfig)
  }
  if (length(lcol) == 0) {
    lcol <- 1:nfig
  } else {
    lcol <- rep(lcol, length.out=nfig)
  }

  # point style, point width, point color
  if (length(pty) == 0) {
    pty <- rep(1, nfig)
  } else {
    pty <- rep(pty, length.out=nfig)
  }
  if (length(pwd) == 0) {
    pwd <- rep(1, nfig)
  } else {
    pwd <- rep(pwd, length.out=nfig)
  }
  if (length(pcol) == 0) {
    pcol <- lcol
  } else {
    pcol <- rep(pcol, length.out=nfig)
  }

  # CS shade, CS color
  if (length(CSshade) == 0) {
    CSshade <- rep(0.2, nfig)
  } else {
    CSshade <- rep(CSshade, length.out=nfig)
  }
  if (length(CScol) == 0) {
    CScol <- lcol
  } else {
    CScol <- rep(CScol, length.out=nfig)
  }

  # legend
  # Changes made by Xinwei
  if (length(legendTitle) == 0) {
    legendTitle <- ""
  } else {
    legendTitle <- legendTitle[1]
  }
  if (length(legendGroups) > 0) {
    legendGroups <- rep(legendGroups, length.out=nfig)
    legend_default <- FALSE
  } else {
    legend_default <- TRUE
  }

  ########################################
  # initializing plot
  ########################################
  temp_plot <- ggplot() + theme_bw() #+ theme(legend.position="none")

  CI_l <- CI_r <- CB_l <- CB_r <- tau.cl <- eval <- Sname <- mid <- tau <- NULL

  ########################################
  # looping over input models
  ########################################
  ### all colors
  col_all <- lty_all <- pty_all <- c()
  for (i in 1:nfig) {
    data_x <- data.frame(x[[i]]$Estimate[, c("X1", "tau.cl", "tau.bc", "se.cl", "se.rb")])
    colnames(data_x)[1] <- "eval"
    if (CS != "cb") {
      z_val <- qnorm(1 - alpha[i]/2)
      data_x$CI_l <- data_x$tau.bc - z_val * data_x$se.rb
      data_x$CI_r <- data_x$tau.bc + z_val * data_x$se.rb
    }
    if (CS != "ci") {
      c_val <- x[[i]]$sup.cval
      if (is.na(c_val)) stop("critical value not available.\n")
      data_x$CB_l <- data_x$tau.bc - c_val * data_x$se.rb
      data_x$CB_r <- data_x$tau.bc + c_val * data_x$se.rb
    }

    ## disconnect pieces and prepare binscatter data
    method <- x[[i]]$opt$method
    if (method == "Piecewise Polynomial") {
      knot <- x[[i]]$knot[[1]][2:(length(x[[i]]$knot[[1]])-1)]
      data_x <- data_x[order(data_x$eval),]
      left.index <- sapply(knot, function(z) sum(data_x$eval < z))
      left.index <- left.index[left.index > 0 & left.index < nrow(data_x)]
      if (length(left.index) > 0) {
        left.eval <- (data_x$eval[left.index] + data_x$eval[left.index+1]) / 2
        add.na <- cbind(left.eval, matrix(NA, length(left.eval), ncol(data_x)-1))
        colnames(add.na) <- names(data_x)
        data_x <- rbind(data_x, add.na)
      }
      if (type[i] == "binscatter") {
        if (x[[i]]$opt$m != 1) stop("basis order incorrectly specified.\n")
        knot <- x[[i]]$knot[[1]]
        if (length(knot)-1 != length(left.index)+1) stop("evaluation points insufficient.\n")
        bins <- data.frame(mid = (knot[-length(knot)] + knot[-1]) / 2,
                           tau = c(data_x$tau.cl[left.index],
                                   data_x$tau.cl[left.index[length(left.index)]+1]))
      }
    } else {
      if (type[i] == "binscatter") stop("type mismatch.\n")
    }

    # changes made by Xinwei
    if (legend_default) {
      data_x$Sname <- paste("Series", i, sep=" ")
      legendGroups <- c(legendGroups, data_x$Sname)
      if (type[i] == "binscatter") bins$Sname <- paste("Series", i, sep=" ")
    } else {
      data_x$Sname <- legendGroups[i]
      if (type[i] == "binscatter") bins$Sname <- legendGroups[i]
    }

    ########################################
    # add CS regions to the plot
    if (CStype[i]%in%c("region", "all")) {
      if (CS == "ci") temp_plot <- temp_plot + geom_ribbon(data=data_x, aes(x=eval, ymin=CI_l, ymax=CI_r),
                                                           alpha=CSshade[i], fill=CScol[i])
      if (CS == "cb") temp_plot <- temp_plot + geom_ribbon(data=data_x, aes(x=eval, ymin=CB_l, ymax=CB_r),
                                                           alpha=CSshade[i], fill=CScol[i])
    }
    if (CBtype[i] == "all")
      temp_plot <- temp_plot + geom_ribbon(data=data_x, aes(x=eval, ymin=CB_l, ymax=CB_r),
                                           alpha=CSshade[i], fill=CScol[i])

    ########################################
    # add CI lines to the plot
    if (CStype[i]%in%c("line", "all")) {
      if (CS == "ci") temp_plot <- temp_plot + geom_line(data=data_x, aes(x=eval, y=CI_l),
                                                         linetype=2, alpha=CSshade[i], col=CScol[i]) +
                                               geom_line(data=data_x, aes(x=eval, y=CI_r),
                                                         linetype=2, alpha=CSshade[i], col=CScol[i])
      if (CS == "cb") temp_plot <- temp_plot + geom_line(data=data_x, aes(x=eval, y=CB_l),
                                                         linetype=2, alpha=CSshade[i], col=CScol[i]) +
                                               geom_line(data=data_x, aes(x=eval, y=CB_r),
                                                         linetype=2, alpha=CSshade[i], col=CScol[i])
    }
    if (CBtype[i] == "all")
      temp_plot <- temp_plot + geom_line(data=data_x, aes(x=eval, y=CB_l),
                                         linetype=2, alpha=CSshade[i], col=CScol[i]) +
                               geom_line(data=data_x, aes(x=eval, y=CB_r),
                                         linetype=2, alpha=CSshade[i], col=CScol[i])
    ########################################
    # add error bars to the plot
    if (CStype[i]%in%c("ebar", "all"))
      temp_plot <- temp_plot + geom_errorbar(data=data_x[complete.cases(data_x[c("eval", "CI_l", "CI_r")]),],
                                             aes(x=eval, ymin=CI_l, ymax=CI_r),
                                             alpha=CSshade[i], col=CScol[i], linetype=1)

    ########################################
    # add lines to the plot
    if (type[i]%in%c("line", "both", "binscatter")) {
      temp_plot <- temp_plot + geom_line(data=data_x, aes(x=eval, y=tau.cl, colour=Sname,
                                                          linetype=Sname), size=lwd[i])
    }

    ########################################
    # add points to the plot
    if (type[i]%in%c("points", "both")) {
      temp_plot <- temp_plot + geom_point(data=data_x[complete.cases(data_x[c("eval", "tau.cl")]),],
                                          aes(x=eval, y=tau.cl, colour=Sname, shape=Sname), size=pwd[i])
    }

    if (type[i] == "binscatter") {
      temp_plot <- temp_plot + geom_point(data=bins, aes(x=mid, y=tau, colour=Sname, shape=Sname), size=pwd[i])
    }

    #######################################
    if (type[i] == "line") {
      col_all <- c(col_all, lcol[i])
      lty_all <- c(lty_all, lty[i])
      pty_all <- c(pty_all, NA)
    } else if (type[i] == "both" | type[i] == "binscatter") {
      col_all <- c(col_all, lcol[i])
      lty_all <- c(lty_all, lty[i])
      pty_all <- c(pty_all, pty[i])
    } else if (type[i] == "points") {
      col_all <- c(col_all, pcol[i])
      lty_all <- c(lty_all, NA)
      pty_all <- c(pty_all, pty[i])
    } else {
      col_all <- c(col_all, NA)
      lty_all <- c(lty_all, NA)
      pty_all <- c(pty_all, NA)
    }
  }

  ########################################
  # change color, line type and point shape back, and customize legend
  ########################################

  if (all(type%in%c("line", "points", "both", "binscatter"))) {
    index <- sort.int(legendGroups, index.return=TRUE)$ix
    temp_plot <- temp_plot + scale_color_manual(values = col_all[index]) +
                 scale_linetype_manual(values = lty_all[index]) +
                 scale_shape_manual(values = pty_all[index]) +
                 guides(colour=guide_legend(title=legendTitle)) +
                 guides(linetype=guide_legend(title=legendTitle)) +
                 guides(shape=guide_legend(title=legendTitle))
  }

  ########################################
  # add title, x and y labs
  ########################################
  temp_plot <- temp_plot + labs(x=xlabel, y=ylabel) + ggtitle(title)

  ########################################
  # return the plot
  ########################################
  return (temp_plot)
}
