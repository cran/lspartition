# Uniform knot list (including xmin and xmax as boundaries)
genKnot.u <- function(x.min, x.max, d, n) {
  knot <- list()
  for (j in 1:d) {
    knot[[j]] <- seq(x.min[j], x.max[j], length.out = n[j]+2)
  }
  return(knot)
}

# Quantile knot list (including xmin and xmax as boundaries)
genKnot.q <- function(x, d, n) {
  knot <- list()
  for (j in 1:d) {
    knot[[j]] <- quantile(x[,j], seq(0, 1, 1/(n[j]+1)), names = F)
  }
  return(knot)
}

# Generate index set
genIndex <- function(m, d, method) {
  if (method == "bs" | method=="wav") index.m <- diag(m+1, d)
  if (method == "pp")                 index.m <- t(xsimplex(d, m+1))
  return(index.m)
}

# Locate the closest knots to the left of evaluating points and interval length
locate.j <- function(eval.j, knot.j) {
  sorted <- sort.int(eval.j, index.return = T)
  tmp    <- knot.j[-length(knot.j)][knot.j[-length(knot.j)] <= max(eval.j)]
  count  <- sapply(tmp, function(z) which.max(sorted$x >= z))
  count  <- c(diff(count), length(eval.j)-count[length(count)]+1)
  tx     <- rep(tmp, count)
  hx     <- rep(diff(c(tmp, knot.j[length(tmp)+1])), count)
  tx     <- tx[order(sorted$ix)]
  hx     <- hx[order(sorted$ix)]
  return(list(tx=tx, hx=hx))
}

locate.x <- function(eval, d, knot) {
  raw <- sapply(1:d, function(j) unlist(locate.j(eval[,j], knot[[j]]), use.names = F))
  tx  <- raw[1:nrow(eval),, drop=F]
  hx  <- raw[-(1:nrow(eval)),, drop=F]
  return(list(tx=tx, hx=hx))
}

#Check n of eval
duplicates.drop <- function(data, names) {
  data <- dplyr::distinct(data, names, .keep_all = T)
  return(data)
}

gen.eN <- function(eval, x, d, knot) {
  neval <- nrow(eval)
  range <- locate.x(eval, d, knot)
  range <- cbind(1:neval, range$tx, range$tx+range$hx)
  range <- as.data.frame(cbind(1:neval, range[do.call(order, as.list(as.data.frame(range[,2:(d+1), drop=F]))),]))
  #range <- c(list(range), lapply(as.list(names(range[,3:(d+2)])), as.name), list(.keep_all=T))
  #uniq  <- do.call(distinct, range)
  uniq <- duplicates.drop(range, lapply(as.list(names(range[,3:(d+2)])), as.name))
  count <- c(diff(uniq[,1]),  neval-uniq[nrow(uniq), 1]+1)
  eN.1 <- apply(uniq[, 3:(d+2), drop=F], 1, function(z) rowSums(sweep(x, 2, z) >= 0))
  eN.2 <- apply(uniq[, -(1:(d+2)), drop=F], 1, function(z) rowSums(sweep(x, 2, z) <= 0))
  eN   <- colSums(eN.1 + eN.2 == 2*d)
  eN   <- rep(eN, count)[order(range[,2])]
  return(eN)
}

# Daubechies scaling function
dbDesign <- function(eval, knot, m, J) {
  m <- m+1
  x.max <- knot[length(knot)]; x.min <- knot[1]
  n <- length(eval); k <- 2^J; h <- (x.max - x.min)/ k;
  ext.knot <- c(rep(x.min, m), seq(from = x.min + h, by = h, length.out = 2^J-2*m),
                rep(x.max-(2*m-1)*h, m))
  ext.knot <- matrix(ext.knot, nrow = n, ncol = k, byrow = TRUE)
  precision <- 2^11
  index <- round((eval - ext.knot) / h * precision) + 2
  index[index < 2 | index > precision * (2*m-1) + 2] <- 1
  P <- matrix(NA, nrow = n, ncol = k)
  db <- filelist[[paste0("db", m, ".txt")]]
  for (i in 1:m) P[,i] <- db[,i][index[,i]]
  if (m+1 <= k-m) P[, (m+1):(k-m)] <- matrix(db[,m+1][index[, (m+1):(k-m)]], nrow = n, ncol = k - 2*m)
  for (i in (k-m+1):k) P[,i] <- db[,i-k+m+m+1][index[,i]]
  return(P)
}

# Marginal design for bs and wav
design.margin <- function(eval, method, m, J, smooth, knot, deriv) {
   if (method == "bs") {
      if (length(knot)-1>=2) {
        ext.knot <- c(rep(knot[1], m+1), rep(knot[2:(length(knot)-1)], each = m-smooth), rep(knot[length(knot)], m+1))
      } else {
        ext.knot <- c(rep(knot[1], m+1), rep(knot[length(knot)], m+1))
      }
      P <- splineDesign(knots = ext.knot, eval, ord = m+1, derivs = rep(deriv, length(eval)))
   } else if (method == "wav") {
      P <- dbDesign(eval, knot, m, J)
   }
   return(P)
}

# Calculate power series matrix
genPower <- function(x, d, k) {
  P <- 1
  for (j in 1:d) {
    P <- P * outer(x[,j], k[j,], "^")
  }
  return(P)
}

# Modified poly function
gpoly <- function(eval, m, deriv) {
  if (m >= 1) {
    if (nrow(eval) == 1) {
      P.level <- poly(rbind(eval, eval), degree = m, raw = T, simple = T)[1,, drop=F]
    } else {
      P.level <- poly(eval, degree = m, raw = T, simple = T)
    }
    if (!all(deriv == rep(0, ncol(eval)))) {
      d       <- ncol(eval)
      name    <- colnames(P.level)
      de      <- as.matrix(sapply(name, function(x) as.numeric(unlist(strsplit(x, "[.]")))))
      if (d == 1) de <- t(de)
      de.low  <- sweep(de, 1, deriv)
      index <- which(colSums(de.low < 0) > 0)
      P.level[,index]    <- 0
      if (length(index) != ncol(P.level)) {
        if (length(index) == 0) {
          index <- 1:ncol(P.level)
        } else {
          index <- (1:ncol(P.level))[-index]
        }
        power     <- de[,index, drop=F]
        power.low <- de.low[,index, drop=F]
        P.level[,index] <- genPower(eval, d, power.low)
        P.level[,index] <- sweep(P.level[,index, drop=F], 2,
                                 colProds(factorial(power)/factorial(power.low)),
                                 FUN = "*")
      }
      P.level <- cbind(rep(0, nrow(eval)), P.level)
    } else {
      P.level <- cbind(rep(1, nrow(eval)), P.level)
    }
  } else if (m == 0) {
    P.level <- as.matrix(rep(1, nrow(eval)))
  }
  return(P.level)
}

# Computation function used in locpoly
compute.locpoly <- function(eval, coord.tx, coord.hx, coord.ehx, d.loc, m, deriv) {
  P <- matrix(0, nrow(eval), d.loc)
  index <- scale(eval, center=coord.tx, scale=coord.ehx)
  index <- which(rowSums((index>=0) + (index<1)) == 2*ncol(eval))
  eval  <- scale(eval, center=coord.tx, scale=coord.hx)
  if (length(index) != 0) {
    P[index, ] <- gpoly(eval[index,, drop=F], m, deriv) / prod(coord.hx^deriv)
  }
  return(P)
}

# Partitioned polynomial basis
locpoly <- function(eval, m, knot, deriv) {
  d   <- ncol(eval)
  coord.tx  <- as.matrix(expand.grid(lapply(knot, function(s) s[-length(s)])))
  coord.hx  <- as.matrix(expand.grid(lapply(knot, function(s) diff(s))))
  coord.ehx <- as.matrix(expand.grid(lapply(knot, function(s)
                         diff(c(s[1:(length(s)-1)], s[length(s)]+0.1)))))
  d.loc <- choose(d + m, m)
  P <- sapply(seq_len(nrow(coord.tx)), function(i) compute.locpoly(eval, coord.tx[i,], coord.hx[i,],
                                                   coord.ehx[i,], d.loc, m, deriv))
  P <- matrix(P, nrow(eval), nrow(coord.tx) * d.loc)
  return(P)
}

# Generate tensor-product basis or partitioning basis
genDesign <- function(eval, method, m, J, smooth, knot, deriv) {
  if (method == "bs" | method == "wav") {
    S <- list()
    for (j in 1 : ncol(eval)) {
      S[[j]] <- design.margin(eval[,j], method, m, J[j], smooth, knot[[j]], deriv[j])
    }
    P <- tensor.prod.model.matrix(S)
  } else if (method == "pp") {
    P <- locpoly(eval, m, knot, deriv)
  }
  return(P)
}

# Shifted univariate legendre polynomial
legen.j <- function(degree, x.j) {
  P <- sapply(0:degree, function(k) choose(degree, k) * choose(degree+k, k) * (-x.j)^k)
  if (length(x.j) == 1) {
    P <- (-1)^degree * sum(P)
  } else {
    P <- (-1)^degree * rowSums(P)
  }
  return(P)
}

# Multivariate legendre polynomial
legpoly <- function(degree, x) {
  L <- sapply(seq_len(length(degree)), function(j) legen.j(degree[j], x[,j]))
  if (nrow(x) == 1) {
    L <- prod(L)
  } else {
    L <- rowProds(L)
  }
  return(L)
}

# Generate bias vector
genBias <- function(eval, method, m, q, knot, bknot, bsmooth, deriv) {
  d <- ncol(eval)
  pos <- locate.x(eval, d, knot)
  B <- 0
  if (method == "bs") {
    index.m <- diag(m+1, d)
    index <- which(rowSums(sweep(index.m, 2, deriv) >= 0) == d)
    if (length(index) != 0) {
      for (j in index) {
          order <- m + 1 - deriv[j]
          deriv.bias <- index.m[j,]
          B <- B + bernoulli(order, (eval[,j] - pos$tx[,j]) / pos$hx[,j]) *
               pos$hx[,j]^order / factorial(order) *
               genDesign(eval=eval, method=method, m=q, smooth=bsmooth,
                         knot=bknot, deriv=deriv.bias)
      }
    }
  } else if (method == "pp") {
    index.m <- t(xsimplex(d, m+1))
    index <- which(rowSums(sweep(index.m, 2, deriv) >= 0) == d)
    if (length(index) != 0) {
      for (j in index) {
        deriv.bias <- index.m[j,]; order <- deriv.bias - deriv
        B <- B + legpoly(order, (eval - pos$tx) / pos$hx) * colProds(t(pos$hx)^order) /
                 prod(factorial(order)) / prod(choose(2 * order, order)) *
                 genDesign(eval=eval, method=method, m=q, knot=bknot, deriv=deriv.bias)
      }
    }
  }
  return(B)
}

# Generate residuals
lsprobust.res <- function(y, m, hii, vce, d) {
  n <- length(y)
  res <- matrix(NA,n,1)
  if (vce == "hc0") w = 1
  else if (vce == "hc1") w = sqrt(n/(n-d))
  else if (vce == "hc2") w = sqrt(1/(1-hii))
  else                   w =      1/(1-hii)
  res[,1] = w * (y-m[,1])
  return(res)
}

lsprobust.vce <- function(X, res) {
  M = crossprod(c(res) * X)
  return(M)
}

lsprobust.cov <- function(X.p, X.q, res) {
  M <- crossprod(c(res) * X.p, c(res) * X.q)
  return(M)
}

qrXXinv <- function(x, ...) {
  inv <- try(chol2inv(chol(crossprod(x))), silent = T)
  if (inherits(inv, "try-error")) {
    warning('Gram is nearly singular')
    inv <- ginv(crossprod(x))
  }
  return(inv)
}

# Generate ROT with global polynomials (for rot, deriv option not provided,
# since rates are the same as levels)
lspkselect.imse.rot <- function(y, x, m, method) {
  N <- nrow(x); d <- ncol(x)
  x.max <- colMaxs(x); x.min <- colMins(x)
  x <- scale(x, center=x.min, scale=x.max-x.min)

  p   <- m + 4
  ind <- genIndex(m, d, method)
  x.p <- gpoly(x, p, rep(0, d))
  est <- lm(y ~ x.p - 1)
  beta <- est$coefficients; est <- est$fitted.values
  coef.ind <- which(!is.na(beta))

  #bias constant
  imse.b <- 0
  if (method == "bs") {
    cons.B <- abs(bernoulli(2*m+2, 0)) / factorial(2*m+2)
  } else if (method == "wav") {
    cons.B <- 1/factorial(m+1)^2 * filelist$cwav[m]
  } else {
    cons.B <-c()
  }
  if (method == "pp") {
    for (j in 1:nrow(ind)) {
      cons.B <- prod(1 / (2*ind[j,] + 1) / factorial(ind[j,])^2 / choose(2*ind[j,], ind[j,])^2)
      imse.b <- imse.b + cons.B * mean((gpoly(x, p, ind[j,])[,coef.ind, drop=F] %*% beta[coef.ind])^2)
    }
  } else {
    for (j in 1:nrow(ind)) {
      imse.b <- imse.b + cons.B * mean((gpoly(x, p, ind[j,])[,coef.ind, drop=F] %*% beta[coef.ind])^2)
    }
  }

  # variance constant
  s2 <- mean(lm(y^2 ~ x.p - 1)$fitted.values - est^2)

  if (method == "pp")  cons.V <- choose(m+d, m)
  else                 cons.V <- 1
  imse.v <- cons.V * s2

  k.rot <- ceiling((imse.b*2*(m+1)/(d*imse.v))^(1/(2*m+2+d)) * N^(1/(2*m+2+d)))

  return(k.rot)
}

getBeta <- function(y, x, method, m, smooth, knot, d) {
  P <- genDesign(x, method, m=m, J=NULL, smooth=smooth, knot=knot, deriv=rep(0, d))
  inv <- qrXXinv(P)
  beta <- inv %*% crossprod(P, y)
  return(beta)
}

genB <- function(y, x, d, method, m, deriv, knot, proj) {
  pos <- locate.x(x, d, knot)
  B <- 0
  if (method == "bs") {
    index.m <- diag(m+1, d)
    index <- which(rowSums(sweep(index.m, 2, deriv) >= 0) == d)
    beta <- getBeta(y, x, method="bs", m=m+1, smooth=m, knot=knot, d=d)
    if (length(index) != 0) {
      for (j in index) {
        order <- m + 1 - deriv[j]
        B <- B + bernoulli(order, (x[,j] - pos$tx[,j]) / pos$hx[,j]) * pos$hx[,j]^order / factorial(order) *
             genDesign(x, method="bs", m=m+1, J=NULL, smooth=m, knot=knot, deriv=index.m[j,]) %*% beta
      }
    }
    if (proj) {
      bias.0 <- 0
      for (j in 1:d) {
        bias.0 <- bias.0 + bernoulli(m+1, (x[,j] - pos$tx[,j]) / pos$hx[,j]) * pos$hx[,j]^(m+1) / factorial(m+1) *
                  genDesign(x, method="bs", m=m+1, J=NULL, smooth=m, knot=knot, deriv=index.m[j,]) %*% beta
      }
      beta <- getBeta(bias.0, x, method="bs", m=m, smooth=m-1, knot=knot, d=d)
      B <- B - genDesign(x, method="bs", m=m, J=NULL, smooth=m-1, knot=knot, deriv=deriv) %*% beta
    }
  } else if (method == "pp") {
    index.m <- t(xsimplex(d, m+1))
    index <- which(rowSums(sweep(index.m, 2, deriv) >= 0) == d)
    beta <- getBeta(y, x, method="pp", m=m+1, smooth=NULL, knot=knot, d=d)
    if (length(index) != 0) {
      for (j in index) {
        order <- index.m[j,] - deriv
        B <- B + legpoly(order, (x - pos$tx) / pos$hx) * colProds(t(pos$hx)^order) /
                 prod(factorial(order)) / prod(choose(2 * order, order)) *
                 genDesign(x, method="pp", m=m+1, J=NULL, smooth=NULL, knot=knot, deriv=index.m[j,]) %*% beta
      }
    }
    if (proj) {
      bias.0 <- 0
      for (j in 1:nrow(index.m)) {
        order0 <- index.m[j,]
        bias.0 <- bias.0 + legpoly(order0, (x - pos$tx) / pos$hx) * colProds(t(pos$hx)^order0) /
                           prod(factorial(order0)) / prod(choose(2 * order0, order0)) *
                           genDesign(x, method="pp", m=m+1, J=NULL, smooth=NULL, knot=knot, deriv=order0) %*% beta
      }
      beta <- getBeta(bias.0, x, method="pp", m=m, smooth=NULL, knot=knot, d=d)
      B <- B - genDesign(x, method="pp", m=m, J=NULL, smooth=NULL, knot=knot, deriv=deriv) %*% beta
    }
  }
  B <- mean(B^2)
  return(B)
}

genV <- function(y, x, d, method, m, J, knot, deriv, vce) {
  P <- genDesign(x, method, m, J, smooth=m-1, knot, rep(0, d))
  invG.p <- qrXXinv(P)
  basis.p <- genDesign(x, method, m, J, smooth=m-1, knot, deriv)
  hii.p <- rowSums((P %*% invG.p) * P); hii.p[hii.p >= 0.9] <- 0.9
  hii.p[hii.p < 0] <- 0; d.p <- ncol(P)
  predict.p <- P %*% invG.p %*% crossprod(P, y)
  res.p <- lsprobust.res(y, predict.p, hii.p, vce, d.p)
  se.cl <- sqrt(rowSums((basis.p %*% invG.p %*% lsprobust.vce(P, res.p) %*% invG.p) * basis.p))
  return(se.cl)
}

# DPI selector
lspkselect.imse.dpi <- function(y, x, m, method, ktype, vce, deriv, proj) {
  k.rot <- lspkselect.imse.rot(y, x, m, method)
  N <- nrow(x); d <- ncol(x); q <- sum(deriv)
  x.max <- colMaxs(x); x.min <- colMins(x)
  x <- scale(x, center=x.min, scale=x.max-x.min)
  if (ktype == "uni") {
    knot <- genKnot.u(rep(0, d), rep(1, d), d, rep(k.rot, d))
  } else {
    knot <- genKnot.q(x, d, rep(k.rot, d))
  }
  if (method == "wav") {
    J <- pmax(ceiling(log2(rep(k.rot, d)+1)), ceiling(log2(2*(m+2))))
  } else {
    J <- NULL
  }

  # bias constant
  # estimate deriv
  if (method == "wav") {
    ind <- genIndex(m, d, "wav")
    beta <- genDesign(x, method="bs", m=m+1, J=NULL, smooth=m, knot=knot, deriv=rep(0, d))
    beta <- lm(y ~ beta-1)$coefficients; index <- which(!is.na(beta))
    imse.b <- 0
    for (j in 1:nrow(ind)) {
      imse.b <- imse.b + mean((genDesign(x, method="bs", m=m+1, J=NULL, smooth=m,
                         knot=knot, deriv=ind[j,])[,index, drop=F] %*% beta[index])^2)
    }
    imse.b <- imse.b / factorial(m+1)^2 * filelist$cwav[m]
  } else {
    imse.b <- genB(y, x, d, method, m, deriv, knot, proj) * (k.rot+1)^(2*(m-q+1))
  }

  # variance constant
  imse.v <- mean((genV(y, x, d, method, m, J, knot, deriv, vce))^2) * N * (k.rot+1)^((-d-2*q)/d)

  k.dpi <- ceiling((imse.b*2*(m-q+1)/((d+2*q)*imse.v))^(1/(2*m+2+d)) * N^(1/(2*m+2+d)))

  return(k.dpi)
}

lssqrtm <- function(A) {
  decomp <- svd(A)
  rootA  <- decomp$u %*% diag(sqrt(decomp$d)) %*% t(decomp$u)
  return(rootA)
}

lsprobust.sup.pl <- function(num, denom, N, B, level) {
  temp.sup <- rep(NA, B)
  for (i in 1:B) {
    eps         <- matrix(rnorm(N, 0, 1), ncol = 1)
    tx          <- (num %*% eps) / denom
    temp.sup[i] <- max(abs(tx))
  }
  q <- quantile(temp.sup, level/100, na.rm=T)
  return(q)
}

lsprobust.sup.wb <- function(num1, num2, denom, res, N, B, level) {
  temp.sup <- rep(NA, B)
  for (i in 1:B) {
    eta <- rbinom(N, 1, 0.5) * 2 - 1
    res.wb <- eta * res
    tx <- (num1 %*% crossprod(num2, res.wb)) / denom
    temp.sup[i] <- max(abs(tx))
  }
  q <- quantile(temp.sup, level/100, na.rm=T)
  return(q)
}
