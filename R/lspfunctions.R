#Uniform knot list
genKnot.u <- function(x.min, x.max, d, n) {
  knot <- list()
  for (j in 1:d) {
    if (n[j] == 0) knot[[j]] <- NULL
    else           knot[[j]] <- seq(x.min[j], x.max[j], length.out = n[j]+2)[2 : (n[j]+1)]
  }
  knot[d+1] <- "NULL"
  return(knot)
}

#Quantile knot list
genKnot.q <- function(x, d, n) {
  knot <- list()
  for (j in 1:d) {
    if (n[j] == 0) knot[[j]] <- NULL
    else           knot[[j]] <- quantile(x[,j], seq(0, 1, 1/(n[j]+1)))[2: (n[j]+1)]
  }
  knot[d+1] <- "NULL"
  return(knot)
}

#Generate index set
genIndex <- function(m, d, method) {
  if (method == "bs" | method=="wav") index.m <- diag(m+1, d)
  if (method == "pp")                index.m <- t(xsimplex(d, m+1))
  return(index.m)
}

#Locate the closest knots to the left of evaluating points and the length of interval
pos.x <- function(eval, x, knot) {
  d <- ncol(x); tx <- c(); hx <- c()
  for (j in 1:d) {
    knot[[j]] <- c(min(x[,j]), knot[[j]], max(x[,j]))
    ind <- sum(knot[[j]] <= eval[j])
    if (ind == length(knot[[j]])) ind <- ind - 1
    tx  <- c(tx, knot[[j]][ind])
    hx  <- c(hx, knot[[j]][ind+1] - knot[[j]][ind])
  }
  return(list(tx = tx, hx = hx))
}

#slightly change pos.x, for locpoly only, output length only, 0.1 is arbitrary
locate.h <- function(coord, d, x.min, x.max, knot) {
  hx <- c(); ehx <- c()
  for (j in 1:d) {
    knot[[j]] <- c(x.min[j], knot[[j]], x.max[j])
    ind <- sum(knot[[j]] <= coord[j])
    hx  <- c(hx, knot[[j]][ind+1] - knot[[j]][ind])
    if (ind == length(knot[[j]]) - 1) {
      ehx <- c(ehx, knot[[j]][ind+1]+0.1 - knot[[j]][ind])
    } else {
      ehx <- c(ehx, knot[[j]][ind+1] - knot[[j]][ind])
    }
  }
  return(list(hx = hx, ehx = ehx))
}

#check n of eval
gen.eN <- function(R, x, d) {
  eN <- sum(colAlls(rbind(apply(x, 1, ">=", R[1:d]), apply(x, 1, "<=", R[(d+1):(2*d)]))))
  return(eN)
}

#Daubechies scaling function
dbDesign <- function(eval, x, m, J) {
  m <- m+1
  x.max <- max(x); x.min <- min(x); n <- length(eval); k <- 2^J; h <- (x.max - x.min)/ k;
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

#generate marginal design
design.margin <- function(eval, x, method, m, J, smooth, knot, deriv) {
   if (method == "bs") {
      ext.knot <- c(rep(min(x), m+1), rep(knot, each = m-smooth), rep(max(x), m+1))
      P <- splineDesign(knots = ext.knot, eval, ord = m+1, derivs = rep(deriv, length(eval)))
   } else if (method == "wav") {
      P <- dbDesign(eval, x, m, J)
   }
   return(P)
}

#calculate power series matrix
genPower <- function(A, k) {
  B <- rowProds(sweep(A, 2, k, FUN="^"))
  return(B)
}


#modify poly function
gpoly <- function(eval, m, deriv) {
  if (nrow(eval) == 1) {
    P.level <- t(as.matrix(poly(rbind(eval, eval), degree = m, raw = T, simple = T)[1,]))
  } else {
    P.level <- poly(eval, degree = m, raw = T, simple = T)
  }
  if (!all(deriv == rep(0, ncol(eval)))) {
    d       <- ncol(eval)
    name    <- colnames(P.level)
    de      <- as.matrix(sapply(name, function(x) as.numeric(unlist(strsplit(x, "[.]")))))
    if (d == 1) de <- t(de)
    de.low  <- sweep(de, 1, deriv)
    de.low[de.low < 0] <- NA
    index              <- which(is.na(colSums(de.low)))
    P.level[,index]    <- 0
    if (length(index) != ncol(P.level)) {
      if (length(index) == 0) {
        index <- 1:ncol(P.level)
      } else {
        index <- (1:ncol(P.level))[-index]
      }
      power     <- as.matrix(de[,index])
      power.low <- as.matrix(de.low[,index])
      if (d == 1) {
        power     <- t(power)
        power.low <- t(power.low)
      }
      P.level[,index]  <- apply(power.low, 2, function(z) genPower(eval, z))
      P.deriv <- as.matrix(P.level[,index])
      if (nrow(eval) == 1) P.deriv <- t(P.deriv)
      P.level[,index]  <- sweep(P.deriv, 2,
                                colProds(factorial(power)/factorial(power.low)),
                                FUN = "*")
    }
    P.level <- cbind(rep(0, nrow(eval)), P.level)
  } else P.level <- cbind(rep(1, nrow(eval)), P.level)
  return(P.level)
}

#computation function used in locpoly
compute.locpoly <- function(coord, eval, d, d.loc, x.min, x.max, knot, m, deriv) {
  P <- matrix(0, nrow(eval), d.loc)
  mesh <- locate.h(coord, d, x.min, x.max, knot)
  index <- x.input <- sweep(eval, 2, coord)
  x.input <- sweep(x.input, 2, mesh$hx, FUN = "/")
  index   <- sweep(index, 2, mesh$ehx, FUN = "/")
  index[index < 0 | index >= 1] <- NA
  index   <- which(!is.na(rowSums(index)))
  if (length(index) != 0) {
    P[index, ] <- gpoly(matrix(x.input[index,], length(index), d), m, deriv) / prod(mesh$hx^deriv)
  }
  return(P)
}

#partitioned polynomial basis
locpoly <- function(eval, x, m, knot, deriv) {
  d   <- ncol(eval); x.max <- colMaxs(x); x.min <- colMins(x)
  seq <- list()
  for (j in 1:d) seq[[j]] <- c(min(x[,j]), knot[[j]])
  coord <- as.matrix(expand.grid(seq))
  d.loc <- choose(d + m, m)
  P <- apply(coord, 1, function(v) compute.locpoly(v, eval, d, d.loc, x.min, x.max, knot, m, deriv))
  P <- matrix(P, nrow(eval), nrow(coord) * d.loc)
  return(P)
}

#Generate tensor-product basis or partitioning basis
genDesign <- function(eval, x, method, m, J, smooth, knot, deriv) {
  if (method == "bs" | method == "wav") {
    S <- list()
    for (j in 1 : ncol(x)) {
      S[[j]] <- design.margin(eval[,j], x[,j], method, m, J[j], smooth, knot[[j]], deriv[j])
    }
    P <- tensor.prod.model.matrix(S)
  } else if (method == "pp") {
    P <- locpoly(eval, x, m, knot, deriv)
  }
  return(P)
}

#Shifted multivariate legendre polynomial
legpoly <- function(degree, x) {
  L <- 1
  for (j in 1:length(degree)) {
    if (degree[j] != 0) L <- L * sapply(x[,j], function(z) legendre(degree[j], 2*z-1)[1,])
  }
  return(L)
}

#Generate bias vector
genBias <- function(eval, x, method, m, q, knot, bknot, bsmooth, deriv) {
  d <- ncol(eval)
  pos <- matrix(unlist(apply(eval, 1, function(z) pos.x(z, x, knot))), ncol = 2 * d, byrow = T)
  pos.tx <- pos[,1:d, drop=F]
  pos.hx <- pos[,(d+1):(2*d), drop=F]
  if (method == "bs") {
    index.m <- diag(m+1, d)
    subset  <- sweep(index.m, 2, deriv)
    subset[subset < 0] <- NA
    index <- which(!is.na(rowSums(subset)))
    B <- 0
    if (length(index) != 0) {
      for (j in index) {
          order <- m + 1 - deriv[j]
          deriv.bias <- rep(0, d); deriv.bias[j] <- m+1;
          B <- B + bernoulli(order, (eval[,j] - pos.tx[,j]) / pos.hx[,j]) *
               pos.hx[,j]^order / factorial(order) *
               genDesign(eval=eval, x=x, method=method, m=q, smooth=bsmooth,
                         knot=bknot, deriv=deriv.bias)
      }
    }
  } else if (method == "pp") {
    index.m <- t(xsimplex(d, m+1))
    subset  <- sweep(index.m, 2, deriv)
    subset[subset < 0] <- NA
    index <- which(!is.na(rowSums(subset)))
    B <- 0
    if (length(index) != 0) {
      for (j in index) {
        deriv.bias <- index.m[j,]; order <- subset[j,];
        B <- B + legpoly(order, (eval - pos.tx) / pos.hx) * colProds(t(pos.hx)^order) /
                 prod(factorial(order)) / prod(choose(2 * order, order)) *
                 genDesign(eval=eval, x=x, method=method, m=q, knot=bknot, deriv=deriv.bias)
      }
    }
  }
  return(B)
}

#Generate squared bias, only for dpi function
genB <- function(y, x, xmin, xmax, method, m, kappa, deriv, ktype, vce, proj, proj.bias) {
  d <- ncol(x)
  nknot <- rep(kappa, d)
  if (ktype == "uni") knot <- genKnot.u(xmin, xmax, d, nknot)
  if (ktype == "qua") knot <- genKnot.q(x, d, nknot)
  pos <- matrix(unlist(apply(x, 1, function(z) pos.x(z, x, knot))), ncol = 2 * d, byrow = T)
  pos.tx <- pos[,1:d, drop=F]
  pos.hx <- pos[,(d+1):(2*d), drop=F]
  if (method == "bs") {
    index.m <- diag(m+1, d)
    subset  <- sweep(index.m, 2, deriv)
    subset[subset < 0] <- NA
    index <- which(!is.na(rowSums(subset)))
    B <- 0
    if (length(index) != 0) {
      for (j in index) {
        order <- m + 1 - deriv[j]
        deriv.bias <- rep(0, d); deriv.bias[j] <- m+1
        mu.bias.deriv <- lsprobust(y, x, x, method = method, m = m+2, deriv = deriv.bias, vce = vce,
                                    knot = knot, same = T)$Estimate[,"tau.cl"]
        B <- B + bernoulli(order, (x[,j] - pos.tx[,j]) / pos.hx[,j]) *
                 pos.hx[,j]^order / factorial(order) * mu.bias.deriv
      }
    }
    if (proj == T) {
      bias.0 <- 0
      for (j in 1:d) {
        order0 <- index.m[j,]
        mu.bias.0 <- lsprobust(y, x, x, method = method, m = m+2, deriv = order0, vce = vce,
                               knot = knot, same = T)$Estimate[,"tau.cl"]
        bias.0 <- bias.0 + bernoulli(m+1, (x[,j] - pos.tx[,j]) / pos.hx[,j]) *
                           pos.hx[,j]^(m+1) / factorial(m+1) * mu.bias.0
      }
      B <- B - proj.bias %*% bias.0
    }
  } else if (method == "pp") {
    index.m <- t(xsimplex(d, m+1))
    subset  <- sweep(index.m, 2, deriv)
    subset[subset < 0] <- NA
    index <- which(!is.na(rowSums(subset)))
    B <- 0
    if (length(index) != 0) {
      for (j in index) {
        deriv.bias <- index.m[j,]; order <- subset[j,]
        mu.bias.deriv <- lsprobust(y, x, x, method = method, m = m+2, deriv = deriv.bias, vce = vce,
                                    knot = knot, same = T)$Estimate[,"tau.cl"]
        B <- B + legpoly(order, (x - pos.tx) / pos.hx) * colProds(t(pos.hx)^order) /
                 prod(factorial(order)) / prod(choose(2 * order, order)) * mu.bias.deriv
      }
    }
    if (proj == T) {
      bias.0 <- 0
      for (j in 1:nrow(index.m)) {
        order0 <- index.m[j,]
        mu.bias.0 <- lsprobust(y, x, x, method = method, m = m+2, deriv = order0, vce = vce,
                                knot = knot, same = T)$Estimate[,"tau.cl"]
        bias.0 <- bias.0 + legpoly(order0, (x - pos.tx) / pos.hx) * colProds(t(pos.hx)^order0) /
                           prod(factorial(order0)) / prod(choose(2 * order0, order0)) * mu.bias.0
      }
      B <- B - proj.bias %*% bias.0
    }
  }
  B <- mean(B^2)
  return(B)
}

#Generate residuls
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

#Generate ROT with global polynomials (for rot, do not provide deriv option,
#since rates are the same as levels)
lspkselect.imse.rot <- function(y, x, m, method) {
  N <- nrow(x); d <- ncol(x)
  x.max <- colMaxs(x); x.min <- colMins(x)
  z <- sweep(sweep(x, 2, x.min), 2, (x.max - x.min), FUN="/")

  #(m+1)th deriv of g
  p   <- m + 4
  ind <- genIndex(m, d, method)
  z.p <- gpoly(x, p, rep(0, d))
  beta <- lm(y ~ z.p - 1)
  coef.ind <- which(!is.na(beta$coefficients))
  g.m.hat <- matrix(NA, N, nrow(ind))
  cons.B <- c()
  for (j in 1:nrow(ind)) {
    g.m.hat[,j] <- gpoly(x, p, ind[j,])[,coef.ind] %*% beta$coeff[coef.ind]
    if (method == "pp") cons.B[j]  <- prod(1 / (2*ind[j,] + 1) / factorial(ind[j,])^2 /
                                             choose(2*ind[j,], ind[j,])^2)
  }
  if (method == "bs")  cons.B <- abs(bernoulli(2*m+2, 0)) / factorial(2*m+2)
  if (method == "wav") cons.B <- 1/factorial(m+1)^2 * filelist$cwav[m]

  #bias constant
  imse.b <- sum(colMeans(g.m.hat^2) * cons.B)

  #variance constant
  beta2 <- lm(y^2 ~ z.p - 1)
  s2 <- mean(beta2$fitted.values - (beta$fitted.values)^2)

  if (method == "pp") cons.V <- choose(m+d, m)
  else                 cons.V <- 1
  imse.v <- cons.V * s2

  k.rot <- ceiling((imse.b*2*(m+1)/(d*imse.v))^(1/(2*m+2+d)) * N^(1/(2*m+2+d)))
  return(k.rot)
}

lspkselect.imse.dpi <- function(y, x, m, method, ktype, vce, deriv, proj) {
  k.rot <- lspkselect.imse.rot(y, x, m, method)
  N <- nrow(x); d <- ncol(x)
  x.max <- colMaxs(x); x.min <- colMins(x)
  z <- sweep(sweep(x, 2, x.min), 2, (x.max - x.min), FUN="/")

  #estimate deriv
  if (method == "wav") {
    ind <- genIndex(m, d, "wav")
    g.m.hat <- matrix(NA, N, nrow(ind))
    for (j in 1:nrow(ind)) {
      g.m.hat[,j] <- lsprobust(y, z, z, method = "bs", m = m+2, deriv = ind[j,], vce = vce,
                               nknot = rep(k.rot, d), ktype = ktype, same = T)$Estimate[,"tau.cl"]
    }
  }

  g.deriv.se <- lsprobust(y, z, z, method = method, m = m+1, deriv = deriv, vce = vce,
                          nknot = rep(k.rot, d), ktype = ktype, same = T)$Estimate[,"se.cl"]

  q <- sum(deriv)

  #bias constant
  if (method != "wav") {
    proj.bias <- NULL
    if (proj == T) {
      if (ktype == "uni") knot <- genKnot.u(rep(0,d), rep(1,d), d, rep(k.rot, d))
      if (ktype == "qua") knot <- genKnot.q(z, d, rep(k.rot, d))
      P <- genDesign(z, z, method=method, m=m, smooth=m-1, knot=knot, deriv=rep(0, d))
      invG.p <- qrXXinv(P)
      basis.p <- genDesign(z, z, method=method, m=m, smooth=m-1, knot=knot, deriv=deriv)
      proj.bias <- basis.p %*% invG.p %*% t(P)
    }
    imse.b <- genB(y, z, rep(0, d), rep(1, d), method, m, k.rot, deriv, ktype, vce,
                   proj=proj, proj.bias=proj.bias) * (k.rot+1)^(2*(m-q+1))
  } else  {
    imse.b <- sum(colMeans(g.m.hat^2)) / factorial(m+1)^2 * filelist$cwav[m]
  }

  #variance constant
  imse.v <- mean(g.deriv.se^2) * N * (k.rot+1)^((-d-2*q)/d)

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
    eps    <- matrix(rnorm(N, 0, 1), ncol = 1)
    tx     <- (num %*% eps) / denom
    temp.sup[i] <- max(abs(tx))
  }
  q <- quantile(temp.sup, level/100, na.rm=T)
  return(q)
}

lsprobust.sup.wb <- function(num, denom, res, N, B, level) {
  temp.sup <- rep(NA, B)
  for (i in 1:B) {
    eta <- rbinom(N, 1, 0.5) * 2 - 1
    res.wb <- eta * res
    tx <- (num %*% res.wb) / denom
    temp.sup[i] <- max(abs(tx))
  }
  q <- quantile(temp.sup, level/100, na.rm=T)
  return(q)
}
