
frb <- function(lmrob.object, nboot=1000){
  thecall <- match.call()

  lmrob.Chi <- Mchi
  lmrob.Psi <- Mpsi
  # lmrob.Psi <- tukeyPsi1
  # lmrob.Chi <- tukeyChi
  co <- lmrob.object$control
  tuning.psi <- co$tuning.psi
  tuning.chi <- co$tuning.chi
  beta.s <- lmrob.object$init.S$coefficients
  beta.mm <- coef(lmrob.object)
  re.mm <- as.vector( residuals(lmrob.object) )
  uu <- model.frame(lmrob.object)
  yy <- as.vector( model.extract(uu, 'response') )
  xx <- model.matrix(lmrob.object)
  n <- (dd <- attr(xx, 'dim'))[1]
  p <- dd[2]
  attributes(xx) <- NULL
  xx <- matrix(xx, n, p )
  re.s <- as.vector( yy - xx %*% beta.s )
  scale <- lmrob.object$scale
  # w.p <- Psi'(r/sigma)
  w.p <- lmrob.Psi(x = re.mm / scale, psi='bisquare', cc=tuning.psi, deriv=1)
  # a more efficient way of doing this?
  a <- t(xx) %*% diag(w.p) %*% xx
  # w <- Psi(r/sigma)/r
  w <- lmrob.Psi(x= re.mm / scale, psi='bisquare', cc=tuning.psi, deriv=0) / re.mm
  # a more efficient way of doing this?
  b <- t(xx) %*% diag(w) %*% xx
  # w.pp <- Psi'(r/sigma)*r/sigma
  w.pp <- w.p * re.mm / scale
  d <- as.vector( t(xx) %*% w.pp )
  # e <- \sum Chi'(r/sigma)*r/sigma
  e <- sum( lmrob.Chi(re.s / scale, psi='bisquare', cc=tuning.chi, deriv=1) * re.s / scale )
  d <- d / 2 * (n - p) * scale / e;
  x2 <- solve(a)
  x3 <- x2 %*% b * scale
  v2 <- x2 %*% d # / scale
  # correction matrix is in x3
  # correction vector is in v2
  # set.seed(seed)
  ss <- lmrob.Chi( re.s/scale, psi='bisquare', cc=tuning.chi, deriv=0 )
  # ss <- Chi(re.s/scale)
  boot.beta <- matrix(0, nboot, p)
  a <- .C('R_frb', as.double(xx), as.double(yy), as.double(w), as.integer(n),
          as.integer(p), as.double(beta.mm), as.double(scale),
          as.double(ss), bb=as.double(boot.beta), as.integer(nboot),
          as.double(x3), as.double(v2), PACKAGE="FRB")$bb

  a <- matrix(a, nboot, p)
  a <- t(t(a) + coef(lmrob.object))
  colnames(a) <- names(coef(lmrob.object))

  res <- list(call = thecall, B = nboot, param = a, lmrob.object = lmrob.object)

  class(res) <- "frb"

  res
}

print.frb <- function(x, ...){
  print(x$call)

  cat("\nOriginal estimates:\n")
  print(coef(x$lmrob.object))

  B <- x$B
  x <- x$param

  m <- apply(x, 2, quantile)
  co <- cov(x)

  cat(paste("\nSummary of", B, "bootstrap distribution:\n"))
  print(m)

  cat("\nCovariance:\n")
  print(co)

  invisible()
}

summary.frb <- function(object, prob=c(.025, .25, .5, .75, .975), ...){
  thecall <- match.call()

  x <- object$param

  co <- cov(x)
  corr <- cor(x)

  qu <- apply(x, 2, quantile, prob=prob)

  res <- list(call = thecall, quantiles = qu, cov = co, cor = corr,
              B = object$B, lmrob.object = object$lmrob.object)
  class(res) <- "summary.frb"

  res
}

print.summary.frb <- function(x, ...){
  print(x$call)

  o <- t(coef(summary(x$lmrob.object))[, 1:2])
  cat("\nOriginal estimates:\n")
  print(o)

  cat(paste("\nSummaries of", x$B, "bootstrap samples:\n"))
  cat("\nQuantiles:\n")
  print(t(x$quantiles))

  cat("\nCovariance:\n")
  print(x$cov)

  cat("\nCorrelation:\n")
  print(x$cor)

}