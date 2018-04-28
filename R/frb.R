#' Function to run the Fast and Robust Bootstrap
#'
#' This function runs a Fast and Robust Bootstrap for robust regression
#' estimators (MM) as computed by \code{robustbase::lmrob}.
#'
#' The fast and robust bootstrap as described in Salibian-Barrera, M. and
#' Zamar, R.H. (2002), and Salibian-Barrera, M., Van Aels, S. and Willems, G.
#' (2008).
#'
#' @param lmrob.object An object of class \code{lmrob} as returned by
#' \code{robustbase::lmrob}
#' @param nboot An integer number. The number of bootstrap samples to be used.
#' @return An object of class 'frb' which includes a matrix of the bootstrap
#'   parameter estimates.
#' @note See the github repository \link{https://github.com/msalibian/FRB}
#' @author Matias Salibian-Barrera, matias@@stat.ubc.ca
#' @references Salibian-Barrera, M. and Zamar, R.H. (2002). Bootstrapping
#' robust estimates of regression. The Annals of Statistics, 30, 556-582.
#' http://dx.doi.org/10.1214/aos/1021379865
#'
#' Salibian-Barrera, M., Van Aels, S. and Willems, G. (2008). Fast and robust
#' bootstrap. Statistical Methods and Applications 17, 41-71.
#' http://dx.doi.org/10.1007/s10260-007-0048-6
#' @keywords robustness robust regression bootstrap
#' @examples
#'
#' library(robustbase)
#' a <- lmrob(LNOx ~ LNOxEm + sqrtWS, data=NOxEmissions)
#' tmp <- frb(lmrob.object=a, nboot=1000, return.coef=FALSE)
#' # Estimated SE's for estimated regression coefficients
#' sqrt(diag(tmp))
#' # [1] 0.056422169 0.007782671 0.012662991
#' # compare with SE's based on the asymptotic approximation
#' sqrt(diag(summary(a)$cov))
#' # (Intercept)      LNOxEm      sqrtWS
#' # 0.054256788 0.007482346 0.013222502
#'
#' @export frb
frb <- function(lmrob.object, nboot=1000){
  thecall <- match.call()

  lmrob.Chi <- robustbase::Mchi
  lmrob.Psi <- robustbase::Mpsi

  co <- lmrob.object$control
  
  tuning.psi <- co$tuning.psi
  tuning.chi <- co$tuning.chi
  
  psi <- co$psi

  uu <- model.frame(lmrob.object)
  yy <- as.vector( model.extract(uu, 'response') )
  xx <- model.matrix(lmrob.object)
  n <- nrow(xx)
  p <- ncol(xx)
  attributes(xx) <- NULL
  xx <- matrix(xx, nrow = n, ncol = p)
  
  scale <- lmrob.object$scale

  beta.s <- lmrob.object$init.S$coefficients
  beta.mm <- lmrob.object$coefficients
  re.mm <- as.vector(lmrob.object$residuals)
  re.s <- as.vector( yy - xx %*% beta.s )
  
  ss <- lmrob.Chi( re.s/scale, psi = psi, cc=tuning.chi, deriv=0 )
  # ss <- Chi(re.s/scale)

  # w <- Psi(r/sigma)/r
  w <- lmrob.Psi(x= re.mm / scale, psi = psi, cc=tuning.psi, deriv=0) / re.mm
  
  # w.p <- Psi'(r/sigma)
  w.p <- lmrob.Psi(x = re.mm / scale, psi = psi, cc=tuning.psi, deriv=1)
  # a more efficient way of doing this?
  a <- t(xx) %*% diag(w.p) %*% xx
  # a more efficient way of doing this?
  b <- t(xx) %*% diag(w) %*% xx
  # w.pp <- Psi'(r/sigma)*r/sigma
  w.pp <- w.p * re.mm / scale
  d <- as.vector( t(xx) %*% w.pp )
  # e <- \sum Chi'(r/sigma)*r/sigma
  e <- sum( lmrob.Chi(re.s / scale, psi = psi, cc=tuning.chi, deriv=1) * re.s / scale )
  d <- d / 2 * (n - p) * scale / e;
  x2 <- solve(a)
  x3 <- x2 %*% b * scale
  v2 <- x2 %*% d # / scale
  # correction matrix is in x3
  # correction vector is in v2
  # set.seed(seed)

  boot.beta <- matrix(0, nrow = nboot, ncol = p)
  a <- .C('R_frb', as.double(xx), as.double(yy), as.double(w), as.integer(n),
          as.integer(p), as.double(beta.mm), as.double(scale),
          as.double(ss), bb=as.double(boot.beta), as.integer(nboot),
          as.double(x3), as.double(v2), PACKAGE="FRB")$bb

  a <- matrix(a, nrow = nboot, ncol = p)
  a <- t(t(a) + lmrob.object$coefficients)
  colnames(a) <- names(lmrob.object$coefficients)

  res <- list(call = thecall, B = nboot, param = a, lmrob.object = lmrob.object)

  class(res) <- "frb"

  res
}

#' @method print frb
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

#' @method summary frb
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

#' @method print summary.frb
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

  invisible(x)
}

#' @method plot frb
plot.frb <- function(x, y=NULL, ...){
  x <- x$param
  np <- ncol(x)
  cn <- colnames(x)
  
  oldPar <- par(no.readonly=TRUE)
  on.exit(par(oldPar))
  
  par(mfrow=c(1, np))
  
  for(i in 1:np){
    hist(x[, i], xlab=cn[i], main="")
  }
  
  invisible()
}

#' @method ggplot frb
ggplot.frb <- function(data, mapping=NULL, fill="blue", color="blue", bins=30, xlab="", ylab="count", ..., environment=parent.frame()){
  x <- data$param
  nms <- rep(colnames(x), each=nrow(x))
  d <- data.frame(value = c(x), names = nms)
  
  ggplot(d, aes(value)) +
    geom_histogram(fill=fill, color=color, bins=bins) +
    facet_wrap(~names, scales="free")
}


