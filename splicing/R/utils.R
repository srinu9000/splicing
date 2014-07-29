  ## We have our own version of this, so that we don't need
  ## to depend on another package, and more importantly, we
  ## can zap the small (negative) eigenvalues in the decomposition
  rmvnorm <- function (n, mean = rep(0, nrow(sigma)),
                       sigma = diag(length(mean))) {
    
    assert_that((isSymmetric(sigma, tol = sqrt(.Machine$double.eps), 
                             check.attributes = FALSE)))
    assert_that((length(mean) == nrow(sigma)))

    sigma1 <- sigma
    dimnames(sigma1) <- NULL

    ev <- eigen(sigma, symmetric = TRUE)
    ev$values <- zapsmall(ev$values)
    if (!all(ev$values >=
             -sqrt(.Machine$double.eps) * abs(ev$values[1]))) {
      warning("sigma is numerically not positive definite")
    }
    retval <- ev$vectors %*% diag(sqrt(ev$values), length(ev$values)) %*% 
      t(ev$vectors)

    retval <- matrix(rnorm(n * ncol(sigma)), nrow = n, byrow = TRUE) %*% 
      retval
    retval <- sweep(retval, 2, mean, "+")
    colnames(retval) <- names(mean)
    retval
  }

check_for_package <- function(pkg, message = "") {
  res <- try(do.call("library", list(pkg)), silent = TRUE)
  if (inherits(res, "try-error")) {
    stop("Cannot load package ", pkg, message, call. = FALSE)
  }
}

pkg_fun <- getExportedValue
