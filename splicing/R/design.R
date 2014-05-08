
rdirichlet <- function (n, alpha) {
  l <- length(alpha)
  x <- matrix(rgamma(l * n, alpha), ncol = l, byrow = TRUE)
  sm <- x %*% rep(1, l)
  x / as.vector(sm)
}

sim_expr_prof <- function(n, expr, bio_var=.2) {
  bio_var <- log((1+sqrt(1+4*(bio_var/expr)^2))/2)
  x <- exp(t(rmvnorm(n, sigma=diag(bio_var))))

  x <- t(expr * x)
  x / rowSums(x)
}
