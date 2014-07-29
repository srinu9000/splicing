
power_pooled <- function(gene_structure, gene=1, isoform, coverage, no_reads,
                        dom_expr=2/3, supp_expr=1/3, read_len=33,
                        no_trials=100, sig_level=0.05) {

  assert_that(is.count(gene))
  assert_that(is.count(isoform))
  assert_that(!missing(coverage) || !missing(no_reads))
  if (!missing(coverage)) {
    assert_that(is.numeric(coverage), length(coverage) == 1,
                coverage > 0)
  }
  if (!missing(no_reads)) {
    assert_that(is.count(no_reads), no_reads >= 1)
  }
  assert_that(length(dom_expr) == 1, is.numeric(dom_expr),
              dom_expr > 0, dom_expr < 1)
  assert_that(length(supp_expr) == 1, is.numeric(supp_expr),
              supp_expr > 0, supp_expr < 1)
  assert_that(is.count(read_len))
  assert_that(is.count(no_trials))

  if (!missing(no_reads) && !missing(coverage)) {
    warning("'no_reads' given, so 'coverage' is ignored")
  }

  mygene <- selectGenes(gene_structure, gene)
  noiso <- noIso(mygene)
  assert_that(noiso > 1, isoform <= noiso)

  if (missing(no_reads)) {
    no_reads <- coverage_to_no_reads(gene_structure, gene=gene,
                                     exp=rep(1/noiso, noiso),
                                     coverage=coverage,
                                     readLength=read_len)
  }


  fail <- list(p=rep(1, no_trials), power=0)

  mat <- assignmentMatrix(mygene, readLength=read_len)
  elen <- rowSums(mat)
  assert_that(are_equal(elen, isoLength(mygene)[[1]] - read_len + 1))

  mat <- t(mat / elen)
  noc <- nrow(mat)

  ## Expression profiles
  expr_dom <- rep((1-dom_expr) / (noiso-1), noiso)
  expr_dom[isoform] <- dom_expr
  assert_that(all.equal(sum(expr_dom), 1))

  expr_supp <- rep((1-supp_expr) / (noiso-1), noiso)
  expr_supp[isoform] <- supp_expr
  assert_that(all.equal(sum(expr_supp), 1))

  ## Check if the given isoform has unique reads.
  ## If it does not, then we don't have much to do
  uniq_m <- which(rowSums(mat != 0) == 1)
  uniq_iso <- apply(mat[uniq_m, , drop=FALSE] != 0, 1, which)
  if (!isoform %in% uniq_iso) { return(fail) }

  matching <- match(names(uniq_iso)[isoform == uniq_iso], rownames(mat))
  assert_that(length(matching) == 1)

  ## We also need to check if there is a class that is
  ## definitely not matching the selected isoform.
  if (all(mat[, isoform] > 0)) { return(fail) }

  ## The odds for sampling (1) unique reads from the isoform
  ## (2) reads definitely not matching the isoform, and
  ## (3) ambiguous reads
  p_dom <- mat %*% expr_dom
  p_supp <- mat %*% expr_supp

  non_matching <- which(mat[, isoform] == 0)
  sample_odds_dom <- c(p_dom[matching],
                       sum(p_dom[non_matching]),
                       sum(p_dom[setdiff(seq_len(noc),
                                         c(matching, non_matching))]))
  assert_that(all.equal(sum(sample_odds_dom), 1))

  sample_odds_supp <- c(p_supp[matching],
                        sum(p_supp[non_matching]),
                        sum(p_supp[setdiff(seq_len(noc),
                                           c(matching, non_matching))]))
  assert_that(all.equal(sum(sample_odds_supp), 1))

  ## OK, now sample 'no_trials' contingency tables and use
  ## Fisher's exact test on the matching/non-matching part

  con_table <- function(odds_dom, odds_supp) {
    dom  <- sample(1:3, no_reads, replace=TRUE, prob=odds_dom)
    supp <- sample(1:3, no_reads, replace=TRUE, prob=odds_supp)
    t_dom <- table(factor(dom, levels=1:3))
    t_supp <- table(factor(supp, levels=1:3))
    cbind(t_dom, t_supp)[1:2,]
  }

  pvals <- sapply(1:no_trials, function(i) {
    tab <- con_table(sample_odds_dom, sample_odds_supp)
    fisher.test(tab)$p.value
  })

  list(p = pvals, power = sum(pvals < sig_level) / length(pvals))
}

power_pooled_all_iso <- function(gene_structure, gene=1, ...) {
  noiso <- noIso(gene_structure)[gene]
  lapply(1:noiso, function(i) {
    power_pooled(gene_structure=gene_structure, gene=gene, isoform=i, ...)
  })
}

## This will be (probably) easy. We just simulate "solutions"
## for the gene: first we sample an expression profile, using
## the assumed biological variance. Then we use this expression
## profile to generate means and variances for the solution,
## using our CR bound. We assume that the distribution of the
## estimator (for a given expression profile) is normal around the
## truth, with the distribution given by the CR bound. So we just
## sample from this normal.

power_rep <- function(gene_structure, gene=1, isoform, coverage, no_reads,
                     dom_expr=2/3, supp_expr=1/3, read_len=33,
                     no_trials=100, sig_level=0.05, bio_var=.2,
                     n_dom=3, n_supp=3, n_sim=100, n_sim_sampling=10) {

  assert_that(is.count(gene))
  assert_that(is.count(isoform))
  assert_that(!missing(coverage) || !missing(no_reads))
  if (!missing(coverage)) {
    assert_that(is.numeric(coverage), length(coverage) == 1,
                coverage > 0)
  }
  if (!missing(no_reads)) {
    assert_that(is.count(no_reads), no_reads >= 1)
  }
  assert_that(length(dom_expr) == 1, is.numeric(dom_expr),
              dom_expr > 0, dom_expr < 1)
  assert_that(length(supp_expr) == 1, is.numeric(supp_expr),
              supp_expr > 0, supp_expr < 1)
  assert_that(is.count(read_len))
  assert_that(is.count(no_trials))
  assert_that(is.numeric(sig_level), length(sig_level) == 1,
              is.finite(sig_level), sig_level > 0, sig_level < 1)
  assert_that(is.numeric(bio_var), length(bio_var) == 1,
              is.finite(bio_var), bio_var >= 0)
  assert_that(is.count(n_dom), n_dom > 1)
  assert_that(is.count(n_supp), n_supp > 1)
  assert_that(is.count(n_sim), n_sim >= 1)
  assert_that(is.count(n_sim_sampling), n_sim_sampling >= 1)

  if (!missing(no_reads) && !missing(coverage)) {
    warning("'no_reads' given, so 'coverage' is ignored")
  }

  mygene <- selectGenes(gene_structure, gene)
  noiso <- noIso(mygene)
  assert_that(noiso > 1, isoform <= noiso)

  if (missing(no_reads)) {
    no_reads <- coverage_to_no_reads(gene_structure, gene=gene,
                                     exp=rep(1/noiso, noiso),
                                     coverage=coverage,
                                     readLength=read_len)
  }

  ## Cache un-normalized assignment matrix
  amat <- assignmentMatrix(gene_structure, gene = gene,
                          readLength = read_len, overHang = 1L,
                          paired = FALSE, fast = FALSE,
                          fragmentProb = NULL, fragmentStart = 0L,
                          normalMean = NA, normalVar = NA, numDevs = 4)

  ## Expression profiles
  expr_dom <- rep((1-dom_expr) / (noiso-1), noiso)
  expr_dom[isoform] <- dom_expr
  assert_that(all.equal(sum(expr_dom), 1))

  expr_supp <- rep((1-supp_expr) / (noiso-1), noiso)
  expr_supp[isoform] <- supp_expr
  assert_that(all.equal(sum(expr_supp), 1))

  do_var <- function(bio_expr) {
    cr <- lapply(seq_len(nrow(bio_expr)), function(i) {
      crMatrix(gene_structure, gene=gene, readLength=read_len,
               expr=bio_expr[i,], assignmentMatrix = amat)
    })
    vv <- lapply(seq_len(nrow(bio_expr)), function(i) {
      rmvnorm(n_sim_sampling, mean=bio_expr[i,], sigma=cr[[i]] / no_reads)
    })
    apply(do.call(rbind, vv), 2, var)
  }

  bio_expr_dom  <- sim_expr_prof(n_sim, expr_dom,  bio_var=bio_var)
  bio_expr_supp <- sim_expr_prof(n_sim, expr_supp, bio_var=bio_var)

  var_dom  <- do_var(bio_expr_dom)[isoform]
  var_supp <- do_var(bio_expr_supp)[isoform]

  ## Power of a two-sample t-test with unequal variance
  dW1 <- (var_dom / n_dom + var_supp / n_supp) ^2
  dW2 <- (var_dom^2 / n_dom^2 / (n_dom-1) +
          var_supp^2 / n_supp^2 / (n_supp-1))
  dW <- dW1 / dW2
  lambdaW <- (dom_expr - supp_expr) / sqrt(var_dom/n_dom + var_supp/n_supp)

  qu <- qt(sig_level / 2, df = dW, lower.tail = FALSE)
  (pt(qu, df = dW, ncp = lambdaW, lower.tail = FALSE) +
   pt(-qu, df = dW, ncp = lambdaW, lower.tail = TRUE))
}

power_rep_all_iso <- function(gene_structure, gene = 1, ...) {
  noiso <- noIso(gene_structure)[gene]
  sapply(1:noiso, function(i) {
    power_rep(gene_structure=gene_structure, gene=gene, isoform=i, ...)
  })
}

## Trend in time series. We do a linear regression against time steps.
## The question is whether to do a weighted regression or not.
## For weighted regression it is likely that we need to use
## simulations. For the regular regression, there are some
## analytical results, mostly from the Cohen book.

## Let's do the simple regression first.

## TODO: we can also do multiple 'n's at once

power_ts <- function(gene_structure, gene=1, isoform, coverage, no_reads,
                    n=3, no_timepoints=5, iso_expr=seq(1/4, 3/4,
                                            length.out=no_timepoints),
                    read_len=33, sig_level=0.05, bio_var=.2,
                    n_sim=100, n_sim_sampling=10) {

  assert_that(is.count(gene))
  assert_that(is.count(isoform))
  assert_that(!missing(coverage) || !missing(no_reads))
  if (!missing(coverage)) {
    assert_that(is.numeric(coverage), length(coverage) == 1,
                coverage > 0)
  }
  if (!missing(no_reads)) {
    assert_that(is.count(no_reads), no_reads >= 1)
  }
  assert_that(is.count(n))
  assert_that(is.count(no_timepoints), no_timepoints >= 2)
  assert_that(is.numeric(iso_expr), length(iso_expr) == no_timepoints,
              all(iso_expr >= 0), all(iso_expr <= 1))
  assert_that(is.count(read_len))
  assert_that(is.numeric(sig_level), length(sig_level) == 1,
              is.finite(sig_level), sig_level > 0, sig_level < 1)
  assert_that(is.numeric(bio_var), length(bio_var) == 1,
              is.finite(bio_var), bio_var >= 0)

  if (!missing(no_reads) && !missing(coverage)) {
    warning("'no_reads' given, so 'coverage' is ignored")
  }

  mygene <- selectGenes(gene_structure, gene)
  noiso <- noIso(mygene)
  assert_that(noiso > 1, isoform <= noiso)

  if (missing(no_reads)) {
    no_reads <- coverage_to_no_reads(gene_structure, gene=gene,
                                     exp=rep(1/noiso, noiso),
                                     coverage=coverage,
                                     readLength=read_len)
  }

  ## Cache un-normalized assignment matrix
  amat <- assignmentMatrix(gene_structure, gene = gene,
                          readLength = read_len, overHang = 1L,
                          paired = FALSE, fast = FALSE,
                          fragmentProb = NULL, fragmentStart = 0L,
                          normalMean = NA, normalVar = NA, numDevs = 4)

  ## Expression profiles
  expr <- sapply(iso_expr, function(e) {
    ee <- rep ((1 - e) / (noiso - 1), noiso)
    ee[isoform] <- e
    ee
  })
  assert_that(all.equal(colSums(expr), rep(1, no_timepoints)))

  ## Estimate the pooled variance
  do_var <- function(bio_expr, expr) {
    cr <- lapply(seq_len(nrow(bio_expr)), function(i) {
      crMatrix(gene_structure, gene=gene, readLength=read_len,
               expr=bio_expr[i,], assignmentMatrix = amat)
    })
    vv <- lapply(seq_len(nrow(bio_expr)), function(i) {
      rmvnorm(n_sim_sampling, mean=bio_expr[i, ], sigma=cr[[i]] / no_reads)
    })
    t(t(do.call(rbind, vv)) - expr)
  }

  bio_expr <- lapply(seq_len(ncol(expr)), function(i) {
    sim_expr_prof(n_sim, expr[, i], bio_var=bio_var)
  })

  var_ts <- lapply(seq_along(bio_expr),
                   function(i) do_var(bio_expr[[i]], expr[, i] ))
  var <- apply(do.call(rbind, var_ts), 2, var)[isoform]

  ## Correlation between time steps and expression,
  ## as estimated by Spearman's correction
  var_iso <- var(iso_expr)
  r <- sqrt(var_iso/ (var_iso + var))
  f2 <- r * r / ( 1 - r * r )

  ## Power, acccoding to Cohen
  u <- 1
  v <- n * no_timepoints- u - 1

  lambda <- f2 * (u + v + 1)
  pf(qf(sig_level, u, v, lower.tail = FALSE), u, v, lambda,
     lower.tail = FALSE)
}

power_ts_all_iso <- function(gene_structure, gene=1, ...) {
  noiso <- noIso(gene_structure)[gene]
  sapply(1:noiso, function(i) {
    power_ts(gene_structure=gene_structure, gene=gene, isoform=i, ...)
  })
}
