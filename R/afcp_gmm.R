#' Estimates AFCPs by GMM under the restriction that direct and indirect preferences agree
#'
#' For each non-baseline level `a` of `attribute`, this estimates \eqn{AFCP(a, b)} (with `b` the baseline)
#' by over-identified efficient GMM. In addition to the direct moment condition that \eqn{\theta} equals the
#' mean choice in tasks comparing `a` and `b`, it adds one moment condition for each other level `c` imposing
#' \eqn{\theta - 1/2 = AFCP(a,c) - AFCP(b,c)}. When these restrictions hold, the GMM estimator is typically more
#' precise than the simple AFCP estimator returned by [afcp()]; when they fail it is biased. The Hansen J-test of
#' the over-identifying restrictions is reported for each level.
#'
#' Every moment condition is linear in \eqn{\theta} and has the form \eqn{\theta - m_k}, where \eqn{m_k} is a
#' target built from the pairwise AFCP estimates: \eqn{m_0 = \widehat{AFCP}(a,b)} and
#' \eqn{m_k = 1/2 + \widehat{AFCP}(a,c_k) - \widehat{AFCP}(b,c_k)}. Each pairwise mean enters the task-level
#' moments through its influence function \eqn{(n/n_{ac}) 1[ac]_i (Y_i - \widehat{AFCP}(a,c))}, so the covariance of
#' the targets, \eqn{\hat V}, is the same cluster-robust (by respondent) covariance used by [afcp()]. The
#' efficient weight matrix \eqn{W = \hat V^{-1}} does not depend on \eqn{\theta}, so the efficient GMM estimator is
#' available in closed form (equivalently, optimal minimum distance): \eqn{\hat\theta = 1'W\hat m / 1'W1} with
#' standard error \eqn{(1'W1)^{-1/2}}. The J-statistic \eqn{(\hat m - \hat\theta 1)'W(\hat m - \hat\theta 1)} is
#' numerically identical to the joint Wald statistic `wald_stat_all` from [afcp()].
#'
#' Estimates are differences in means within tasks, computed from the data in `cjointobj`; the `design` and
#' `weights` passed to `cjoint::amce()` are not stored in `cjointobj` and so are not used, and survey weights are not
#' supported. The over-identifying restrictions assume `attribute` is randomized independently of the
#' other attributes. If its randomization is restricted (e.g. some of its levels cannot appear with certain levels of
#' another attribute), each pairwise AFCP averages over a different distribution of the other attributes, and direct
#' and indirect preferences can differ even when preferences are transitive.
#'
#' @importFrom stats pnorm pchisq qnorm lm coef
#'
#' @inheritParams afcp
#'
#' @return A list containing
#' - `attribute`, `baseline` - The (cleaned) attribute and baseline level
#' - `afcp` - A dataframe containing GMM estimates of the AFCPs relative to the selected baseline level
#' - `jtest` - A dataframe containing the Hansen J-test of the over-identifying restrictions (direct = indirect) for each level
#' @export
#'
afcp_gmm <- function(cjointobj, respondent.id, task.id, profile.id, attribute, baseline = NULL, ci = .95, vcov_type = "HC2"){

  # Shared input checks and cleaning (see utils.R)
  setup <- afcp_setup(cjointobj, respondent.id, task.id, profile.id, attribute, baseline, ci)
  baseline <- setup$baseline

  gmm_results <- list()
  j_tests <- list()

  for (level in setup$attr_use){
    # Make the data wide - one row per task containing `level` or `baseline`
    wide_data <- make.wide.data(indata=setup$data, attr_var = setup$attribute, level_a = level, level_b = baseline,
                                respondentID = setup$respondent.id, choice = setup$choice.outcome,
                                qid = setup$task.id, option = setup$profile.id)

    fit <- gmm_fit_influence(wide_data, level, baseline, vcov_type)

    out_results <- data.frame(level = level, baseline = baseline, afcp = fit$theta, se = fit$se)
    out_results$zstat <- (out_results$afcp - .5)/out_results$se
    out_results$pval <- 2*pnorm(-abs(out_results$zstat))
    out_results$conf_high <- fit$theta + abs(qnorm((1-ci)/2))*fit$se
    out_results$conf_low <- fit$theta - abs(qnorm((1-ci)/2))*fit$se
    out_results$conf_level <- ci
    gmm_results[[level]] <- out_results

    j_tests[[level]] <- data.frame(level = level, baseline = baseline, L_other = fit$L_other,
                                   j_stat = fit$j_stat,
                                   j_p = if (fit$L_other > 0) pchisq(fit$j_stat, fit$L_other, lower.tail = FALSE) else NA_real_)
  }

  afcp_estimates <- bind_rows(gmm_results)
  jtest_batch <- bind_rows(j_tests)
  rownames(afcp_estimates) <- NULL
  rownames(jtest_batch) <- NULL

  return(list(attribute = setup$attribute, baseline = baseline, afcp = afcp_estimates, jtest = jtest_batch))
}


#' Pair indicators for the wide data: a-b, then a-c and b-c for each other level c (sorted)
#'
#' make.wide.data() orients tasks so level_a is first, then level_b; other levels c are second.
#' @noRd
gmm_pairs <- function(wide_data, level_a, level_b){
  v1 <- as.character(wide_data$val1)
  v2 <- as.character(wide_data$val2)
  other_levels <- sort(setdiff(unique(c(v1, v2)), c(level_a, level_b)))
  pairs <- c(list(c(level_a, level_b)),
             lapply(other_levels, function(l) c(level_a, l)),
             lapply(other_levels, function(l) c(level_b, l)))
  ind <- sapply(pairs, function(p) v1 == p[1] & v2 == p[2])
  ind <- matrix(ind, nrow = nrow(wide_data))
  if (any(is.na(ind)) || any(rowSums(ind) != 1)){
    stop("Error: some tasks do not map to exactly one level comparison")
  }
  empty <- colSums(ind) == 0
  if (any(empty)){
    p <- pairs[[which(empty)[1]]]
    stop(paste("Error: no tasks compare ", p[1], " and ", p[2], sep=""))
  }
  list(ind = ind, L_other = length(other_levels))
}

#' Efficient GMM with influence-function moments (optimal minimum distance)
#'
#' Targets: m_0 = AFCP(a,b); m_k = 1/2 + AFCP(a,c_k) - AFCP(b,c_k). With V the cluster-robust covariance of the
#' pairwise means (as in afcp()), W = (C V C')^{-1}.
#' @noRd
gmm_fit_influence <- function(wide_data, level_a, level_b, vcov_type){
  pairs <- gmm_pairs(wide_data, level_a, level_b)
  L_other <- pairs$L_other

  # Cell-means regression: coefficients are the pairwise AFCPs
  pair_ind <- pairs$ind * 1 # numeric: lm() treats a logical matrix as a factor
  fit <- lm(wide_data$choose ~ pair_ind - 1)
  V <- sandwich::vcovCL(fit, cluster = as.character(wide_data$respid), type = vcov_type) # as.character: vcovCL's HC2 is wrong when the cluster is a factor with unused levels

  C <- matrix(0, nrow = L_other + 1, ncol = 2*L_other + 1)
  C[1, 1] <- 1
  for (k in seq_len(L_other)){
    C[1 + k, 1 + k] <- 1
    C[1 + k, 1 + L_other + k] <- -1
  }
  m <- c(C %*% coef(fit)) + c(0, rep(.5, L_other))
  W <- solve_moment_vcov(C %*% V %*% t(C))
  ones <- rep(1, L_other + 1)

  info <- c(crossprod(ones, W %*% ones))
  theta <- c(crossprod(ones, W %*% m))/info
  resid <- m - theta
  j_stat <- if (L_other > 0) c(crossprod(resid, W %*% resid)) else NA_real_

  list(theta = theta, se = sqrt(1/info), j_stat = j_stat, L_other = L_other)
}

#' @noRd
solve_moment_vcov <- function(S){
  tryCatch(solve(S), error = function(e){
    stop("Error: covariance matrix of the GMM moments is singular -- some level comparisons may have too few tasks")
  })
}
