#' AMCE estimates for every combination of held-out attribute levels ("specification curve")
#'
#' Implements the sensitivity analysis for the inclusion of indirect comparisons. For an attribute with
#' \eqn{D} levels, there are \eqn{2^{D-2}} ways to choose which of the levels other than `level_a` and `level_b`
#' to hold out. For each combination, every task (not just every profile) containing a held-out level is dropped
#' and the AMCE of `level_a` relative to `level_b` is re-estimated by regressing the choice outcome on the attribute,
#' with cluster-robust standard errors clustered on respondent. Holding out no levels gives the AMCE; holding out
#' all other levels gives an estimate that only uses direct comparisons of `level_a` and `level_b`
#' (whose estimand equals the centered AFCP when `attribute` is independently randomized).
#'
#' @details
#' The estimates come from a regression of the choice outcome on `attribute` alone; the `design` and `weights`
#' passed to `cjoint::amce()` are not stored in `cjointobj` and so are not used. This is a valid estimator of the
#' AMCE when `attribute` is randomized independently of the other attributes and from the same distribution in
#' both profiles of a task (its levels need not be equally likely). It does not account for randomization
#' restrictions involving `attribute`. Because the other attributes in the `amce()` formula are left out, the
#' AMCE estimate will generally differ in finite samples from the one reported by `cjoint::amce()`, though both
#' are unbiased under independent randomization. Survey weights can be supplied through `weights`.
#'
#' A held-out set is not estimable if, after dropping its tasks, `level_a` or `level_b` no longer appears
#' against a different level in any task (for example, when every task containing `level_b` also contains a
#' held-out level, or when `level_a` and `level_b` never appear in the same task). Those rows are returned with
#' `NA` estimates, with a warning.
#'
#' @importFrom stats lm coef pnorm qnorm
#'
#' @param cjointobj Object of class `amce` returned by the `amce()` function in the `cjoint` package
#' @param respondent.id Character denoting the column identifying the unique respondent in the dataset from `cjointobj`
#' @param task.id Character denoting the column identifying the task/question in the dataset from `cjointobj`
#' @param attribute Character denoting the attribute of interest
#' @param level_a Character denoting the level whose AMCE is estimated
#' @param level_b Character denoting the reference level for the AMCE
#' @param ci Numeric between 0 and 1 denoting the size of the confidence interval to report for each estimate
#' @param vcov_type Character passed to the `type` argument of `sandwich::vcovCL()` for the cluster-robust (by respondent) variance. The default `"HC2"` is the CR2 (Bell-McCaffrey) estimator, as in `estimatr::lm_robust()`; `"HC1"` matches the standard errors from `cjoint::amce()`.
#' @param weights Character denoting a column of non-negative weights in the dataset from `cjointobj`, or `NULL` (the default) for no weights. This is not taken from `cjointobj`, so it must be supplied even if weights were passed to `cjoint::amce()`. With weights, the `"HC2"` adjustment treats them as inverse-variance weights (as `clubSandwich::vcovCR(..., inverse_var = TRUE)` does), so it can differ slightly from `estimatr::lm_robust()`, which does not.
#'
#' @return A dataframe of class `spec_curve` sorted by the estimate, with columns
#' - `held_out` - The held-out levels, separated by "; " (empty for the AMCE)
#' - `n_held_out` - The number of held-out levels
#' - `type` - "AMCE" (no levels held out), "AFCP" (all other levels held out) or "Held-out"
#' - `estimate`, `se`, `pval`, `conf_low`, `conf_high`, `conf_level` - The AMCE estimate and inference
#' - `rank` - Position of the estimate in the sorted curve (`NA` for held-out sets that are not estimable)
#' @export
#'
spec_curve <- function(cjointobj, respondent.id, task.id, attribute, level_a, level_b, ci = .95, vcov_type = "HC2", weights = NULL){

  # Shared input checks and cleaning (see utils.R) - tasks with no or multiple choices are allowed here
  setup <- afcp_setup(cjointobj, respondent.id, task.id, profile.id = NULL, attribute = attribute,
                      baseline = level_b, ci = ci, require_single_choice = FALSE)
  level_b <- setup$baseline
  level_a <- cjoint:::clean.names(level_a)

  if (!(level_a %in% setup$attr_use)){
    stop(paste("Error: 'level_a', ", level_a,  ", not in levels of 'attribute' or equal to 'level_b'", sep=""))
  }

  data <- setup$data
  y <- data[[setup$choice.outcome]]
  attr_values <- data[[setup$attribute]]
  respondent <- data[[setup$respondent.id]]
  w <- spec_curve_weights(data, weights)
  task_key <- task_keys(respondent, data[[setup$task.id]])

  # Whether each task has level_a (level_b) against some different level - needed for an estimate to be defined
  vs_other <- function(l){
    is_l <- !is.na(attr_values) & attr_values == l
    n_l <- stats::ave(as.numeric(is_l), task_key, FUN = sum)
    n_task <- stats::ave(rep(1, length(task_key)), task_key, FUN = sum)
    n_l > 0 & n_l < n_task
  }
  a_vs_other <- vs_other(level_a)
  b_vs_other <- vs_other(level_b)

  # Indicator for whether each task contains each of the other levels (in any profile)
  others <- setdiff(setup$attr_levels, c(level_a, level_b))
  task_has <- sapply(others, function(l) task_key %in% task_key[attr_values == l])
  task_has <- matrix(task_has, ncol = length(others))

  # Every subset of held-out levels, from none (the AMCE) to all (the AFCP)
  held_out_sets <- unlist(lapply(0:length(others), function(k) utils::combn(length(others), k, simplify = FALSE)),
                          recursive = FALSE)

  z <- abs(qnorm((1-ci)/2))
  results <- lapply(held_out_sets, function(drop){
    keep <- if (length(drop) == 0) rep(TRUE, length(y)) else rowSums(task_has[, drop, drop = FALSE]) == 0
    held_out <- paste(others[drop], collapse = "; ")
    if (!any(keep & a_vs_other) || !any(keep & b_vs_other)){
      return(data.frame(held_out = held_out, n_held_out = length(drop), estimate = NA_real_, se = NA_real_,
                        pval = NA_real_, conf_low = NA_real_, conf_high = NA_real_, conf_level = ci))
    }
    remaining <- setdiff(setup$attr_levels, c(level_a, level_b, others[drop]))
    treat <- factor(attr_values[keep], levels = c(level_b, level_a, remaining))
    w_keep <- w[keep]
    fit <- lm(y[keep] ~ treat, weights = w_keep)
    vcov_cl <- sandwich::vcovCL(fit, cluster = as.character(respondent[keep]), type = vcov_type) # as.character: vcovCL's HC2 is wrong when the cluster is a factor with unused levels
    coef_a <- paste0("treat", level_a) # by name: lm() drops empty levels, so position 2 isn't guaranteed to be level_a
    est <- unname(coef(fit)[coef_a])
    se <- sqrt(vcov_cl[coef_a, coef_a])
    data.frame(held_out = held_out, n_held_out = length(drop),
               estimate = est, se = se, pval = 2*pnorm(-abs(est/se)),
               conf_low = est - z*se, conf_high = est + z*se, conf_level = ci)
  })
  out <- bind_rows(results)

  not_estimable <- is.na(out$estimate)
  if (any(not_estimable)){
    sets <- ifelse(out$held_out[not_estimable] == "", "(none)", out$held_out[not_estimable])
    warning(paste("After dropping tasks, '", level_a, "' or '", level_b, "' no longer appears against a different level, ",
                  "so the estimate is NA for held-out sets: ", paste(sets, collapse = ", "), sep = ""), call. = FALSE)
  }

  out$type <- "Held-out"
  out$type[out$n_held_out == 0] <- "AMCE"
  out$type[out$n_held_out == length(others)] <- "AFCP"

  out <- out[order(out$estimate), c("held_out", "n_held_out", "type", "estimate", "se", "pval",
                                    "conf_low", "conf_high", "conf_level")]
  out$rank <- NA_integer_
  out$rank[!is.na(out$estimate)] <- seq_len(sum(!is.na(out$estimate)))
  rownames(out) <- NULL

  attr(out, "attribute") <- setup$attribute
  attr(out, "level_a") <- level_a
  attr(out, "level_b") <- level_b
  class(out) <- c("spec_curve", "data.frame")
  out
}


#' Validate and extract the weights for spec_curve()
#' @noRd
spec_curve_weights <- function(data, weights){
  if (is.null(weights)) return(rep(1, nrow(data)))
  weights <- cjoint:::clean.names(weights)
  if (!(weights %in% names(data))){
    stop("Error: 'weights' not a column in the dataset from 'cjointobj'")
  }
  w <- data[[weights]]
  if (!is.numeric(w) || anyNA(w) || any(w < 0)){
    stop("Error: 'weights' must be numeric, non-negative and not missing")
  }
  w
}


#' Plot a specification curve of held-out AMCE estimates
#'
#' Plots the estimates from [spec_curve()] in increasing order with their confidence intervals, highlighting the
#' AMCE (no levels held out) and the AFCP (all other levels held out).
#'
#' @importFrom ggplot2 .data
#'
#' @param x Object of class `spec_curve` returned by [spec_curve()]
#' @param colors Named character vector with the highlight colors for the `AMCE` and/or `AFCP` estimates, e.g.
#'   `c(AFCP = "darkgreen")`. Any element not supplied keeps its default (`AMCE = "blue"`, `AFCP = "red"`).
#' @param ylab Label for the y-axis
#' @param ... Unused
#'
#' @return A `ggplot` object
#' @export
#'
plot.spec_curve <- function(x, colors = c(AMCE = "blue", AFCP = "red"), ylab = "AMCE", ...){

  # Fill in any color not supplied with its default
  if (is.null(names(colors)) || !all(names(colors) %in% c("AMCE", "AFCP"))){
    stop("Error: `colors` must be a named vector with elements 'AMCE' and/or 'AFCP'")
  }
  colors <- c(colors, c(AMCE = "blue", AFCP = "red")[setdiff(c("AMCE", "AFCP"), names(colors))])

  x <- as.data.frame(x)
  x <- x[!is.na(x$estimate), ] # held-out sets that are not estimable
  amce <- x[x$type == "AMCE", ]
  afcp <- x[x$type == "AFCP", ]

  ggplot2::ggplot(x, ggplot2::aes(x = .data$rank, y = .data$estimate)) +
    ggplot2::geom_point(size = 1) +
    ggplot2::geom_errorbar(ggplot2::aes(ymin = .data$conf_low, ymax = .data$conf_high),
                           width = 0, linewidth = 0.2, color = "gray25") +
    ggplot2::geom_point(data = amce, color = colors[["AMCE"]], size = 3) +
    ggplot2::geom_errorbar(ggplot2::aes(ymin = .data$conf_low, ymax = .data$conf_high), data = amce,
                           width = 0, linewidth = 0.6, color = colors[["AMCE"]]) +
    ggplot2::geom_point(data = afcp, color = colors[["AFCP"]], size = 3) +
    ggplot2::geom_errorbar(ggplot2::aes(ymin = .data$conf_low, ymax = .data$conf_high), data = afcp,
                           width = 0, linewidth = 0.6, color = colors[["AFCP"]]) +
    ggplot2::geom_hline(yintercept = 0) +
    ggplot2::labs(x = "", y = ylab) +
    ggplot2::theme_minimal() +
    ggplot2::theme(axis.text.x = ggplot2::element_blank(),
                   axis.ticks.x = ggplot2::element_blank())
}
