#' AMCE estimates for every combination of held-out attribute levels ("specification curve")
#'
#' Implements the sensitivity analysis for the inclusion of indirect comparisons. For an attribute with
#' \eqn{D} levels, there are \eqn{2^{D-2}} ways to choose which of the levels other than `level_a` and `level_b`
#' to hold out. For each combination, every task (not just every profile) containing a held-out level is dropped
#' and the AMCE of `level_a` relative to `level_b` is re-estimated by regressing the choice outcome on the attribute,
#' with cluster-robust standard errors clustered on respondent. Holding out no levels gives the AMCE; holding out
#' all other levels gives an estimate that only uses direct comparisons of `level_a` and `level_b`
#' (whose estimand equals the centered AFCP under complete randomization).
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
#'
#' @return A dataframe of class `spec_curve` sorted by the estimate, with columns
#' - `held_out` - The held-out levels, separated by "; " (empty for the AMCE)
#' - `n_held_out` - The number of held-out levels
#' - `type` - "AMCE" (no levels held out), "AFCP" (all other levels held out) or "Held-out"
#' - `estimate`, `se`, `pval`, `conf_low`, `conf_high`, `conf_level` - The AMCE estimate and inference
#' - `rank` - Position of the estimate in the sorted curve
#' @export
#'
spec_curve <- function(cjointobj, respondent.id, task.id, attribute, level_a, level_b, ci = .95, vcov_type = "HC2"){

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
  task_key <- paste(respondent, data[[setup$task.id]], sep = "_")

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
    remaining <- setdiff(setup$attr_levels, c(level_a, level_b, others[drop]))
    treat <- factor(attr_values[keep], levels = c(level_b, level_a, remaining))
    fit <- lm(y[keep] ~ treat)
    vcov_cl <- sandwich::vcovCL(fit, cluster = as.character(respondent[keep]), type = vcov_type) # as.character: vcovCL's HC2 is wrong when the cluster is a factor with unused levels
    est <- unname(coef(fit)[2])
    se <- sqrt(vcov_cl[2, 2])
    data.frame(held_out = paste(others[drop], collapse = "; "), n_held_out = length(drop),
               estimate = est, se = se, pval = 2*pnorm(-abs(est/se)),
               conf_low = est - z*se, conf_high = est + z*se, conf_level = ci)
  })
  out <- bind_rows(results)

  out$type <- "Held-out"
  out$type[out$n_held_out == 0] <- "AMCE"
  out$type[out$n_held_out == length(others)] <- "AFCP"

  out <- out[order(out$estimate), c("held_out", "n_held_out", "type", "estimate", "se", "pval",
                                    "conf_low", "conf_high", "conf_level")]
  out$rank <- seq_len(nrow(out))
  rownames(out) <- NULL

  attr(out, "attribute") <- setup$attribute
  attr(out, "level_a") <- level_a
  attr(out, "level_b") <- level_b
  class(out) <- c("spec_curve", "data.frame")
  out
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
