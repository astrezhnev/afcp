test_that("spec curve has one estimate per subset of held-out levels with AMCE and AFCP endpoints", {
  dat <- sim_conjoint(N = 300, K = 5, lvls = c("A", "B", "C", "D", "E"),
                      afcps = c(AB = .6, AC = .4, BC = .6))
  cj <- fit_cjoint(dat)
  sc <- spec_curve(cj, "unit", "task", "attribute", level_a = "A", level_b = "B")

  expect_s3_class(sc, "spec_curve")
  expect_equal(nrow(sc), 2^3)
  expect_equal(sum(sc$type == "AMCE"), 1)
  expect_equal(sum(sc$type == "AFCP"), 1)
  expect_equal(sc$rank, seq_len(8))
  expect_false(is.unsorted(sc$estimate))

  # No held-out levels: the bivariate AMCE
  amce_fit <- lm(choice ~ relevel(attribute, "B"), data = dat)
  expect_equal(sc$estimate[sc$type == "AMCE"], unname(coef(amce_fit)[2]))

  # All held out: only tasks with just A and B
  task_key <- paste(dat$unit, dat$task)
  bad <- unique(task_key[!dat$attribute %in% c("A", "B")])
  sub <- dat[!task_key %in% bad, ]
  expect_equal(sc$estimate[sc$type == "AFCP"],
               unname(coef(lm(choice ~ relevel(droplevels(attribute), "B"), data = sub))[2]))
})

test_that("spec curve SEs do not depend on unused levels of a factor respondent id", {
  dat <- sim_conjoint(N = 200, K = 3, lvls = c("A", "B", "C", "D"), afcps = c(AB = .6))
  dat_f <- dat
  dat_f$unit <- factor(dat_f$unit, levels = 1:500)   # levels with no observations
  sc <- spec_curve(fit_cjoint(dat), "unit", "task", "attribute", "A", "B")
  sc_f <- spec_curve(fit_cjoint(dat_f), "unit", "task", "attribute", "A", "B")
  expect_equal(sc$se, sc_f$se)
})

test_that("plot method returns a ggplot", {
  dat <- sim_conjoint(N = 100, K = 5, lvls = c("A", "B", "C", "D"), afcps = c(AB = .6))
  sc <- spec_curve(fit_cjoint(dat), "unit", "task", "attribute", "A", "B")
  expect_s3_class(plot(sc), "ggplot")
  expect_error(plot(sc, colors = c(AMCE = "red", Other = "green")), "colors")
  expect_error(plot(sc, colors = "red"), "colors")

  # Colors not supplied keep their defaults
  point_colors <- function(p) unname(vapply(p$layers[c(3, 5)], function(l) l$aes_params$colour, character(1)))
  expect_equal(point_colors(plot(sc)), c("blue", "red"))
  expect_equal(point_colors(plot(sc, colors = c(AFCP = "darkgreen"))), c("blue", "darkgreen"))
  expect_equal(point_colors(plot(sc, colors = c(AFCP = "black", AMCE = "orange"))), c("orange", "black"))
})

test_that("spec curve checks its levels", {
  dat <- sim_conjoint(N = 100, K = 5, lvls = c("A", "B", "C"), afcps = c(AB = .6))
  cj <- fit_cjoint(dat)
  expect_error(spec_curve(cj, "unit", "task", "attribute", "Z", "B"), "level_a")
  expect_error(spec_curve(cj, "unit", "task", "attribute", "B", "B"), "level_a")
})

# Conjoint data from a given matrix of task pairs (one row per task, one column per profile)
sim_from_pairs <- function(pairs, lvls, K = 5, p = .6, seed = 1){
  set.seed(seed)
  n <- nrow(pairs)
  choice_1 <- rbinom(n, 1, p)
  data.frame(unit = rep(rep(seq_len(n / K), each = K), each = 2),
             task = rep(rep(seq_len(K), times = n / K), each = 2),
             profile = rep(c(1, 2), times = n),
             attribute = factor(c(t(pairs)), levels = lvls),
             choice = c(rbind(choice_1, 1L - choice_1)))
}
random_pairs <- function(n, lvls, seed = 1){
  set.seed(seed)
  cbind(sample(lvls, n, replace = TRUE), sample(lvls, n, replace = TRUE))
}

test_that("held-out sets where level_b disappears are NA, not a contrast against another level", {
  lvls <- c("A", "B", "C", "D")
  P <- random_pairs(4000, lvls)
  has <- function(l) rowSums(P == l) > 0
  P <- P[!has("B") | has("C"), ][1:1500, ]   # every task with B also has C
  dat <- sim_from_pairs(P, lvls)

  expect_warning(sc <- spec_curve(fit_cjoint(dat), "unit", "task", "attribute", "A", "B"), "C; D")
  expect_true(all(is.na(sc$estimate[grepl("C", sc$held_out)])))
  expect_true(all(is.na(sc$rank[grepl("C", sc$held_out)])))
  expect_false(anyNA(sc$estimate[!grepl("C", sc$held_out)]))
  expect_equal(sort(sc$rank), 1:2)
  expect_s3_class(plot(sc), "ggplot")
})

test_that("held-out sets where level_a disappears are NA", {
  lvls <- c("A", "B", "C", "D")
  P <- random_pairs(4000, lvls)
  has <- function(l) rowSums(P == l) > 0
  P <- P[!has("A") | has("D"), ][1:1500, ]   # every task with A also has D
  dat <- sim_from_pairs(P, lvls)

  expect_warning(sc <- spec_curve(fit_cjoint(dat), "unit", "task", "attribute", "A", "B"), "level")
  expect_true(all(is.na(sc$estimate[grepl("D", sc$held_out)])))
  expect_false(anyNA(sc$estimate[!grepl("D", sc$held_out)]))
})

test_that("the AFCP is NA when level_a and level_b never appear in the same task", {
  lvls <- c("A", "B", "C", "D")
  P <- random_pairs(4000, lvls)
  P <- P[!(P[, 1] %in% c("A", "B") & P[, 2] %in% c("A", "B") & P[, 1] != P[, 2]), ][1:1500, ]
  dat <- sim_from_pairs(P, lvls)

  expect_warning(sc <- spec_curve(fit_cjoint(dat), "unit", "task", "attribute", "A", "B"), "C; D")
  expect_true(is.na(sc$estimate[sc$type == "AFCP"]))
  expect_equal(sum(is.na(sc$estimate)), 1)
})

test_that("held-out tasks are matched on respondent and task, not on a pasted id", {
  dat <- sim_conjoint(N = 100, K = 5, lvls = c("A", "B", "C", "D"), afcps = c(AB = .6, AC = .7))
  # Respondent "a" task "b_1" and respondent "a_b" task "1" would both paste to "a_b_1"
  dat_s <- dat
  dat_s$unit <- c("a", "a_b", as.character(3:100))[dat$unit]
  dat_s$task <- ifelse(dat$unit == 1, paste0("b_", dat$task), as.character(dat$task))
  sc <- spec_curve(fit_cjoint(dat), "unit", "task", "attribute", "A", "B")
  sc_s <- spec_curve(fit_cjoint(dat_s), "unit", "task", "attribute", "A", "B")
  expect_equal(sc_s$estimate, sc$estimate)
  expect_equal(sc_s$se, sc$se)
})

test_that("weighted estimates and SEs match estimatr (HC1) and clubSandwich (CR2), and weights are validated", {
  skip_if_not_installed("estimatr")
  skip_if_not_installed("clubSandwich")
  dat <- sim_conjoint(N = 200, K = 5, lvls = c("A", "B", "C", "D"), afcps = c(AB = .6, AC = .7))
  dat$svy_wt <- runif(200, .5, 2)[dat$unit]
  dat$ones <- 1
  cj <- fit_cjoint(dat)

  expect_equal(spec_curve(cj, "unit", "task", "attribute", "A", "B", weights = "ones"),
               spec_curve(cj, "unit", "task", "attribute", "A", "B"))

  # HC1 is the Stata-style CR1 in lm_robust
  sc1 <- spec_curve(cj, "unit", "task", "attribute", "A", "B", vcov_type = "HC1", weights = "svy_wt")
  ref <- estimatr::lm_robust(choice ~ relevel(attribute, "B"), data = dat, weights = svy_wt,
                             clusters = unit, se_type = "stata")
  expect_equal(sc1$estimate[sc1$type == "AMCE"], unname(coef(ref)[2]))
  expect_equal(sc1$se[sc1$type == "AMCE"], unname(ref$std.error[2]))

  # sandwich's weighted HC2 is the CR2 with an inverse-variance working model; lm_robust uses an identity one
  sc2 <- spec_curve(cj, "unit", "task", "attribute", "A", "B", weights = "svy_wt")
  wfit <- lm(choice ~ relevel(attribute, "B"), data = dat, weights = svy_wt)
  vcr <- clubSandwich::vcovCR(wfit, cluster = dat$unit, type = "CR2", inverse_var = TRUE)
  expect_equal(sc2$se[sc2$type == "AMCE"], sqrt(vcr[2, 2]))

  dat$neg_wt <- -1
  expect_error(spec_curve(fit_cjoint(dat), "unit", "task", "attribute", "A", "B", weights = "neg_wt"), "weights")
  expect_error(spec_curve(cj, "unit", "task", "attribute", "A", "B", weights = "nope"), "weights")
})
