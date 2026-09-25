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
