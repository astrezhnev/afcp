# Pr(level x chosen) in tasks comparing x and y, computed directly from the long data
pairwise_afcp <- function(dat, x, y){
  task_key <- paste(dat$unit, dat$task)
  has <- function(l) task_key %in% task_key[dat$attribute == l]
  mean(dat$choice[has(x) & has(y) & dat$attribute == x])
}

test_that("afcp() handles levels whose names contain another level", {
  # "1" is a substring of "10": grepl() grouping used to put "2, 10" with the "1" comparisons
  dat <- sim_conjoint(N = 300, K = 5, lvls = c("1", "2", "10"), afcps = c(`12` = .6))
  cj <- fit_cjoint(dat, baseline = "2")
  out <- afcp(cj, "unit", "task", "profile", "attribute")

  expect_equal(out$afcp$afcp[out$afcp$level == "1"], pairwise_afcp(dat, "1", "2"))
  expect_equal(out$afcp$afcp[out$afcp$level == "10"], pairwise_afcp(dat, "10", "2"))

  di <- out$direct_indirect[out$direct_indirect$level_a == "1", ]
  expect_equal(di$level_c, "10")
  expect_equal(di$afcp_ac, pairwise_afcp(dat, "1", "10"))
  expect_equal(di$afcp_bc, pairwise_afcp(dat, "2", "10"))

  expect_equal(afcp_gmm(cj, "unit", "task", "profile", "attribute")$afcp$level, c("1", "10"))
})

test_that("direct and indirect AFCPs line up with their levels when there are several other levels", {
  dat <- sim_conjoint(N = 300, K = 5, lvls = c("A", "B", "C", "D", "E"),
                      afcps = c(AB = .6, AC = .7, BC = .55, AD = .4, BD = .45, AE = .65))
  out <- afcp(fit_cjoint(dat), "unit", "task", "profile", "attribute")
  di <- out$direct_indirect[out$direct_indirect$level_a == "A", ]
  expect_equal(di$level_c, c("C", "D", "E"))
  expect_equal(di$afcp_ac, sapply(di$level_c, function(l) pairwise_afcp(dat, "A", l)), ignore_attr = TRUE)
  expect_equal(di$afcp_bc, sapply(di$level_c, function(l) pairwise_afcp(dat, "B", l)), ignore_attr = TRUE)
  expect_equal(out$wald$L_other, rep(3, 4))
})

test_that("afcp() errors when an a-c comparison has no matching b-c tasks", {
  dat <- sim_conjoint(N = 300, K = 5, lvls = c("A", "B", "C", "D"), afcps = c(AB = .6))
  task_key <- paste(dat$unit, dat$task)
  b_d <- intersect(task_key[dat$attribute == "B"], task_key[dat$attribute == "D"])
  dat <- dat[!task_key %in% b_d, ]
  expect_error(afcp(fit_cjoint(dat), "unit", "task", "profile", "attribute"), "no tasks compare B and D")
})

test_that("tasks with neither or both profiles chosen are caught even when the counts offset", {
  dat <- sim_conjoint(N = 100, K = 5, lvls = c("A", "B", "C"), afcps = c(AB = .6))
  dat$choice[1:2] <- 1L   # task 1: both chosen
  dat$choice[3:4] <- 0L   # task 2: neither chosen
  expect_equal(sum(dat$choice), nrow(dat) / 2)
  cj <- fit_cjoint(dat)
  expect_error(afcp(cj, "unit", "task", "profile", "attribute"), "neither or both")
  expect_error(afcp_gmm(cj, "unit", "task", "profile", "attribute"), "neither or both")
  expect_no_error(spec_curve(cj, "unit", "task", "attribute", "A", "B"))   # allowed in spec_curve()
})

test_that("missing id columns give a clear error", {
  cj <- fit_cjoint(sim_conjoint(N = 100, K = 5, lvls = c("A", "B", "C"), afcps = c(AB = .6)))
  expect_error(afcp(cj, "resp", "task", "profile", "attribute"), "'respondent.id', resp")
  expect_error(afcp(cj, "unit", "task", "prof", "attribute"), "'profile.id', prof")
  expect_error(spec_curve(cj, "unit", "tsk", "attribute", "A", "B"), "'task.id', tsk")
})
