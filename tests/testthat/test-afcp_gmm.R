test_that("GMM recovers the AFCP and is more precise when direct and indirect preferences agree", {
  # Transitive: AFCP(A,B) - .5 = AFCP(A,C) - AFCP(B,C) and similarly for D
  dat <- sim_conjoint(N = 2000, K = 5, lvls = c("A", "B", "C", "D"),
                      afcps = c(AB = .6, AC = .7, BC = .6, AD = .65, BD = .55, CD = .45))
  cj <- fit_cjoint(dat)
  g <- afcp_gmm(cj, "unit", "task", "profile", "attribute", baseline = "B")
  a <- afcp(cj, "unit", "task", "profile", "attribute", baseline = "B")

  expect_equal(g$afcp$level, a$afcp$level)
  expect_equal(g$afcp$afcp[g$afcp$level == "A"], .6, tolerance = .03)
  expect_true(all(g$afcp$se < a$afcp$se))
  expect_equal(g$jtest$L_other, c(2, 2, 2))
  expect_true(all(g$jtest$j_p > .01))
  expect_equal(g$afcp$conf_high - g$afcp$afcp, qnorm(.975) * g$afcp$se)
})

test_that("J-test rejects when there is a preference cycle", {
  dat <- sim_conjoint(N = 2000, K = 5, lvls = c("A", "B", "C"),
                      afcps = c(AB = .6, AC = .4, BC = .6))
  cj <- fit_cjoint(dat)
  g <- afcp_gmm(cj, "unit", "task", "profile", "attribute", baseline = "B")
  expect_lt(g$jtest$j_p[g$jtest$level == "A"], .01)
})

test_that("input checks match afcp()", {
  dat <- sim_conjoint(N = 100, K = 5, lvls = c("A", "B", "C"), afcps = c(AB = .5, AC = .5, BC = .5))
  cj <- fit_cjoint(dat)
  expect_error(afcp_gmm(cj, "unit", "task", "profile", "attribute", baseline = "Z"), "not in levels")
  expect_error(afcp_gmm(cj, "unit", "task", "profile", "attribute", ci = 1.5), "ci")
})

test_that("influence-function J-test equals the joint Wald test from afcp()", {
  dat <- sim_conjoint(N = 400, K = 5, lvls = c("A", "B", "C", "D"),
                      afcps = c(AB = .6, AC = .4, BC = .6), seed = 3)
  cj <- fit_cjoint(dat)
  g <- afcp_gmm(cj, "unit", "task", "profile", "attribute", baseline = "B")
  a <- afcp(cj, "unit", "task", "profile", "attribute", baseline = "B")
  expect_equal(g$jtest$j_stat, a$wald$wald_stat_all)
  expect_equal(g$jtest$j_p, a$wald$wald_p_all)
})


test_that("GMM is the inverse-variance weighted combination of the direct and indirect estimates", {
  dat <- sim_conjoint(N = 400, K = 5, lvls = c("A", "B", "C"),
                      afcps = c(AB = .55, AC = .6, BC = .55), seed = 2)
  cj <- fit_cjoint(dat)
  g <- afcp_gmm(cj, "unit", "task", "profile", "attribute", baseline = "B")

  # One restriction: combine AFCP(A,B) with 1/2 + AFCP(A,C) - AFCP(B,C) using afcp()'s regression covariance
  wd <- make.wide.data(cj$data, "attribute", "A", "B", "unit", "choice", "task", "profile")
  fit <- lm(choose ~ treatment, data = wd)
  V <- sandwich::vcovCL(fit, cluster = as.character(wd$respid), type = "HC2")
  C <- rbind(c(1, 0, 0), c(0, 1, -1))
  m <- c(C %*% coef(fit)) + c(0, .5)
  W <- solve(C %*% V %*% t(C))
  expect_equal(g$afcp$afcp[g$afcp$level == "A"], sum(W %*% m) / sum(W))
  expect_equal(g$afcp$se[g$afcp$level == "A"], sqrt(1 / sum(W)))
})

test_that("integer or reordered profile ids are accepted", {
  dat <- sim_conjoint(N = 100, K = 5, lvls = c("A", "B", "C"), afcps = c(AB = .6))
  a <- afcp(fit_cjoint(dat), "unit", "task", "profile", "attribute", baseline = "B")
  dat_int <- dat
  dat_int$profile <- as.integer(dat_int$profile)
  dat_int <- dat_int[c(2, 1, 3:nrow(dat_int)), ]
  a_int <- afcp(fit_cjoint(dat_int), "unit", "task", "profile", "attribute", baseline = "B")
  expect_equal(a_int$afcp, a$afcp)
})
