test_that("vcov_type = 'HC1' reproduces cjoint's clustered standard errors", {
  dat <- sim_conjoint(N = 200, K = 4, lvls = c("A", "B", "C", "D"), afcps = c(AB = .6))
  cj <- fit_cjoint(dat)
  sc <- spec_curve(cj, "unit", "task", "attribute", "A", "B", vcov_type = "HC1")
  cj_se <- summary(cj)$amce
  expect_equal(sc$se[sc$type == "AMCE"], cj_se$`Std. Err`[cj_se$Level == "A"])
})

test_that("default vcov_type is HC2 and is passed through to sandwich::vcovCL", {
  dat <- sim_conjoint(N = 100, K = 5, lvls = c("A", "B", "C"), afcps = c(AB = .6, AC = .4, BC = .6))
  cj <- fit_cjoint(dat)
  se <- sapply(c("HC0", "HC1", "HC2"), function(s) afcp(cj, "unit", "task", "profile", "attribute", baseline = "B", vcov_type = s)$afcp$se[1])
  expect_equal(afcp(cj, "unit", "task", "profile", "attribute", baseline = "B")$afcp$se[1], unname(se["HC2"]))
  expect_true(se["HC0"] < se["HC1"])
  expect_error(afcp(cj, "unit", "task", "profile", "attribute", vcov_type = "stata"))
})

test_that("the GMM J-statistic equals the Wald statistic under every vcov_type", {
  dat <- sim_conjoint(N = 300, K = 5, lvls = c("A", "B", "C", "D"), afcps = c(AB = .6, AC = .4, BC = .6), seed = 5)
  cj <- fit_cjoint(dat)
  for (s in c("HC0", "HC1", "HC2", "HC3")){
    g <- afcp_gmm(cj, "unit", "task", "profile", "attribute", baseline = "B", vcov_type = s)
    a <- afcp(cj, "unit", "task", "profile", "attribute", baseline = "B", vcov_type = s)
    expect_equal(g$jtest$j_stat, a$wald$wald_stat_all)
  }
})
