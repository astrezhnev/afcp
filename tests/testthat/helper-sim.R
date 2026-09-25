# Simulate a forced-choice conjoint with one attribute and fixed pairwise AFCPs
# afcps: named vector of Pr(row level chosen over column level), e.g. c(AB = .6, AC = .7, BC = .6, ...)
sim_conjoint <- function(N, K, lvls, afcps, seed = 1){
  set.seed(seed)
  n_tasks <- N * K
  lookup <- matrix(.5, length(lvls), length(lvls), dimnames = list(lvls, lvls))
  for (nm in names(afcps)){
    l <- strsplit(nm, "")[[1]]
    lookup[l[1], l[2]] <- afcps[[nm]]
    lookup[l[2], l[1]] <- 1 - afcps[[nm]]
  }
  attr_1 <- sample(lvls, n_tasks, replace = TRUE)
  attr_2 <- sample(lvls, n_tasks, replace = TRUE)
  choice_1 <- rbinom(n_tasks, 1, lookup[cbind(attr_1, attr_2)])
  data.frame(unit = rep(rep(1:N, each = K), each = 2),
             task = rep(rep(1:K, times = N), each = 2),
             profile = rep(c(1, 2), times = n_tasks),
             attribute = factor(c(rbind(attr_1, attr_2)), levels = lvls),
             choice = c(rbind(choice_1, 1L - choice_1)))
}

fit_cjoint <- function(dat, baseline = "B"){
  cjoint::amce(choice ~ attribute, data = dat, respondent.id = "unit",
               baselines = c(attribute = baseline))
}
