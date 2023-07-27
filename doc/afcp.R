## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----setup--------------------------------------------------------------------
library(cjoint)
library(afcp)

## -----------------------------------------------------------------------------
data("immigrationconjoint")

## -----------------------------------------------------------------------------
table(immigrationconjoint$`Language Skills`)

## -----------------------------------------------------------------------------
amce_results <- amce(Chosen_Immigrant ~ `Language Skills`, data=immigrationconjoint, cluster=T, respondent.id = "CaseID", design = "uniform")
summary(amce_results)

## -----------------------------------------------------------------------------
afcp_results <- afcp(amce_results, respondent.id = "CaseID", task.id = "contest_no", profile.id = "profile", attribute = "Language Skills")

## -----------------------------------------------------------------------------
afcp_results$afcp

## -----------------------------------------------------------------------------
afcp_results$wald

## -----------------------------------------------------------------------------
afcp_results$direct_indirect

