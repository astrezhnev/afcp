## cjoint - load temp data
library(cjoint)
library(estimatr)

# Immigration Choice Conjoint Experiment Data from Hainmueller et. al. (2014).
data("immigrationconjoint")
data("immigrationdesign")

# Run AMCE estimator using all attributes in the design
results <- amce(ChosenImmigrant ~  `Prior Entry`, data=immigrationconjoint,
                cluster=TRUE, respondent.id="CaseID", design="uniform")
# Print summary
summary(results)

afcp_est <- afcp(results, respondent.id="CaseID", task.id = "contestno", profile.id="profile", attribute = "PriorEntry")

#' Estimates AFCPs from a cjoint object
#'
afcp <- function(cjointobj, respondent.id, task.id, profile.id, attribute, baseline = NULL){

  # Sanity checks

  # Is cjointobj an amce object?
  if (class(cjointobj) != "amce"){
    stop("Error: 'cjointobj' not of class 'amce'")
  }

  # Is attribute in cjointobj
  if (!(attribute %in% names(cjointobj$attributes))){
    stop("Error: 'attribute' not in list of attributes in 'cjointobj'")
  }

  # Get dataset
  data <- cjointobj$data

  # Get formula and outcome
  choice.outcome =  cjoint:::clean.names(as.character(cjointobj$formula)[2])

  # Get the list of levels
  attr_levels <- cjointobj$attributes[[attribute]]

  # Get the baseline
  # If null, choose existing baseline
  if (is.null(baseline)){
    baseline <- cjointobj$baselines[[attribute]]
  }

  # Error check - is baseline in the list of levels
  if (!(baseline %in% attr_levels)){
    stop(paste("Error: 'baseline', ", baseline,  ", not in levels of 'attribute'", sep=""))
  }

  # Drop baseline among usable levels
  attr_use <- attr_levels[attr_levels != baseline]

  # Clean the data columns using cjoint's clean.names() function
  data[[attribute]] <- cjoint:::clean.names(as.character(data[[attribute]]))

  # For each non-baseline attribute, estimate the AFCP relative to the baseline.
  afcp_list <- list()
  for (level in attr_use){
    # Make the data wide
    wide_data <- make.wide.data(indata=data, amce_var = attribute, level_a = level, level_b = baseline,
                                respondentID = respondent.id, choice = choice.outcome, qid = task.id, option = profile.id)

    # Get estimates
    estimates <- lm_robust(choose ~ treatment, data=wide_data, cluster=respid) # respid = respondent ID from wide_data

    # Variance-covariance matrix
    var_cov <- vcov(estimates)

    # Get AFCPs from model
    afcp_est = c(estimates$coefficients[1])

    # Get standard errors
    afcp_se = sqrt(c(var_cov[1,1]))

    # Generate results matrix
    out_results <- data.frame(name = paste(level, baseline, sep = ", "), afcp = afcp_est, afcp_centered = afcp_est - .5, se = afcp_se)
    out_results$afcp_zstat <- (out_results$afcp - .5)/out_results$se
    out_results$afcp_pval <- 2*pnorm(-abs(out_results$afcp_zstat))

    rownames(out_results) <- NULL

    # Wald test for AFCP-transitivity
    # AFCP(a,b) - 1/2 = AFCP(a,c) - AFCP(b,c)
    # \beta_0 - \beta_1 + \beta_2 = 1/2

    # Number of *other* levels
    L_other <- (length(estimates$coefficients) - 1)/2 #

    # If there's more than two levels
    if (L_other > 0){
      # Constraint Mat
      CMat <- matrix(nrow=L_other, ncol=length(estimates$coefficients)) # L-2 constraints
      for (k in 1:L_other){
        CMat[k,] <- 0
        CMat[k,1] <- 1
        CMat[k,1 + k] <- -1
        CMat[k,1 + L_other + k] <- 1 # This works because of how we've arranged the levels - First level_a - level_c then level_b - level_c (same order of the Cs)
      }
      colnames(CMat) <- names(estimates$coefficients)
      # Const equality
      c <- rep(1/2, L_other)

      # Flag problematic divergences
      direct_v_indirect = c(CMat%*%estimates$coefficients - c)
      names(direct_v_indirect) =  gsub(paste0("treatment",level,", "), "", names(estimates$coefficients)[2:(L_other+1)])

      # Construct the wald test statistic
      wald_stat <- t(CMat%*%estimates$coefficients - c)%*%solve((CMat)%*%var_cov%*%t((CMat)))%*%(CMat%*%estimates$coefficients - c)

      # Get a p-value
      wald_p <- pchisq(wald_stat, L_other, lower.tail=F)
    }else{
      direct_v_indirect <- NA
      wald_stat <- NA
      wald_p <- NA
    }

    afcp_list[[level]] <- list(model_afcp = estimates, summary = out_results, L_other = L_other, wald_stat = wald_stat, wald_p = wald_p, direct_indirect = direct_v_indirect, CMat = CMat)
  }

  return(afcp_list)
}


#' Make the dataset wide from the default long to estimate AFCPs - drops all observations w/o level_a or level_b
#' @param indata Dataframe containing the "long" dataset from `cjoint()`
#' @param amce_var Character name of the attribute to split on
#' @param level_a Character denoting the name of the first attribute level of interest
#' @param level_b Character denoting the name of the second attribute level of interest
#' @param respondentID Character denoting the column identifying the unique respondent
#' @param choice Character denoting the column identifying the binary choice outcome variable
#' @param qid Character denoting the question identifier
#' @param option Character denoting the identifier for each option/choice

#' @return A dataframe with the following columns

make.wide.data <- function(indata, amce_var, level_a, level_b, respondentID, choice, qid, option){

  # Get attribute, question, respondent, and answers
  sub <- indata[,c(amce_var, respondentID, choice, qid, option)]

  # Sanity checks
  # There must only be two options
  if (!identical(unique(sub[[option]]), c(1,2))){
    stop("Error: Conjoint contains tasks with more than two profiles")
  }

  # Rearrange
  sub <- sub %>% arrange(respondentID, qid, option)

  # Split
  sub1 <- sub[sub[[option]]==1,]
  sub2 <- sub[sub[[option]]==2,]
  names(sub1) <- c("val1", "respid", "choose1", "qid", "option")
  names(sub2) <- c("val2", "respid", "choose2", "qid", "option")

  # Merge sub2 to sub1
  sub_merge <- merge(sub1, sub2, by = c("respid", "qid"))

  # Restructure the data - drop if not contain either level_a or level_b
  sub_merge <- sub_merge %>% filter(val1 %in% c(level_a, level_b) | val2 %in% c(level_a, level_b))

  # Drop duplicated tasks
  sub_merge <- sub_merge %>% filter(val1 != val2)

  # Get all of the possible levels of the amce var ordered alphabetically
  level_list = sort(unique(c(sub_merge$val1, sub_merge$val2)))

  # Flip so that we get no AFCPs that are duplicated
  sub_final <- sub_merge

  # Flip so level_a is first and level_b is second
  sub_final[as.character(sub_merge$val2) == level_a & as.character(sub_merge$val1) == level_b,]$val1 <- sub_merge[as.character(sub_merge$val2) == level_a & as.character(sub_merge$val1) == level_b,]$val2
  sub_final[as.character(sub_merge$val2) == level_a & as.character(sub_merge$val1) == level_b,]$val2 <- sub_merge[as.character(sub_merge$val2) == level_a & as.character(sub_merge$val1) == level_b,]$val1
  sub_final[as.character(sub_merge$val2) == level_a & as.character(sub_merge$val1) == level_b,]$choose1 <- sub_merge[as.character(sub_merge$val2) == level_a & as.character(sub_merge$val1) == level_b,]$choose2
  sub_final[as.character(sub_merge$val2) == level_a & as.character(sub_merge$val1) == level_b,]$choose2 <- sub_merge[as.character(sub_merge$val2) == level_a & as.character(sub_merge$val1) == level_b,]$choose1

  # Flip so any level_a is first
  sub_final[as.character(sub_merge$val2) == level_a & as.character(sub_merge$val1) != level_b,]$val1 <- sub_merge[as.character(sub_merge$val2) == level_a & as.character(sub_merge$val1) != level_b,]$val2
  sub_final[as.character(sub_merge$val2) == level_a & as.character(sub_merge$val1) != level_b,]$val2 <- sub_merge[as.character(sub_merge$val2) == level_a & as.character(sub_merge$val1) != level_b,]$val1
  sub_final[as.character(sub_merge$val2) == level_a & as.character(sub_merge$val1) != level_b,]$choose1 <- sub_merge[as.character(sub_merge$val2) == level_a & as.character(sub_merge$val1) != level_b,]$choose2
  sub_final[as.character(sub_merge$val2) == level_a & as.character(sub_merge$val1) != level_b,]$choose2 <- sub_merge[as.character(sub_merge$val2) == level_a & as.character(sub_merge$val1) != level_b,]$choose1

  # Flip so any level_b is first if level_a is not first
  sub_final[as.character(sub_merge$val2) == level_b & as.character(sub_merge$val1) != level_a,]$val1 <- sub_merge[as.character(sub_merge$val2) == level_b & as.character(sub_merge$val1) != level_a,]$val2
  sub_final[as.character(sub_merge$val2) == level_b & as.character(sub_merge$val1) != level_a,]$val2 <- sub_merge[as.character(sub_merge$val2) == level_b & as.character(sub_merge$val1) != level_a,]$val1
  sub_final[as.character(sub_merge$val2) == level_b & as.character(sub_merge$val1) != level_a,]$choose1 <- sub_merge[as.character(sub_merge$val2) == level_b & as.character(sub_merge$val1) != level_a,]$choose2
  sub_final[as.character(sub_merge$val2) == level_b & as.character(sub_merge$val1) != level_a,]$choose2 <- sub_merge[as.character(sub_merge$val2) == level_b & as.character(sub_merge$val1) != level_a,]$choose1

  # Make a joint treatment
  sub_final$treatment <- paste(sub_final$val1, sub_final$val2, sep=", ")

  # what are the treatment conditions?
  third_treatment_conditions <- unique(sub_final$treatment)
  third_treatment_conditions <- third_treatment_conditions[third_treatment_conditions != paste(level_a, level_b, sep=", ")]
  third_treatment_conditions <- third_treatment_conditions[order(third_treatment_conditions)]

  # Reorder the joint treatment so A v. B is first, then all A third-comparisons and all B third-comparisons
  level_order <- c(paste(level_a, level_b, sep=", "),
                   third_treatment_conditions[grepl(level_a, third_treatment_conditions)],
                   third_treatment_conditions[grepl(level_b, third_treatment_conditions)])

  # Force into a factor
  sub_final$treatment <- factor(sub_final$treatment, levels = level_order)

  # Return it
  return(sub_final %>% select(respid, qid, choose = choose1, treatment, val1, val2))

}

#' Estimate FCPs
#' amce_var: Name of variable
#' cluster: Cluster variable (respondent ID)
#' level_a: First level of variable
#' level_b: Second level of variable
#' level_c: Alternate FCP
fcp.est.3 <- function(indata, amce_var, level_a, level_b, level_c){

  # Make the data wide
  data_wide <- make.wide.data(indata, amce_var, level_a, level_b)

  # Subset out level_c
  data_wide <- data_wide %>% filter((val1 == level_a&val2 == level_b)|val2 == level_c)

  # Get estimates
  estimates <- lm_robust(choose ~ treatment, data=data_wide, cluster=respid) # hardcoding respid cluster for now

  afcp_names = c(paste(level_a, level_b, sep=", "), levels(data_wide$treatment)[grepl(level_c, levels(data_wide$treatment))])

  # Variance-covariance matrix
  var_cov <- vcov(estimates)

  # Get AFCPs from model
  afcp_est = c(estimates$coefficients[1], estimates$coefficients[2] + estimates$coefficients[1], estimates$coefficients[3] + estimates$coefficients[1])

  # Get standard errors
  afcp_se = sqrt(c(var_cov[1,1], var_cov[1,1] + var_cov[2,2] + 2*var_cov[1,2], var_cov[1,1] + var_cov[3,3] + 2*var_cov[1,3]))

  # Generate results matrix
  out_results <- data.frame(name = afcp_names, afcp = afcp_est, se = afcp_se)
  out_results$zstat <- (out_results$afcp - .5)/out_results$se
  out_results$pval <- 2*pnorm(-abs(out_results$zstat))

  rownames(out_results) <- NULL

  # Wald test for AFCP-transitivity
  # AFCP(a,b) - 1/2 = AFCP(a,c) - AFCP(b,c)
  # \beta_0 - \beta_1 + \beta_2 = 1/2

  # Constraint Mat
  CMat <- matrix(nrow=1, ncol=3)
  CMat[1,] <- c(1, -1, 1)
  # Const equality
  c <- 1/2

  # Construct the wald test statistic
  # Construct the wald test statistic
  wald_stat <- t(CMat%*%estimates$coefficients - c)%*%solve(CMat%*%var_cov%*%t(CMat))%*%(CMat%*%estimates$coefficients - c)

  # Get a p-value
  wald_p <- pchisq(wald_stat, 1, lower.tail = F)

  return(list(model = estimates, summary = out_results, wald_stat = wald_stat, wald_p = wald_p))

}
