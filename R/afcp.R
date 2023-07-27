
#' Estimates AFCPs from a cjoint object
#'
#' This estimates the Average Feature Choice Probability (AFCP) for a given attribute and level combination
#' using the output from an amce object returned by the cjoint package
#'
#' @import cjoint
#' @import dplyr
#' @import sandwich
#'
#' @param cjointobj Object of class `amce` returned by the `amce()` function in the `cjoint` package
#' @param respondent.id Character denoting the column identifying the unique respondent in the dataset from `cjointobj`
#' @param task.id Character denoting the column identifying the task/question in the dataset from `cjointobj`
#' @param profile.id Character denoting the column identifying the profile in each question in the dataset from `cjointobj` - this column should take on either 1 or 2 for a binary choice task
#' @param attribute Character denoting the attribute of interest
#' @param baseline Character denoting the level of the attribute to be used as the baseline. Default of `NULL` will pull the baseline from `cjointobj`.
#' @param ci Numeric between 0 and 1 denoting the size of the confidence interval to report for each afcp
#'
#' @return A list containing three dataframes
#' - `afcp` - A dataframe containing estimated AFCPs relative to the selected baseline level.
#' - `wald` - A dataframe containing the results for the aggregate Wald test of the equivalence of direct and indirect preferences across all other levels
#' - `wald_three_level` - A dataframe containing the direct (AFCP) and indirect (difference in AFCPs) preference estimates across all other levels
#' @export
#'
afcp <- function(cjointobj, respondent.id, task.id, profile.id, attribute, baseline = NULL, ci = .95){

  # Clean respondent.id, task.id, profile.id and attribute in line with what amce() does in cjoint
  respondent.id <- cjoint:::clean.names(respondent.id)
  task.id <- cjoint:::clean.names(task.id)
  profile.id <- cjoint:::clean.names(profile.id)
  attribute <- cjoint:::clean.names(attribute)

  # Sanity checks
  # Is cjointobj an amce object?
  if (class(cjointobj) != "amce"){
    stop("Error: 'cjointobj' not of class 'amce'")
  }

  # Is attribute in cjointobj
  if (!(attribute %in% names(cjointobj$attributes))){
    stop("Error: 'attribute' not in list of attributes in 'cjointobj'")
  }

  # Is CI valid
  if (!(ci < 1&0 < ci)){
    stop("Error: `ci` invalid -- must be between 0 and 1")
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
  }else{
    baseline <- cjoint:::clean.names(baseline)
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
  afcp_results <- list() # Store the estimates + p-values
  wald_tests <- list() # Combined wald tests for cycling
  wald_three_level <- list() # Three-level tests for each pair

  for (level in attr_use){
    # Make the data wide
    wide_data <- make.wide.data(indata=data, attr_var = attribute, level_a = level, level_b = baseline,
                                respondentID = respondent.id, choice = choice.outcome, qid = task.id, option = profile.id)

    # Get estimates
    estimates <- lm(choose ~ treatment, data=wide_data) # respid = respondent ID from wide_data

    # Variance-covariance matrix
    var_cov <- sandwich::vcovCL(estimates, cluster = wide_data$respid, type = "HC2")

    # Get AFCPs from model
    afcp_est = c(estimates$coefficients[1])

    # Get standard errors
    afcp_se = sqrt(c(var_cov[1,1]))

    # Generate results matrix
    out_results <- data.frame(level = level, baseline = baseline, afcp = afcp_est, se = afcp_se)
    out_results$zstat <- (out_results$afcp - .5)/out_results$se
    out_results$pval <- 2*pnorm(-abs(out_results$zstat))


    out_results$conf_high <- afcp_est + abs(qnorm((1-ci)/2))*afcp_se
    out_results$conf_low <- afcp_est - abs(qnorm((1-ci)/2))*afcp_se
    out_results$conf_level <- ci


    rownames(out_results) <- NULL
    # Store the results for this level
    afcp_results[[level]] <- out_results

    # Wald test for AFCP-transitivity
    # AFCP(a,b) - 1/2 = AFCP(a,c) - AFCP(b,c)
    # \beta_0 - \beta_1 + \beta_2 = 1/2

    # Number of *other* levels
    L_other <- (length(estimates$coefficients) - 1)/2 #

    # If there's more than two levels
    if (L_other > 0){

      # Get names of the other levels sorted alphabetically
      other_levels <- sort(attr_use[attr_use != level])

      #####
      # For each other level, test for direct vs. indirect equivalence.
      #####

      wald_threes <- list()
      for (next_L in 1:L_other){

        # Pull the coefficients for each of the three levels
        coef_sub <- estimates$coefficients[c(1, 1 + next_L, 1 + L_other + next_L)]

        # Do the same for the variance-covariance matrix
        vcov_sub <- var_cov[c(1, 1 + next_L, 1 + L_other + next_L),][,c(1, 1 + next_L, 1 + L_other + next_L)]

        # Constraint Mat
        CMat_3 <- matrix(nrow=1, ncol=3)
        CMat_3[1,] <- c(1, -1, 1)
        # Const equality
        c <- 1/2

        # Construct the wald test statistic
        wald_stat_3 <- t(CMat_3%*%coef_sub - c)%*%solve(CMat_3%*%vcov_sub%*%t(CMat_3))%*%(CMat_3%*%coef_sub - c)

        # Get a p-value
        wald_p_3 <- pchisq(wald_stat_3, 1, lower.tail = F)

        # Get AFCPs from model
        afcp_3_est = c(coef_sub[1], coef_sub[2] + coef_sub[1], coef_sub[3] + coef_sub[1])

        # Get standard errors
        afcp_3_se = sqrt(c(vcov_sub[1,1], vcov_sub[1,1] + vcov_sub[2,2] + 2*vcov_sub[1,2], vcov_sub[1,1] + vcov_sub[3,3] + 2*vcov_sub[1,3]))


        wald_threes[[other_levels[next_L]]] <- data.frame(level = level, baseline = baseline, third = other_levels[next_L], afcp_level_baseline = afcp_3_est[1], se_level_baseline = afcp_3_se[1],
                                                                afcp_level_third = afcp_3_est[2], se_level_third = afcp_3_se[2],
                                                                afcp_baseline_third = afcp_3_est[3], se_baseline_third = afcp_3_se[3],
                                                                direct = afcp_3_est[1] - .5, indirect = afcp_3_est[2] - afcp_3_est[3],
                                                                wald_stat = wald_stat_3, wald_p = wald_p_3)

      }

      wald_three_level[[level]] <- bind_rows(wald_threes)
      rownames(wald_three_level[[level]]) <- NULL

      #######
      ### Complete test
      ######

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

      # Construct the wald test statistic
      wald_stat_all <- t(CMat%*%estimates$coefficients - c)%*%solve((CMat)%*%var_cov%*%t((CMat)))%*%(CMat%*%estimates$coefficients - c)

      # Get a p-value
      wald_p_all <- pchisq(wald_stat_all, L_other, lower.tail=F)
      wald_tests[[level]] <-  data.frame(level = level, baseline = baseline, L_other = L_other, wald_stat_all, wald_p_all)



    }
  }

  # Construct the output

  if (length(attr_use) > 1){
    afcp_estimates <- bind_rows(afcp_results)
    wald_batch <- bind_rows(wald_tests)

    return(list(attribute = attribute, baseline = baseline, afcp = afcp_estimates, wald = wald_batch, direct_indirect = bind_rows(wald_three_level)))
  }else{
    afcp_estimates <- bind_rows(afcp_results)
    return(list(afcp = afcp_estimates, wald = NULL, direct_indirect = NULL))

  }

}


#' Make the dataset wide from the default long to estimate AFCPs - drops all observations w/o level_a or level_b
#'
#' @param indata Dataframe containing the "long" dataset from `cjoint()`
#' @param attr_var Character name of the attribute to split on
#' @param level_a Character denoting the name of the first attribute level of interest
#' @param level_b Character denoting the name of the second attribute level of interest
#' @param respondentID Character denoting the column identifying the unique respondent
#' @param choice Character denoting the column identifying the binary choice outcome variable
#' @param qid Character denoting the question identifier
#' @param option Character denoting the identifier for each option/choice

#' @return A dataframe with the following columns:
#' - `respid` - Respondent identifier
#' - `qid` - Task/question identifier
#' - `choose` - An indicator for whether the profile with `val1` is selected in this task.
#' - `treatment` - Character concatenating the two levels of interest `val1` and `val2`
#' - `val1` - Level assigned to the attribute of interest for profile 1
#' - `val2` - Level assigned to the attribute of interest for profile 2
#'
#' @export

make.wide.data <- function(indata, attr_var, level_a, level_b, respondentID, choice, qid, option){

  # Get attribute, question, respondent, and answers
  sub <- indata[,c(attr_var, respondentID, choice, qid, option)]

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
