#' Shared input checks and cleaning for functions that take a `cjoint` amce object
#'
#' Cleans the id/attribute names in the same way `cjoint::amce()` does, validates the
#' inputs and returns the pieces of `cjointobj` needed downstream.
#'
#' @param require_single_choice If `TRUE`, error unless every task has exactly one chosen profile.
#' @return A list with the cleaned names, the data (with the attribute column cleaned),
#'   the outcome name, the attribute levels, the baseline, and the non-baseline levels.
#' @noRd
afcp_setup <- function(cjointobj, respondent.id, task.id, profile.id = NULL, attribute,
                       baseline = NULL, ci = .95, require_single_choice = TRUE){

  # Clean respondent.id, task.id, profile.id and attribute in line with what amce() does in cjoint
  respondent.id <- cjoint:::clean.names(respondent.id)
  task.id <- cjoint:::clean.names(task.id)
  if (!is.null(profile.id)){
    profile.id <- cjoint:::clean.names(profile.id)
  }
  attribute <- cjoint:::clean.names(attribute)

  # Sanity checks
  # Is cjointobj an amce object?
  if (!inherits(cjointobj, "amce")){
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
  choice.outcome =  cjoint:::clean.names(all.vars(cjointobj$formula)[1])

  # Is the choice outcome a number
  if (!all(data[[choice.outcome]] %in% c(0,1))){
    stop("Error: Outcome is not a binary indicator (0 or 1)")
  }

  # Are the id columns in the data
  ids <- c(respondent.id = respondent.id, task.id = task.id, profile.id = profile.id)
  missing_id <- !(ids %in% names(data))
  if (any(missing_id)){
    stop(paste("Error: '", names(ids)[missing_id][1], "', ", ids[missing_id][1], ", not a column in the dataset from 'cjointobj'", sep=""))
  }

  # Does every task have exactly one choice (checked per task: a total count misses offsetting errors)
  if (require_single_choice && any(rowsum(data[[choice.outcome]], task_keys(data[[respondent.id]], data[[task.id]])) != 1)){
    stop("Error: Some tasks have neither or both options selected")
  }

  # Get the list of levels
  attr_levels <- cjointobj$attributes[[attribute]]

  # If fewer then three levels, send error
  if(length(attr_levels) <= 2){
    stop(paste("Error: Not more than 2 levels in `attribute`", sep=""))
  }

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

  list(respondent.id = respondent.id, task.id = task.id, profile.id = profile.id,
       attribute = attribute, data = data, choice.outcome = choice.outcome,
       attr_levels = attr_levels, baseline = baseline, attr_use = attr_use)
}

#' Pair indicators for the wide data: a-b, then a-c and b-c for each other level c (sorted)
#'
#' make.wide.data() orients tasks so level_a is first, then level_b; other levels c are second. Errors if any task
#' doesn't map to exactly one pair, or if any pair has no tasks. Used by afcp() and afcp_gmm().
#' @noRd
level_pairs <- function(wide_data, level_a, level_b){
  v1 <- as.character(wide_data$val1)
  v2 <- as.character(wide_data$val2)
  other_levels <- sort(setdiff(unique(c(v1, v2)), c(level_a, level_b)))
  pairs <- c(list(c(level_a, level_b)),
             lapply(other_levels, function(l) c(level_a, l)),
             lapply(other_levels, function(l) c(level_b, l)))
  ind <- sapply(pairs, function(p) v1 == p[1] & v2 == p[2])
  ind <- matrix(ind, nrow = nrow(wide_data))
  if (any(is.na(ind)) || any(rowSums(ind) != 1)){
    stop("Error: some tasks do not map to exactly one level comparison")
  }
  empty <- colSums(ind) == 0
  if (any(empty)){
    p <- pairs[[which(empty)[1]]]
    stop(paste("Error: no tasks compare ", p[1], " and ", p[2], sep=""))
  }
  list(ind = ind, other_levels = other_levels, L_other = length(other_levels))
}

#' One key per respondent-task
#'
#' Built from integer codes rather than the raw ids, so that e.g. respondent "1_2" task "3" can't match respondent
#' "1" task "2_3".
#' @noRd
task_keys <- function(respondent, task){
  paste(match(respondent, unique(respondent)), match(task, unique(task)))
}
