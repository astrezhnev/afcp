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

  # Does the number of choices equal the number of tasks
  if (require_single_choice && nrow(data)/2 != sum(data[[choice.outcome]])){
    stop("Error: Number of choices doesn't equal number of tasks. Likely some tasks have neither or both options selected")
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
