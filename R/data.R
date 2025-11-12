#' Riesby Data Set
#' 
#' 
#' A longitudinal dataset that examined the Hamilton depression plasma scores and 
#' clinical response
#' 
#' 66 depressed inpatients (37 endogenous and 29 non-endogenous) with 396 total 
#' observations. 
#' 
#' @format A data frame with 396 rows and 6 variables:
#' \describe{
#'   \item{id}{Subject ID}
#'   \item{hamd}{Hamilton depression score}
#'   \item{intcpt}{Intercept term}
#'   \item{week}{Week of study}
#'   \item{endog}{Endogenous depression status}
#'   \item{endweek}{Interaction of endog and week}
#' }
#' @references Reisby, Niels, et al. "Imipramine: clinical effects and pharmacokinetic variability." Psychopharmacology 54.3 (1977): 263-272.
"riesby"

#' Positive Mood Data Set
#'
#' An intensive longitudinal dataset from an EMA study that examines the effect of 
#' various covariates on positive mood variation in adolescents.
#' 
#' 516 participants with a total of 17574 observations are collected.
#'
#' @format A data frame with 17574 rows and 19 variables:
#' \describe{
#'   \item{id}{Subject ID}
#'   \item{posmood}{Positive mood score}
#'   \item{negmood}{Negative mood score}
#'   \item{t1}{time-of-day indicator of 9am-1:59pm}
#'   \item{t2}{time-of-day indicator of 2pm-5:59pm}
#'   \item{t3}{time-of-day indicator of 6pm-9:59pm}
#'   \item{t4}{time-of-day indicator of 10pm-2:59am}
#'   \item{w1}{day-of-week indicator of Monday}
#'   \item{w2}{day-of-week indicator of Tuesday}
#'   \item{w3}{day-of-week indicator of Wednesday}
#'   \item{w4}{day-of-week indicator of Thursday}
#'   \item{w5}{day-of-week indicator of Friday}
#'   \item{w6}{day-of-week indicator of Saturday}
#'   \item{other_bs}{subject mean of the prompt-varying ”with others” indicator}
#'   \item{other_ws}{prompt-varying deviation of the ”with others” indicator minus the subject mean of this variable.}
#'   \item{genderf}{gender, 0 for males, 1 for females}
#'   \item{age15}{Subject's age minus 15}
#'   \item{tirbor}{prompt-level assessments of tired/bored}
#'   \item{frustr}{prompt-level assessments of frustrated/stressed}
#' }
#' 
#' @references Mermelstein, R., Hedeker, D., and Weinstein, S. (2010). Ecological momentary assessment of mood-smoking relationships in adolescent smokers. American Psychological Association
"posmood"