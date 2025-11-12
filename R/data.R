#' Riesby Data Set
#' ...
#' @format A data frame with X rows and 6 variables:
#' \describe{
#'   \item{id}{Subject ID}
#'   \item{hamd}{Hamilton depression score}
#'   \item{intcpt}{Intercept term}
#'   \item{week}{Week of study}
#'   \item{endog}{Endogenous depression status}
#'   \item{endweek}{Interaction of endog and week}
#' }
"riesby"

#' Positive Mood Data Set
#'
#' A dataset about positive mood... (describe it briefly).
#'
#' @format A data frame with X rows and Y variables:
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
"posmood"