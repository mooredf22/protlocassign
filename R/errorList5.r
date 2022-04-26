#' pre-computed errors from CPA on five-fraction data
#' 
#' A set of pre-computed errors from applying CPA from five-fraction 
#'  simulated mixture data
#'  to simulated proteins
#'  resident in varying degrees in two compartments, for all
#'  combinations of pairs of compartments
#'  This is for input into the "mixtureHeatMap" program,
#'  optional argument "errorListIn". See tutorial 5 for details
#'  
#' @docType data
#' @format A data frame of pairwise errors
#'
#' \describe{
#'   \item{RSA}{errors using RSA values}
#'   \item{NSA}{errors using NSA values}
#'   \item{Acup}{errors using Acup values}
#'   \item{log RSA}{errors using log RSA values}
#'   \item{log NSA}{errors using log NSA values}
#'   \item{log Acup}{errors using log Acup values}
#'   }
#' @keywords datasets
#' @usage data(errorList5)
#' @examples
#' data(errorList5)
#' head(errorList5)
"errorList5"