#' Acup reference profiles
#' 
#' Profiles for the subcellular compartments based on Acup profile means of the 
#'    individual reference proteins that represent each compartment. 
#'    This function produces a matrix that has one row for each 
#'    compartment and one column for each fraction that comprises the profile.
#'   Data from the TMT MS2 data, AT5 experiment; see Tannous et al.
#'   Created using AcupFromNSA on refLocProfNSA.
#'   See help file for `refLocProfNSA` and `AcupFromNSA` for more information.
#' 
#'
#' @docType data
#' @format list of profiles
#'
#' \describe{
#'   \item{\code{N}}{a numeric vector}
#'   \item{\code{M}}{a numeric vector}
#'   \item{\code{L1}}{a numeric vector}
#'  \item{\code{L2}}{a numeric vector}
#'   \item{\code{P}}{a numeric vector}
#'   \item{\code{S}}{a numeric vector}
#'   \item{\code{Nyc1}}{a numeric vector}
#'   \item{\code{Nyc2}}{a numeric vector}
#'   \item{\code{Nyc3}}{a numeric vector}
#' }
#' @keywords datasets
#'
#' @references Tannous, A.; Boonen, M.; Zheng, H.; Zhao, C.;
#'  Germain, C. J.; Moore, D. F.;
#'  Sleat, D. E.; Jadot, M.; Lobel, P., Comparative Analysis of
#'  Quantitative Mass Spectrometric Methods for Subcellular Proteomics.
#'  J Proteome Res 2020, 19, (4), 1718-1730.
#'
#' @usage data(refLocProfAcup)
#' @examples
#' data(refLocProfAcup)
#' round(refLocProfAcup, digits=4)

"refLocProfAcup"
