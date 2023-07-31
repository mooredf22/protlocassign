#' Protein-level test data set, Acup profiles
#'
#' Acup profiles.  Acup denotes relative amounts, which is the amount 
#'   of given protein in a fraction / amount of that given protein 
#'   in starting material (see Tutorial 3).  
#'   The test set was derived from TMT MS2 data from experiment AT5 of 
#'   Tannous et al and consists of mean profiles for proteins 
#'   TLN1, TLN2, AIF1 and those listed in “markerListJadot”. 
#' 
#'
#' @docType data
#' @format a test data set
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
#'   \item{\code{Nspectra}}{a numeric vector}
#'   \item{\code{Npep}}{a numeric vector}
#' }
#' @keywords datasets
#' @return a test data set
#'
#' @references Tannous, A.; Boonen, M.; Zheng, H.; Zhao, C.;
#'  Germain, C. J.; Moore, D. F.;
#'  Sleat, D. E.; Jadot, M.; Lobel, P., Comparative Analysis of
#'  Quantitative Mass Spectrometric Methods for Subcellular Proteomics.
#'  J Proteome Res 2020, 19, (4), 1718-1730.
#'
#' @usage data(protAcup_test)
#' @examples
#' data(protAcup_test)
#' str(protAcup_test)
"protAcup_test"
