#' Protein-level test data set, RSA profiles
#'
#' RSA profiles.  RSA refers to data expressed as relative specific amounts, 
#'   which are ratios of two ratios: the numerator is the amount of a 
#'   given protein in a particular fraction divided by the amount of that 
#'   given protein in the starting material while the denominator is 
#'   amount of total protein in a particular fraction divided by 
#'   the amount of total protein in the starting material. 
#'   The RSA describes the fold-enrichment (RSA>1) or depletion (RSA<1) 
#'   of a protein during the fractionation process, and is analogous to 
#'   the relative specific activity term used in classical analytical 
#'   subcellular fractionation (see Tutorial 3).  
#'   The test set was derived from TMT MS2 data from experiment AT5 of 
#'   Tannous et al. and consists of mean profiles for proteins 
#'   TLN1, TLN2, AIF1 and those listed in “markerListJadot”.  
#' 
#'
#' @docType data
#' @format A test data set
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
#'
#' @return A test data set
#' @references Tannous, A.; Boonen, M.; Zheng, H.; Zhao, C.;
#'  Germain, C. J.; Moore, D. F.;
#'  Sleat, D. E.; Jadot, M.; Lobel, P., Comparative Analysis of
#'  Quantitative Mass Spectrometric Methods for Subcellular Proteomics.
#'  J Proteome Res 2020, 19, (4), 1718-1730.
#'
#' @usage data(protRSA_test)
#' @examples
#' data(protRSA_test)
#' str(protRSA_test)
"protRSA_test"
