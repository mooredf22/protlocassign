#' Protein-level test data set, NSA profiles
#' 
#' NSA profiles.  NSA refers to data in the form of normalized 
#'   specific amounts, where a constant amount of total protein from each 
#'   fraction is analyzed and the signal attributed to a specific protein 
#'   in each fraction is scaled so the values sum to one (see Tutorial 3).  
#'   The test data set was derived from TMT MS2 data from experiment 
#'   AT5 of Tannous et al and consists of mean profiles for proteins
#'   TLN1, TLN2, AIF1, and 
#'   the proteins and listed in “markerListJadot”. 
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
#' @references Tannous, A.; Boonen, M.; Zheng, H.; Zhao, C.;
#'  Germain, C. J.; Moore, D. F.;
#'  Sleat, D. E.; Jadot, M.; Lobel, P., Comparative Analysis of
#'  Quantitative Mass Spectrometric Methods for Subcellular Proteomics.
#'  J Proteome Res 2020, 19, (4), 1718-1730.
#'
#' @usage data(protNSA_test)
#' @examples
#' data(protNSA_test)
#' str(protNSA_test)
"protNSA_test"
