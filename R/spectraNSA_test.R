#' Protein test data set with spectral level profiles
#' 
#' Test data set consisting of descriptive information and NSA profiles 
#'   for all peptide-spectrum matches associated with
#'   TLN1, TLN2, AIF1, NSA profiles
#'   and proteins listed in "markerListJadot".`
#'  Data from the TMT MS2 data, AT5 experiment; see Tannous et al.
#'
#' @docType data
#' @format a test data set of descriptive information and NSA 
#'  profiles for spectra
#'
#' \describe{
#'   \item{\code{prot}}{protein name}
#'   \item{\code{peptide}}{protein::peptide}
#'   \item{\code{SpectrumId}}{identification label of spectrum}
#'   \item{\code{PeptidesStartPostionInProtein}}{number of start position}
#'   \item{\code{modification}}{peptide modification list}
#'   \item{\code{N}}{a numeric vector}
#'   \item{\code{M}}{a numeric vector}
#'   \item{\code{L1}}{a numeric vector}
#'   \item{\code{L2}}{a numeric vector}
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
#' @usage data(spectraNSA_test)
#' @examples
#' data(spectraNSA_test)
#' str(spectraNSA_test)
"spectraNSA_test"
