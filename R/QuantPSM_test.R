#' Test data set containing spectral level profiles not yet 
#'    formatted for use with prolocassign package
#' 
#' Test data frame that is the precursor for spectraNSA_test and 
#'    is used as input for the function proteinDataPrep.  
#'    It contains descriptive information and NSA profiles for all 
#'    peptide-spectrum matches associated with all spectra of 
#'    TLN1, TLN2, AIF1 and the proteins listed in "markerListJadot".' 
#'    Data are from the TMT MS2 data, AT5 experiment; see Tannous et al
#'
#' @docType data
#' @format A data frame containing  NSA profiles and descriptive 
#'        information for spectra
#' @return A data frame containing NSA profiles for spectra
#'
#' \describe{
#'   \item{\code{Protein.Gene}}{protein name}
#'   \item{\code{PeptideSequence}}{peptide}
#'   \item{\code{SpectrumID}}{identification label of spectrum}
#'   \item{\code{PeptidesStartPositionInProtein}}{number of start position}
#'   \item{\code{modifications}}{peptide modification list}
#'   \item{\code{N}}{a numeric vector}
#'   \item{\code{M}}{a numeric vector}
#'   \item{\code{L1}}{a numeric vector}
#'   \item{\code{L2}}{a numeric vector}
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
#' @usage data(QuantPSM_test)
#' @examples
#' data(QuantPSM_test)
#' str(QuantPSM_test)
"QuantPSM_test"
