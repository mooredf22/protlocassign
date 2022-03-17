#' Test data set with all spectra of QuantPSM_test profile
#'
#' Drawn from data from the TMT MS2 data, AT5 experiment; see Tannous et al
#'
#' @docType data
#' @format a test data set of NSA spectra
#'
#' \describe{
#'   \item{\code{Protein.Gene}}{protein name}
#'   \item{\code{PeptideSequence}}{peptide}
#'   \item{\code{SpectrumId}}{identification label of spectrum}
#'   \item{\code{PeptidesStartPostionInProtein}}{number of start position}
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
