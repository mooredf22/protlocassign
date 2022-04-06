#' Test data set with spectral level profiles formatted for protlocassign
#' 
#' Test data set that is the output of the function proteinDataPrep 
#'   applied to QuantPSM_test.  Like QuantPSM_test it contains descriptive 
#'   information and NSA profiles for all peptide-spectrum matches 
#'   associated with all spectra of TLN1, TLN2, AIF1 and the proteins 
#'   listed in "markerListJadot".' In addition, peptides names 
#'   (amino acid sequences) are prepended by their associated protein 
#'   name and there are unique protein and peptide numerical identifiers.  
#'   Data from the TMT MS2 data, AT5 experiment; see Tannous et
#'
#' @docType data
#' @format A test data frame
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
#'   \item{\code{protId}}{a numeric vector}
#'   \item{\code{pepId}}{a numeric vector}
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
