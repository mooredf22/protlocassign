#' Formats a data set of spectral profiles for use with protlocassign
#' 
#' The function sorts a data set of peptide-spectrum matches with
#'    associated profiles and descriptive information, 
#'    ensuring that spectra are nested within peptides and peptides 
#'    are nested within proteins so that all peptide and protein data 
#'    are contiguous. Peptide names (typically sequences) are prepended with 
#'    the associated protein name. The function assigns 
#'    unique numerical identifiers to each protein (protId) and 
#'    each peptide (pepId). The input data frame must be formatted 
#'    so that the first column contains protein names and the second 
#'    peptide names.  
#'    This is followed by a series of columns containing descriptive 
#'    information and profile data.  

#' @param protClass data frame containing protein and peptide names, 
#'        annotation data and profiles
#' @param numRefCols number of reference columns  
#'        (includes protein and peptide names)
#' @param numDataCols number of fractions in each profile 
#' @return Data frame properly formatted for protlocassign
#' @examples
#' data(QuantPSM_test)
#' specNSA_out <- proteinDataPrep(QuantPSM_test, numRefCols=5,numDataCols=9)
#' head(specNSA_out)
#' @export

proteinDataPrep <- function(protClass, numRefCols,
    numDataCols) {
    names(protClass)[seq_len(2)] <- c("prot", "peptide")
    # remove '(prot NAME)'
    protNamesUpper <- toupper(protClass$prot)
    protClass$prot <- protNamesUpper

    # strip any whitespace before or after
    # protein name
    protClass$prot <- trimws(protClass$prot)

    # ensure that there are no extra columns
    # after the data columns
    protClassOrig <- protClass
    protClass <- protClassOrig[, seq_len((numRefCols +
        numDataCols))]

    # # # # # # # # # # get unique peptides by
    # pasting the protein and peptide names
    protSeqProteinModifTemp <- paste(protClass$prot,
        protClass$peptide, sep = "::")
    # list of unique peptides
    uniquePeptideList <- unique(protSeqProteinModifTemp)
    uniqueProtList <- unique(protClass$prot)

    # order proteins, assuring that the
    # associated peptides travel with them
    uniquePeptideOrderInd <- order(protSeqProteinModifTemp)
    protClassSort <- protClass[uniquePeptideOrderInd,
        ]

    protPep <- paste(protClassSort$prot, protClassSort$peptide,
        sep = "::")

    protId <- cumsum(!duplicated(protClassSort$prot))
       # gives a unique number to each protein
    pepId <- cumsum(!duplicated(protPep))  # unique number to each peptide

    # replace plain peptide with protPep, in
    # column 2, to protect against non-unique
    # peptide names
    protClassSort[, 2] <- protPep
    protClassExtend <- data.frame(protClassSort, protId,
        pepId)
    protClassExtend
}




