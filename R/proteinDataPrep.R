#' The program ensures that spectra are nested within peptides, and these
#'  are nested within proteins, and all peptide and protein data are
#'  contiguous. That is, not split up into separate locations
#' This program should be run first on a comma-separeted data file
#'
#' The first column must contain prot name, and the second must be the peptide
#' The columns after 'numRefCols' contain relative abundance levels
#' There should be 'numDataCols' abundance level columns

#' @param protClass frame of protein and peptide names and
#'         data as described above
#' @param numRefCols number of reference columns
#' @param numDataCols numer of data columns
#' @return data frame like the input, but with spectra and peptides
#'        nested within proteins in same location
#' @examples
#' data(TLN1_test)
#' TLN1fixed <- proteinDataPrep(TLN1_test, numRefCols=5,numDataCols=9)
#' head(TLN1fixed)
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




