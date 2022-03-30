#' constrained CPA
#' 
#' Carry out constrained proportional assignment for protein
#'  number i in profile
#'  with constraint being the number of compartments that are allowed to
#'  vary freely; other compartment proportions are fixed at zero
#' @param profile one-row data frame of a protein name (in the row name) and
#'     relative abundance levels.
#'               Nspectra and Npep, if present, will be removed
#' @param refLocationProfiles A matrix giving the abundance level profiles
#'    of the subcellular locations
#' @param numDataCols Number of channels of abundance levels
#' @param startProps starting values for proportional assignments;
#'    set equal if this is null (default)
#' @param showProgress   default is T
#' @param maxit maximum number of iterations (default is 10000)
#' @param nCPAcomparts number of compartments to fit restricted CPA;
#'        remaining proportions are fixed at zero
#' @return data frame with CPA estimates for every combination
#'        ordered by value, with best fit listed first

#' @examples
#' data(protNSA_test)
#' data(markerListJadot)
#' nTestProts <- nrow(protNSA_test)
#' refLocationProfilesNSA <- locationProfileSetup(profile=protNSA_test,
#'     markerList=markerListJadot, numDataCols=9)
#' protCPAfromNSA_test <- fCPAsubsets(profile=protNSA_test[1,],
#'                          refLocationProfiles=refLocationProfilesNSA,
#'                          numDataCols=9, startProps=NULL, nCPAcomparts=2)
#' head(protCPAfromNSA_test)
#' @importFrom utils combn
#' @export
#'
#'
fCPAsubsets <- function(profile, refLocationProfiles,
    numDataCols, startProps = NULL, showProgress = TRUE,
    maxit = 10000, nCPAcomparts = 2) {
    # maxit and assignPRobsStart must be
    # specified assignProbsStart must be NULL or
    # have a column 'protName' and assignment
    # probabilities to use as starting values use
    # the spg function (in package BB) to assign
    # proportionate assignments to compartments
    if (nrow(profile) != 1) {
        warning("profile must have one row\n")
        return("")
    }
    # profile <- profile[1:numDataCols]
    profile <- profile[seq_len(numDataCols)]

    n.locs <- nrow(refLocationProfiles)  # number of subcellular compartments
    combins.mat <- combn(x = n.locs, m = nCPAcomparts)
    result.mat <- matrix(NA, nrow = ncol(combins.mat),
        ncol = n.locs + 1)
    locs.vec <- rep(NA, ncol(combins.mat))
    names.locs <- rownames(refLocationProfiles)
    # for (ii in 1:ncol(combins.mat)) {
    for (ii in seq_len(ncol(combins.mat))) {
        # ii=2
        result.mat[ii, ] <- fCPAone(profile, refLocationProfiles,
            numDataCols, startProps = NULL, maxit = maxit,
            ind.vary = combins.mat[, ii], minVal = TRUE)
        locs.vec[ii] <- paste(names.locs[combins.mat[,
            ii]], collapse = "")
    }
    resultAll <- data.frame(result.mat)
    names(resultAll) <- c(rownames(refLocationProfiles),
        "value")
    row.names(resultAll) <- locs.vec
    # result.df <- data.frame(result.df,
    # locs.vec)
    #ind.min <- which.min(resultAll$value)

    #resultMin <- resultAll[ind.min, ]
    rOrder <- order(resultAll$value)
    resultOrder <- resultAll[rOrder, ]


    return(resultOrder)

}


