#' Restricted CPA
#' 
#' Carry out constrained proportional assignment for protein i, 
#'   restricting the number of compartments that are allowed to vary 
#'   freely and fixing the CPA values of the others at 0
#'   
#' @param profile a vector containing a specified protein (row name) 
#'      profile ‘Nspectra’ and ‘Npep’, if present, will be removed
#' @param refLocationProfiles data frame of profiles for the  
#'      reference compartments
#' @param numDataCols number of fractions in each profile
#' @param startProps starting values for proportional assignments;
#'    set equal if this is null (default)
#' @param showProgress   default is TRUE
#' @param maxit maximum number of iterations (default is 10000)
#' @param nCPAcomparts number of compartments to fit restricted CPA;
#'        remaining proportions are fixed at zero
#' @return Data frame with CPA estimates for every combination
#'        ordered by value, with best fit listed first

#' @examples
#' data(protRSA_test)
#' data(refLocProfRSA)
#' protCPAfromRSA_out2 <- fCPAsubsets(profile=protRSA_test[1,],
#'                          refLocationProfiles=refLocProfRSA,
#'                          numDataCols=9, startProps=NULL, nCPAcomparts=2)
#' round(head(protCPAfromRSA_out2), digits=4)
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
    
    rOrder <- order(resultAll$value)
    resultOrder <- resultAll[rOrder, ]


    return(resultOrder)

}


