#' Set up reference profiles for constrained proportional assignment
#'
#' @param profile data frame of protein identifiers (protName)
#'        and their relative abundance in centrifugation fractions.
#' @param markerList List of reference proteins and their subcellular locations
#' @param numDataCols Number of channels of abundance levels
#' @return A matrix refLocationProfiles giving the abundance level profiles of
#'         the subcellular locations
#' @import knitr
#' @import rmarkdown
#' @export
#' @examples
#' data(protNSA_AT5tmtMS2)
#' data(markerListJadot)
#' refLocationProfilesNSA <- locationProfileSetup(profile=protNSA_AT5tmtMS2,
#'   markerList=markerListJadot, numDataCols=9)


locationProfileSetup <- function(profile, markerList,
    numDataCols) {
    # Find the mean profiles of the subcellular
    # locations

    # 'profile' must be a data frame with
    # proteins as row names and these columns:
    # columns 1 - numDataCols: normalized
    # specific amounts or specific amounts Last
    # two columns: Nspectra (number of spectra)
    # and Nseq (number of distince peptide
    # sequences)

    # 'markerList' must be a list of (1) the
    # reference proteins and (2) the
    # corresponding subcelluar locations

    if (!is.data.frame(profile))
        print("profile must be a data frame")
    if (!is.data.frame(markerList))
        print("markerList must be a data frame")
    if (is.null(names(markerList)[1]))
        print("first markerList variable must be 'protName'")
    if (is.null(names(markerList)[2]))
        print("second markerList variable must be 'referenceCompartment'")
    if ({
        names(markerList)[1] != "protName"
    }) {
        print("first markerList variable must be 'protName'")
    }
    if ({
        names(markerList)[2] != "referenceCompartment"
    }) {
        print("second markerList variable must be 'referenceCompartment'")
    }
    # make names all upper case
    protNamesUpper <- toupper(rownames(profile))

    rownames(profile) <- protNamesUpper
    profileWithProts <- data.frame(protNamesUpper,
        profile)
    # this has a column of protnames for merge:
    names(profileWithProts)[1] <- "protName"

    names(markerList)[1] <- "protName"
    # all must by upper case:
    markerList$protName <- toupper(markerList$protName)
    if (nrow(markerList) >= 1) {
        meanReferenceProts <- merge(x = markerList,
            y = profileWithProts, by.x = "protName",
            by.y = "protName", all.x = FALSE, sort = FALSE)
    }
    if (nrow(meanReferenceProts) == 0) {
        print("Error from locationProfileSetup; no protein found")
    }

    # Find mean profiles for each sub-cellular
    # location, using reference proteins

    location.list <-
        as.character(unique(meanReferenceProts$referenceCompartment))
    n.loc <- length(location.list)
    meanProfile <- NULL
    # for (i in 1:n.loc) {
    for (i in seq_len(n.loc)) {
      # i=1
      loc.i <- location.list[i]

      # profile.i <-
      # meanReferenceProts[meanReferenceProts$referenceCompartment
      # == loc.i,2+1:numDataCols]
      profile.i <-
            meanReferenceProts[meanReferenceProts$referenceCompartment ==
            loc.i, 2 + seq_len(numDataCols)]
      meanProfile.i <- apply(profile.i, 2, mean,
            na.rm = TRUE)
      meanProfile <- rbind(meanProfile, meanProfile.i)
    }
    row.names(meanProfile) <- location.list

    # Change name to 'markerLoc' and its
    # transpose 'refLocationProfiles' for
    # compatibility with past Then plot the
    # profiles

    markerLoc <- t(meanProfile)
    refLocationProfiles <- as.data.frame(t(markerLoc))
    refLocationProfiles
}



#' internal function; goodness-of-fit measure
#' @param pvec vector of assignment proportions to
#'             compartments; to be estimated; must sum to 1
#' @param y observed data for amounts of protein in fractions;
#'        must sum to 1
#' @param gmat reference protein amounts
#' @param methodQ either 'sumsquares' (default) or 'sumabsvalue'
#' @export
#' @return value of goodness-of-fit function
#' @examples
#' # suppose there are four fractions and three compartments
#' yy <- c(0.1, 0.2, 0.1, 0.6)
#' gmat <- matrix(c(0,   0,   1,   # ref prot 1, all in compartment 3
#'                  0,   1,   0,   # ref prot 2, all in compartment 2
#'                  0,   0.5, 0.5, # ref prot 3, half in comp 2, half in 3)
#'                  0.5, 0,   0.5),# ref prot 4, half in comp 1, half in 3
#'                  ncol=3, byrow=TRUE )
#' pvec <- c(0, 0.3, 0.7)  # 30%  in comp 2, 70% in comp 3
#' Qfun4(pvec, yy, gmat)
Qfun4 <- function(pvec, y, gmat, methodQ = "sumsquares") {
    # Assign sub-cellular location probabilities
    # to each protein. We use the 'spg' function
    # in package 'BB'.

    resultA <- y - pvec %*% t(gmat)
    if (methodQ == "sumsquares")
        result <- sum(resultA^2)
    if (methodQ == "sumabsvalue")
        result <- sum(abs(resultA))
    result
}
#' internal function; goodness-of-fit measure for
#' subset of varying parameters
#' @param pvec.vary vector of assignment proportions to
#'         compartments; to be estimated
#' @param yy observed data for amounts of protein in
#'           fractions
#' @param gmat reference protein amounts
#' @param methodQ either 'sumsquares' (default) or
#'                'sumabsvalue'
#' @param ind.vary if not NULL, indexes of varying parameters
#' @param ind.fixed indexes of fixed parameters; all parameters
#'                 must be one or the other
#' @param par.fixed values of fixed parameters; typically 0
#' @export
#' @return value of goodness-of-fit function
#' @examples
#' # suppose there are four fractions and three compartments
#' yy <- c(0.1, 0.2, 0.1, 0.6)
#' gmat <- matrix(c(0,   0,   1,   # ref prot 1, all in compartment 3
#'                  0,   1,   0,   # ref prot 2, all in compartment 2
#'                  0,   0.5, 0.5, # ref prot 3, half in comp 2, half in 3)
#'                  0.5, 0,   0.5),# ref prot 4, half in comp 1, half in 3
#'                  ncol=3, byrow=TRUE )
#' pvec <- c(0, 0.3, 0.7)  # 30%  in comp 2, 70% in comp 3
#' pvec.vary <- c(.3, 0.7)  # just the compartments that vary
#' ind.vary <- c(2,3)  # compartments 2 and 3 may vary
#' ind.fixed <- 1   # fix value of compartment 1
#' par.fixed <- 0   # fix that value at 0
#' Qfun4subset(pvec.vary, yy, gmat, ind.vary=ind.vary,
#'            ind.fixed=ind.fixed, par.fixed=par.fixed)
#' # Now allow only compartments 1 and 3 to vary
#' Qfun4subset(pvec.vary, yy, gmat, ind.vary=c(1,3),
#'            ind.fixed=2, par.fixed=0)
#'
Qfun4subset <- function(pvec.vary, yy, gmat,
    methodQ = "sumsquares",
    ind.vary, ind.fixed, par.fixed) {
    # Assign sub-cellular location probabilities
    # to each protein. We use the 'spg' function
    # in package 'BB'.

    # reassemble
    vals.order <- c(ind.vary, ind.fixed)


    pvec.unorder <- c(pvec.vary, par.fixed)
    n.vec <- length(pvec.unorder)
    pvec.orig <- rep(NA, n.vec)
    # for (i in 1:n.vec) {
    for (i in seq_len(n.vec)) {
        pvec.orig[vals.order[i]] <- pvec.unorder[i]
    }

    y <- yy
    pvec <- pvec.orig
    resultA <- y - pvec %*% t(gmat)

    if (methodQ == "sumsquares")
        result <- sum(resultA^2)
    if (methodQ == "sumabsvalue")
        result <- sum(abs(resultA))
    result
}

# internal function; project an n-dim vector y to
# the simplex S_n with sum constrained to 1
# @param y n-dim vector
# @return simplex of
# dimension one lower
# @export
projSimplex <- function(y) {
    # project an n-dim vector y to the
    # simplex S_n S_n = { x | x \in R^n, 0 <=
    # x <= 1, sum(x) = 1} Derived from:
    # Ravi Varadhan, Johns Hopkins
    # University August 8, 2012
    ##
    ## See also:
  # Projection onto a simplex
  # Yunmei Chen, Xiaojing Ye -
  # arXiv preprint arXiv:1101.6081, 2011 - arxiv.org
  #
# Matlab version (Xiaojing Ye):
#  https://www.mathworks.com/matlabcentral/fileexchange/
  #        30332-projection-onto-simplex
    #####

    n <- length(y)
    sy <- sort(y, decreasing = TRUE)
    csy <- cumsum(sy)
    rho <- max(which(sy > (csy - 1)/seq_len(n)))
    theta <- (csy[rho] - 1)/rho
    return(pmax(0, y - theta))
}


#' return index of a protein name, or (if exactMatch=TRUE) indices of
#'    proteins starting with the string given in 'protName'
#' @param protName  name of protein to search for
#' @param profile data frame of protein identifiers (protName) and
#'        their relative abundance in centrifugation fractions
#' @param exactMatch  default is F
#' @export
#' @return the protein name and its index (row in profile)
#' @examples
#' data(protNSA_AT5tmtMS2)
#' protIndex('TLN1', profile=protNSA_AT5tmtMS2)

protIndex <- function(protName, profile, exactMatch = FALSE) {
    # return index of a protein name, or (if
    # exactMatch=TRUE) indices of proteins
    # starting with the string given in
    # 'protName'
    n.prot <- nrow(profile)

    prot.list <- toupper(as.character(rownames(profile)))
    # prot.list must be in column 1
    if (exactMatch)
        inx <- (seq_len(n.prot))[protName == prot.list]
    if (!exactMatch)
        inx <- grep(paste("^", protName, sep = ""),
            prot.list, ignore.case = TRUE)
    if (length(inx) == 0)
        inx <- NA

    if (!is.na(inx[1])) {
        result <- data.frame(inx, prot.list[inx])
        names(result) <- c("Prot index number", "Prot name")
    }
    if (is.na(inx[1])) {
        result <- NA
        cat("protein not found\n")
    }
    result
}

#' Carry out constrained proportional assignment for protein i;
#'    service function for fitCPA
#' @param profile data frame of one protein name (row name) and
#'        relative abundance levels for that protein
#' @param refLocationProfiles A matrix giving the abundance
#'        level profiles of the subcellular locations
#' @param numDataCols Number of channels of abundance levels
#' @param startProps starting valuese for proportional assignements;
#'        set equal if this is null (default)
#' @param maxit maximum number of iterations (default is 10000)
#' @param ind.vary if not NULL, indexes of proportions allowed to vary
#' @param minVal default is false. If true, return minimum value of
#'               goodness of fit
#' @examples
#' data(protNSA_test)
#' data(markerListJadot)
#' nTestProts <- nrow(protNSA_test)
#' refLocationProfilesNSA <- locationProfileSetup(profile=protNSA_test,
#'     markerList=markerListJadot, numDataCols=9)
#' protCPAfromNSA_test <- fCPAone(profile=protNSA_test[1,],
#'                           refLocationProfiles=refLocationProfilesNSA,
#'                           numDataCols=9, startProps=NULL,
#'                           maxit=10000,
#'                          ind.vary=NULL, minVal=FALSE)
#' @export
#' @return assignProbsOut  Data frame of proportionate assignments of
#'    each protein to compartments. also nspectra and npeptides if in
#'    profile input
#'    If ind.vary specified, only the referenced proportion cpa
#'    estimates are returned

# protLocAssign <- function(i, profile,
# refLocationProfiles, numDataCols,
fCPAone <- function(profile, refLocationProfiles, numDataCols,
    startProps = NULL, maxit = 10000,
    ind.vary = NULL, minVal = FALSE) {

    # maxit and assignPRobsStart must be
    # specified assignProbsStart must be NULL or
    # have a column 'protName' and assignment
    # probabilities to use as starting values use
    # the spg function (in package BB) to assign
    # proportionate assignments to compartments

  # extra columns for numbers of spectra and sequences:
    SpectraSeqInd <- TRUE
    if (numDataCols == ncol(profile))
        SpectraSeqInd <- FALSE  # no extra columns


    n.compartments <- nrow(refLocationProfiles)

    # just the data, and only row 1:
    allPeptideProfilesMat <- profile[1, seq_len(numDataCols)]
    # (there should only be one row)

    yyT <- allPeptideProfilesMat



    yy <- as.numeric(yyT)  # leave them alone


    if (SpectraSeqInd) {
        # if these variables are present
        Nspectra.i <- profile$Nspectra[1]  # this is the number
                            # of spectra for a protein
        Npep.i <- profile$Npep[1]  # number of unique sequences
    }
    if (!SpectraSeqInd) {
        Nspectra.i <- NULL
        Nseq.i <- NULL
    }
    channelsMeanProb.i <- matrix(NA, nrow = 1, ncol = length(yy))
    parEstTemp <- channelsMeanProb.i
    # if (!anyNA(yy)) { start with uniform
    # probabilities:
    startProps <- rep(1/n.compartments, n.compartments)
    if (!is.null(startProps)) {
        if (length(startProps) != n.compartments) {
            cat("invalid startProps\n")
        }
        if (length(startProps) == n.compartments)
            startVals <- startProps/sum(startProps)
    }
    # start with uniform probabilities startVals
    # <- as.numeric(assignProbsStart[i,
    # 2:(n.compartments + 1)]) # start with

    # first attempt at minimization: use uniform
    # starting values ordinary procedure; find
    # optimum over all parameters
    if (is.null(ind.vary)) {
        temp <- try(BB::spg(startVals, fn = Qfun4,
            project = projSimplex, y = yy, gmat = t(refLocationProfiles),
            methodQ = "sumsquares", quiet = TRUE, method = 3,
            control = list(maxit = maxit, trace = FALSE,
                ftol = 1e-12, gtol = 1e-07, eps = 1e-09)), silent=TRUE)
    }
    if (!is.null(ind.vary)) {
        # only find optimum over varying
        # parameters
        ind.vary <- sort(ind.vary)  # parameters that may vary
        # here are the parameters that are fixed
        # (complement of ind.vary) ind.fixed <-
        # (1:length(startVals))[!((1:length(startVals))
        # %in% ind.vary)]
        ind.fixed <-
          (seq_len(length(startVals)))[!((seq_len(length(startVals))) %in%
            ind.vary)]
        # yy.vary <- yy[ind.vary] yy.fixed <-
        # yy[ind.fixed] startVals.vary <-
        # startVals[ind.vary] # params that vary
        startVals.vary <- rep(1/length(ind.vary), length(ind.vary))

        if (length(ind.fixed) > 0)
            startVals.fixed <- rep(0, length(ind.fixed))
        if (length(ind.fixed) == 0)
            startVals.fixed <- NULL
        temp <- try(BB::spg(par = startVals.vary, fn = Qfun4subset,
            project = projSimplex, yy = yy, gmat = t(refLocationProfiles),
            methodQ = "sumsquares", ind.vary = ind.vary,
            ind.fixed = ind.fixed, par.fixed = startVals.fixed,
            quiet = TRUE, alertConvergence = FALSE,
            method = 3, control = list(maxit = maxit,
                trace = FALSE, ftol = 1e-12, gtol = 1e-07,
                eps = 1e-09)), silent=TRUE)
    }
    if (!is.atomic(temp))
        convergeInd <- as.numeric((temp$message ==
            "Successful convergence"))

    # If starting values are given
    # (assignProbsStart is not null) and if the
    # first attempt failed, try the starting
    # values if ({convergeInd != 1} &
    # {!is.null(assignProbsStart)}) { temp <-
    # try(BB::spg(startVals, fn=Qfun4,
    # project=projSimplex, y=yy,
    # gmat=t(refLocationProfiles),
    # methodQ='sumsquares', quiet=T,
    # control=list(maxit=maxit, trace=FALSE)))

    # }

    convergeInd <- as.numeric(!inherits(temp, "try-error"))

    if (is.null(ind.vary))
        channelsMeanProb.i <- rep(NA, n.compartments)
    if (!is.null(ind.vary))
        channelsMeanProb.i <- rep(NA, length(ind.vary))

    convCrit.i <- NA
    feval.i <- NA
    value <- NA
    if (!is.atomic(temp)) {
        channelsMeanProb.i <- temp$par
        convCrit.i <- temp$gradient
        feval.i <- temp$feval
        value <- temp$value
        convergeInd <- as.numeric((temp$message ==
            "Successful convergence"))
    }
    # nNoConverge.i <- 0
    if (convergeInd != 1)
        cat(paste("cpa does not converge for a protein",
           "\n", "returning missing values for cpa estimates for that protein",
          "\n"))


    if (!is.null(ind.vary)) {
        vals.order <- c(ind.vary, ind.fixed)
        pvec.vary <- channelsMeanProb.i
        par.fixed <- rep(0, n.compartments - length(ind.vary))
        pvec.unorder <- c(pvec.vary, par.fixed)
        n.vec <- length(pvec.unorder)
        pvec.orig <- rep(NA, n.vec)
        # for (i in 1:n.vec) {
        for (i in seq_len(n.vec)) {
            pvec.orig[vals.order[i]] <- pvec.unorder[i]
        }
        parEstTemp <- pvec.orig

    }
    if (is.null(ind.vary))
        parEstTemp <- channelsMeanProb.i
    # }


    if (SpectraSeqInd)
        parEst <- c(parEstTemp, Nspectra.i, Npep.i)
    if (!SpectraSeqInd)
        parEst <- parEstTemp
    if (minVal)
        parEst <- c(parEst, value)

    return(parEst)
}




#' Carry out constrained proportional assignment on all proteins
#'
#' @param profile a data frame of protein names (row names) and
#'       relative abundance levels.
#' @param refLocationProfiles A matrix giving the abundance level
#'       profiles of the subcellular locations
#' @param numDataCols Number of channels of abundance levels
#' @param startProps starting valuese for proportional assignements;
#'       set equal if this is null (default)
#' @param maxit maximum number of iterations (default is 10000)
#' @param showProgress print out progress if TRUE, the default
#' @param ind.vary if not NULL, indexes of parameters to allow to vary
#' @param minVal default is false. If true, return minimum value of
#'       goodness of fit
#' @export
#' @return assignProbsOut  Data frame of proportional assignments
#'       of each protein to compartments
#' @examples
#' data(protNSA_test)
#' data(markerListJadot)
#' refLocationProfilesNSA <- locationProfileSetup(profile=protNSA_test,
#'        markerList=markerListJadot, numDataCols=9)
#' protCPAfromNSA_test <- fitCPA(profile=protNSA_test,
#'                            refLocationProfiles=refLocationProfilesNSA,
#'                            numDataCols=9)

fitCPA <- function(profile, refLocationProfiles, numDataCols,
    startProps = NULL, maxit = 10000, showProgress = TRUE,
    ind.vary = NULL, minVal = FALSE) {


    # # # # # # # # # # # #
    n.prot <- nrow(profile)
    # indList <- 1:n.prot
    indList <- seq_len(n.prot)
    SpectraSeqInd <- TRUE  # extra columns for numbers of spectra and sequences
    if (numDataCols == ncol(profile))
        SpectraSeqInd <- FALSE  # no extra columns

    # result <- sapply(indList, fCPAone,
    # profile=profile,
    # refLocationProfiles=refLocationProfiles,
    # numDataCols=numDataCols,
    # startProps=startProps, showProgress=FALSE,
    # simplify=T, maxit=maxit, ind.vary, minVal)

    assignProbs <- NULL
    # for (i in 1:n.prot) {
    for (i in seq_len(n.prot)) {
        assignProbsI <- fCPAone(profile = profile[i,
            ], refLocationProfiles = refLocationProfiles,
            numDataCols = numDataCols, startProps = startProps,
            maxit = maxit, ind.vary = ind.vary, minVal = minVal)
        assignProbs <- rbind(assignProbs, assignProbsI)
        if (showProgress) {
            if (i == 500)
                print(paste(i, "profiles fit"))
            if ((i%%1000) == 0)
                print(paste(i, "profiles fit"))
        }
    }
    # browser()
    assignProbs <- data.frame(assignProbs)
    # if (is.null(ind.vary)) ind.vary <-
    # 1:nrow(refLocationProfiles)
    if (is.null(ind.vary))
        ind.vary <- seq_len(nrow(refLocationProfiles))
    # (so the next check works either with or
    # without ind.vary specified) checkCols <-
    # {ncol(assignProbs) >=
    # nrow(refLocationProfiles[ind.vary,]) + 2}

    # if (checkCols) {
    if (SpectraSeqInd) {
        # if (is.null(ind.vary)) ind.vary <-
        # 1:nrow(refLocationProfiles) # for when
        # ind.vary not specified
        # names(assignProbs) <-
        # c(row.names(refLocationProfiles[ind.vary,]),
        # 'Nspectra', 'Npeptides')
        names(assignProbs)[seq_len(nrow(refLocationProfiles) +
            2)] <- c(row.names(refLocationProfiles),
            "Nspectra", "Npeptides")
        protNames <- rownames(profile)  # make sure it is character
        assignProbsOut <- data.frame(assignProbs)
        rownames(assignProbsOut) <- protNames[indList]
        # name minVal column if present
        if (minVal)
            names(assignProbsOut)[ncol(assignProbs)] <- "value"
            # name of last column
    }
    if (!SpectraSeqInd) {
        # names(assignProbs)[1:(nrow(refLocationProfiles))]
        # <- row.names(refLocationProfiles)
        names(assignProbs)[seq_len(nrow(refLocationProfiles))] <-
          row.names(refLocationProfiles)
        protNames <- rownames(profile)
        assignProbsOut <- data.frame(assignProbs)
        rownames(assignProbsOut) <- protNames[indList]
        # name minVal column if present
        if (minVal)
            names(assignProbsOut)[ncol(assignProbs)] <- "value"
    }
    assignProbsOut
}


#' convert single protein name to capitalize first character only
#' @param x a protein name
#' @return single protein name to capitalize first character only
#' @export
#' @examples
#' protCap('tpp1')
protCap <- function(x) {
    # convert single protein name to capitalize
    # first character only Typically, human
    # proteins are all in capitals and rat are
    # capitalized (first character only) with the
    # rest in lower case.
    firstLetter <- toupper(substring(x, 1, 1))
    remainingLetters <- tolower(substring(x, 2))
    paste(firstLetter, remainingLetters, sep = "")
}

#' assign proteins to a subcellular location using CPA estimates
#'
#' @param assignLocProps matrix of proportion estimates for each protein
#' @param cutoff cutoff for assigning a protein to a location
#' @param Locations list of subcellular locations
#' @return protLoc assigned location of protein
#' @examples
#' data(protNSA_test)
#' data(markerListJadot)
#' refLocationProfilesNSA <- locationProfileSetup(profile=protNSA_test,
#'          markerList=markerListJadot, numDataCols=9)
#' protCPAfromNSA_test <- fitCPA(profile=protNSA_test,
#'          refLocationProfiles=refLocationProfilesNSA,
#'          numDataCols=9)
#' Locations <- unique(markerListJadot$referenceCompartment)
#' table(apply(protCPAfromNSA_test[,1:8],1,assignCPAloc,
#'          cutoff=0.8, Locations=Locations))
#' @importFrom stats complete.cases
#' @export

assignCPAloc <- function(assignLocProps, cutoff = 0.8,
    Locations) {
    # assign location to the one with a
    # proportion > cutoff assignLocProps <-
    # assignPropsAll[1,]
    if (cutoff <= 0.5)
        cutoff <- 0.5001
    exceedCutoff <- {
        assignLocProps > cutoff
    }
    if (sum(exceedCutoff) > 0) {
        indLoc <- which(exceedCutoff)
        protLoc <- Locations[indLoc]
    }
    if (sum(exceedCutoff) == 0)
        protLoc <- "unclassified"
    protLoc
}

# unitize functions # # #
#' normalize vector xx to have unit length
#' @param xx vector
#' @return nomalized vector of unit length
#' @examples
#' xx <- c(0.5, 0.1, 0.6, 0.9)
#' vecUnitize(xx)
#' @export
vecUnitize <- function(xx) {
    if (anyNA(xx))
        xx.norm <- rep(NA, length(xx)) else {
        xxsq <- xx^2
        xx.norm <- sqrt(sum(xxsq))
        xx.norm <- xx/xx.norm
    }
    xx.norm
}

#' normalize all rows of a matrix to have unit length
#' @param protMatOrig matrix of profiles
#' @return matrix with all rows having unit length
#' @examples
#' data(protNSA_test)
#' vectorizeAll(protNSA_test[,1:9])
#' @export
vectorizeAll <- function(protMatOrig) {
    protUnitMat <- NULL
    # for (i in 1:nrow(protMatOrig)) {
    for (i in seq_len(nrow(protMatOrig))) {
        # i=1
        temp <- vecUnitize(protMatOrig[i, ])
        protUnitMat <- rbind(protUnitMat, temp)
    }
    row.names(protUnitMat) <- row.names(protMatOrig)
    protUnitMat
}

