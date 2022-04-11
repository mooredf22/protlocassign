#' Find outlier spectra within a protein or pedtide i.
#' 
#' Service function for outlierFind (see help function for outlierFind).  
#' 
#' @param i unique identifier for protein i
#' @param protClass a data frame containing profiles associated 
#'            with either spectra or peptides (see Tutorial 6)
#' @param outlierLevel 'peptide' for outlier spectra within peptides, or
#'                     'protein' for outlier peptides within proteins
#' @param numRefCols number of columns (variables) before data columns
#' @param numDataCols number of fractions in each profile
#' @param outlierMeth boxplot (recommended), scores, or
#'        none (if no outliers are to be reported)
#' @param range the range parameter used for identifying outliers
#'              using the boxplot method
#' @param proba probability to exclude
#'         outlier for scores method
#' @param eps small value to add so that log argument is greater than zero
#' @param randomError TRUE if allow it to be random
#' @return  New data frame for a single protein or peptide with an 
#'       additional column that indicates the number of fractions in 
#'       a profile (peptides or spectra) that are outliers
#' @examples
#' set.seed(63561)
#' eps <- 0.029885209
#' data(spectraNSA_test)
#' uniqueLabel <- spectraNSA_test$pepId[1]
#' flagSpectraBox_i <- outlierFind_i(i=uniqueLabel, protClass=spectraNSA_test,
#'                       outlierLevel='peptide', numRefCols=5, numDataCols=9,
#'                       outlierMeth='boxplot', range=3, proba=NA, eps=eps,
#'                       randomError=TRUE)
#' str(flagSpectraBox_i)
#' 
#' @importFrom graphics boxplot
#' @importFrom outliers scores
#' @importFrom stats runif
#' @export

outlierFind_i <- function(i, protClass, outlierLevel = "peptide",
    numRefCols, numDataCols, outlierMeth, range, proba,
    eps = eps, randomError) {

    # i=9
    # ================================================================
    # boxplotMod to identify outliers in a vector
    # Outliers that are 3 times the interquartile
    # range from either quartile; note that the
    # default for boxplots is 1.5 times the
    # interquartile range @param x: a vector of
    # numbers @param reject.vec.i: pre-specified
    # rejection locations @param range:
    # pre-specified times of IQR range
    # ================================================================

    boxplotMod <- function(x, range = range) {
        # x[reject.vec] <- NA
        outlier.x <- boxplot(x, plot = FALSE, range = range)$out
        outlier.ind <- (x %in% outlier.x)

        return(outlier.ind * 1)
    }

    # ================================================================
    # scoresMod to identify outliers in a vector
    # Outliers as exceeding the normal
    # probability of proba @param x: a vector of
    # numbers @param reject.vec.i: pre-specified
    # rejection locations @param proba:
    # pre-specified normal probability
    # ================================================================

    scoresMod <- function(x, reject.vec = NULL, proba = proba) {
        x[reject.vec] <- NA
        non.na.indices <- (seq_len(length(x)))[!is.na(x)]
        ind <- scores(x[!is.na(x)], prob = proba)
        ind[is.na(ind)] <- FALSE
        outlier.ind <- rep(FALSE, length(x))
        outlier.ind[non.na.indices] <- ind

        return(outlier.ind * 1)
    }

    if (outlierLevel == "peptide")
        uniqueLabel <- protClass$pepId
    if (outlierLevel == "protein")
        uniqueLabel <- protClass$protId
    protClass.i <- protClass[uniqueLabel == i, ]
    nProt <- nrow(protClass.i)  #(or number of peptides)
    nczfMatLog.i <- log2(protClass.i[, (numRefCols +
        1):(numRefCols + numDataCols)] + eps)
    if (nProt > 1) {
        proteinsAllData <- nczfMatLog.i
        n.spectra <- nrow(proteinsAllData)

        outIndMat.i <- data.frame(matrix(NA, ncol = numDataCols,
            nrow = n.spectra))
        # out.boxplotMod.i <- outIndMat.i
        names(outIndMat.i) <- names(nczfMatLog.i)
        for (j in seq_len(numDataCols)) {
            # j=1
            xx <- as.numeric(proteinsAllData[, j])
            pp <- 2^xx - eps  # convert back to original scale
            # add uniform (0, eps) random number,
            # if 'randomError=TRUE'
            if (randomError == TRUE) {
                n.pp <- length(pp)
                eps.random <- runif(n = n.pp, min = 0,
                  max = eps) * (pp < 10^(-5))
                pp.r <- pp + eps.random
                xx <- log2(pp.r + eps)
            }

            # identify outliers that are 3 times
            # the interquartile range from either
            # quartile note that the default for
            # boxplots is 1.5 times the
            # interquartile range


            if (outlierMeth == "boxplot")
                outIndMat.i[, j] <- boxplotMod(xx,
                  range = range)
            if (outlierMeth == "scores")
                outIndMat.i[, j] <- scoresMod(xx, proba = proba)
            if (outlierMeth == "none")
                outIndMat.i[, j] <- 0 * length(xx)
        }
        outlier.num <- apply(outIndMat.i, 1, sum)

    }
    if (nProt == 1) {
        outlier.num <- 0
    }
    if (outlierLevel == "peptide")
        outlier.num.spectra <- outlier.num
    if (outlierLevel == "protein")
        outlier.num.peptides <- outlier.num
    if (nProt >= 1) {
        if (outlierLevel == "peptide")
            result <- cbind(protClass.i[, seq_len(numRefCols)],
                outlier.num.spectra, protClass.i[,
                  (numRefCols + 1):(numRefCols + numDataCols +
                    2)])  # include protId and pepId
        if (outlierLevel == "protein")
            result <- cbind(protClass.i[, seq_len(numRefCols)],
                outlier.num.peptides, protClass.i[,
                  (numRefCols + 1):(numRefCols + numDataCols +
                    4)])  # include protId and pepId
    }
    if (nProt == 0)
        result <- NULL
    names(result)[1] <- "prot"
    return(result)
}



# ================================================================
#' Identify outlier profiles
#' 
#' Identify outlier profiles.  This can be done at the level of 
#'     identifying outlier spectra when calculating peptide profiles or 
#'     identifying outlier peptides when calculating protein profiles. 
#'     See Tutorial 6 for a description of outlier determination methods.
#'
#' @param     protClass a data frame containing profiles associated 
#'            with either spectra or peptides (see Tutorial 6)
#' @param     outlierLevel 'peptide' for outlier spectra within peptides, or
#'                          'protein' for outlier peptides within proteins
#' @param     numRefCols number of columns (variables) before data columns
#' @param     numDataCols number of fractions in each profile
#' @param     outlierMeth boxplot (recommended), scores, or 
#'                   none (if no outliers are to be reported)
#' @param     range the range parameter used for identifying outliers 
#'              using the boxplot method
#' @param     proba probability to exclude outlier for scores method
#' @param     eps small value to add so that log argument is greater than zero
#' @param     randomError  TRUE if allow it to be random
#' @param     setSeed seed for random number generator
#' @param     cpus 1 (default);
#'            if cpus > 1 use BiocParallel with SnowParm(cpus)
#' @return    New data frame with an additional column that indicates 
#'           the number of fractions in a profile (spectra or peptide) 
#'           that are outliers
#' @examples
#' set.seed(17356)
#' eps <- 0.029885209
#' data(spectraNSA_test)
#' flagSpectraBox <- outlierFind(protClass=spectraNSA_test,
#'                               outlierLevel='peptide', numRefCols=5,
#'                               numDataCols=9,
#'                               outlierMeth='boxplot', range=3, eps=eps,
#'                               randomError=TRUE, cpus=2)
#'                               
#' # examine breakdown of spectral according to the number of fractions 
#' #  in their profiles that are outliers
#' table(flagSpectraBox$outlier.num.spectra)
#' 
#' @importFrom BiocParallel bplapply
#' @importFrom BiocParallel SnowParam
#' @export    outlierFind

outlierFind <- function(protClass, outlierLevel = "peptide",
    numRefCols = 5, numDataCols = 9, outlierMeth = "boxplot",
    range = 3, proba = 0.99, eps = eps, randomError = TRUE,
    setSeed=NULL, cpus=1) {


        if (outlierLevel == "peptide")
            uniqueLabel <- protClass$pepId
        if (outlierLevel == "protein")
            uniqueLabel <- protClass$protId

        indList <- unique(uniqueLabel)

    if (cpus > 1) {
        system.time(result <- bplapply(indList, outlierFind_i,
            protClass = protClass, outlierLevel = outlierLevel,
            numRefCols = numRefCols, numDataCols = numDataCols,
            outlierMeth = outlierMeth, range = range,
            proba = proba, eps = eps, randomError = randomError,
            BPPARAM = SnowParam(workers=cpus,RNGseed = setSeed)))
    }
        if (cpus == 1) {
        system.time(result <- lapply(indList, outlierFind_i,
            protClass = protClass, outlierLevel = outlierLevel,
            numRefCols = numRefCols, numDataCols = numDataCols,
            outlierMeth = outlierMeth, range = range,
            proba = proba, eps = eps, randomError = randomError))
    }

    protsCombineCnew <- do.call(what = "rbind", result)
        # convert list of matrices to one matrix
    protsCombineCnew
}
