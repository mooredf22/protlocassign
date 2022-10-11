
#' Heat map of mixture CPA errors
#' 
#' Produce heat map of errors in CPA estimates for fitting simulated 
#'    proteins distributed between all pairs of compartments using 
#'    different types of profiles (RSA, NSA, and Acup; 
#'    linear and log2-transformed). 
#'    Formatted for use with up to eight compartments , 
#'    which produces 28 (8 chose 2) pairs of combinations of compartments. 
#'    Also writes a table of total assignment errors for each type of 
#'    transformation.
#' 
#' @param Acup Acup profiles of reference compartments
#' @param totProt vector of total protein amounts (derived from a given 
#'         amount of starting material) in each of the fractions 
#'         comprising the profile
#' @param NstartMaterialFractions Number of fractions that reconstitute 
#'       the starting material, e.g., a complete set of differential 
#'       centrifugation fractions.  For experiment AT5, it is 6 
#'       ( N, M, L1, L2, P, and S).
#' @param eps small value to add so that log argument is greater than zero
#' @param errorListIn If not NULL, pre-computed errors are passed to the
#'     function to save computing time for display purposes.
#'      The format must be as
#'     in the $errorList component in the output. Default is NULL
#' @param errorListReturn  If TRUE, return data frame of mixture errors
#'     if FALSE (the default) return matrix of total errors
#' @return Plot of heat map of errors of CPA estimates for pairs of 
#'       simulated protein mixtures represented by different types of 
#'       profiles and a table of total errors for each type of transformation.
#'         Errors are displayed in a 3x2 output where the top row of the 
#'         left, middle and right columns represent RSA, NSA, 
#'         and Acup profiles, respectively and the bottom row represents 
#'         a log2-transformation of these profiles. 
#'         By default return a matrix total errors. If errorListShow=TRUE,
#'         return data frame of mixture errors
#'    
#' @examples
#' \donttest{
#' data(totProtAT5)
#' data(refLocProfAcup)
#' par(mfrow=c(1,1))
#' errorResult <- mixtureHeatMap(Acup=refLocProfAcup, totProt=totProtAT5)
#' round(head(errorResult), digits=3)
#' }
#' # Note that the profile of one protein, AIF1, contains missing values
#' # which causes the cpa routine to generate a nonconvergence message 


#' @importFrom graphics par
#' @importFrom graphics text
#' @importFrom viridisLite magma
#' @import plot.matrix
#' @import grid
#' @import gridExtra
#' @export


mixtureHeatMap <- function(Acup, totProt, NstartMaterialFractions = 6,
    eps = 0.001, errorListIn=NULL, errorListReturn=FALSE) {

    numDataCols <- ncol(Acup)
    nCompart <- nrow(Acup)
    compartNames <- rownames(Acup)
    refLocationProfilesAcup <- Acup
    refLocationProfilesNSA <- NSAfromAcup(refLocationProfilesAcup,
        NstartMaterialFractions = NstartMaterialFractions,
        totProt = totProt)
    refLocationProfilesRSA <- RSAfromNSA(NSA = refLocationProfilesNSA,
        NstartMaterialFractions = NstartMaterialFractions,
        totProt = totProt)

    op <- par(mfrow = c(3, 3))
    protAmtList <- c("Relative specific amount\n(RSA)",
        "Normalized specific amount\n(NSA)", "Relative amount\n(Acup)",
        "Log2 RSA", "Log2 NSA", "Log2 Acup")
    protAmtMat <- matrix(protAmtList, nrow = 2, byrow = TRUE)


    x <- c(0, 1)
    y <- c(0, 1)
    nComp <- nrow(Acup)  # number of compartments
    par(mfrow = c(nComp, nComp), mar = c(0.5, 0.5,
        0.5, 0.5), xpd = NA)
    compartmentList <- rownames(Acup)
    # print column labels
    for (k in seq_len(nComp)) {
        plot(y ~ x, type = "n", axes = FALSE, xlab = "",
            ylab = "")
        if (k != 1)
            text(0.5, 0.5, compartmentList[k], cex = 2)
    }

    # start with a 2 by 3 matrix of zeros; this
    # will contain the sum of all errors matrices
    errorMatAll <- matrix(0, nrow = 2, ncol = 3, byrow = TRUE)
    errorList <- NULL
    jj <- 0
    for (i in seq_len((nComp - 1))) {
        for (j in seq_len(nComp)) {


            if (j <= i) {
                plot(y ~ x, type = "n", axes = FALSE,
                    xlab = "", ylab = "")


                if ((j == 1) & (i != nCompart))
                    text(0.5, 0.5, compartmentList[i],
                    cex = 2)
            }
            if (j > i) {
             if (is.null(errorListIn)) {
                # i=1 j=4 create mixture
                mixProtiProtjAcup <- proteinMix(refLocationProfilesAcup,
                    Loc1 = i, Loc2 = j)
                mixProtiProtjRSA <- RSAfromAcup(Acup = mixProtiProtjAcup,
                    NstartMaterialFractions = NstartMaterialFractions,
                    totProt = totProt)

                # find RSA
                refLocationProfilesRSA <-
                    RSAfromNSA(NSA = refLocationProfilesNSA,
                    NstartMaterialFractions = NstartMaterialFractions,
                    totProt = totProt)
                # apply CPA to relative specific
                # amounts (recommended way)
                mixProtiProtjProp <- fitCPA(profile = mixProtiProtjRSA,
                    refLocationProfiles = refLocationProfilesRSA,
                    numDataCols = numDataCols)

                # normalized specific amounts
                # mixProtiProtjSpecAmt <-
                # t(apply(mixProtiProtjRSA,1,
                # function(x) x/sum(x)))
                mixProtiProtjNSA <- NSAfromRSA(mixProtiProtjRSA)

                # apply CPA to normalized
                # specific amounts
                # mixProtiProtjPropSpecAmt <-
                # fitCPA(profile=mixProtiProtjSpecAmt,
                mixProtiProtjPropNSA <- fitCPA(profile = mixProtiProtjNSA,
                    refLocationProfiles = refLocationProfilesNSA,
                    numDataCols = numDataCols)

                # apply CPA to relative amounts
                # (Acup)
                mixProtiProtjPropAcup <- fitCPA(profile = mixProtiProtjAcup,
                    refLocationProfiles = refLocationProfilesAcup,
                    numDataCols = numDataCols)

                # log transformed

                mixProtiProtjPropLog2 <-
                    fitCPA(profile = log2(mixProtiProtjRSA +
                    eps), refLocationProfiles = log2(refLocationProfilesRSA +
                    eps), numDataCols = numDataCols)


                mixProtiProtjPropNSALog2 <-
                    fitCPA(profile = log2(mixProtiProtjNSA +
                    eps), refLocationProfiles = log2(refLocationProfilesNSA +
                    eps), numDataCols = numDataCols)

                mixProtiProtjPropAcupLog2 <-
                    fitCPA(profile = log2(mixProtiProtjAcup +
                    eps), refLocationProfiles = log2(refLocationProfilesAcup +
                    eps), numDataCols = numDataCols)



                # # # #

            ae11 <- mixtureAreaError(mixProtiProtjCPA = mixProtiProtjProp,
                   NstartMaterialFractions = NstartMaterialFractions,
                   Loc1 = i, Loc2 = j, increment.prop = 0.1)
            ae12 <- mixtureAreaError(mixProtiProtjCPA = mixProtiProtjPropNSA,
                    NstartMaterialFractions = NstartMaterialFractions,
                    Loc1 = i, Loc2 = j, increment.prop = 0.1)
            ae13 <- mixtureAreaError(mixProtiProtjCPA = mixProtiProtjPropAcup,
                    NstartMaterialFractions = NstartMaterialFractions,
                    Loc1 = i, Loc2 = j, increment.prop = 0.1)



        ae21 <- mixtureAreaError(mixProtiProtjCPA = mixProtiProtjPropLog2,
                    NstartMaterialFractions = NstartMaterialFractions,
                    Loc1 = i, Loc2 = j, increment.prop = 0.1)
         ae22 <- mixtureAreaError(mixProtiProtjCPA = mixProtiProtjPropNSALog2,
                   NstartMaterialFractions = NstartMaterialFractions,
                    Loc1 = i, Loc2 = j, increment.prop = 0.1)
        ae23 <- mixtureAreaError(mixProtiProtjCPA = mixProtiProtjPropAcupLog2,
                      NstartMaterialFractions = NstartMaterialFractions,
                  Loc1 = i, Loc2 = j, increment.prop = 0.1)
                errorList_ij <- c(i, j, ae11, ae12, ae13, ae21, ae22, ae23)
             
               
                errorMat <- matrix(c(ae11, ae12, ae13,
                    ae21, ae22, ae23), nrow = 2, byrow = TRUE)
             }
                if (!is.null(errorListIn)) {
                  jj <- jj+1
                  errorMat <- matrix(as.numeric(errorListIn[jj,]), nrow=2,
                                     byrow=TRUE)
                  errorList_ij <- as.numeric(c(i, j, errorListIn[jj,]))
                }
                # add up all errors
                errorMatAll <- errorMatAll + errorMat
                errorList <- rbind(errorList, errorList_ij)
             


               errorMat1m <- 2 - errorMat  # so that higher errors are dark red
                # library(plot.matrix)
                requireNamespace("plot.matrix")
                # library('viridis')


                col <- rev(viridis::magma(20))
                plot(errorMat, col = col, breaks = seq(0,
                  2, 0.1), key = NULL, main = "", axis.col = NULL,
                  axis.row = NULL, xlab = "", ylab = "",
                  digits = 2, cex = 0.8)
            }
        }
    }
    errorListTemp <- data.frame(errorList)
    list.i <- as.numeric(errorListTemp[,1])
    list.j <- as.numeric(errorListTemp[,2])
    mixNames <- paste(compartmentList[list.i], compartmentList[list.j])
    errorList <- errorListTemp[,3:8]
    names(errorList) <- c("RSA", "NSA", "Acup",
                         "log RSA", "log NSA", "log Acup")
    rownames(errorList) <- mixNames
    par(op)  # restore par
    if (errorListReturn) result <- errorList
    if (!errorListReturn) result <- errorMatAll
    #result <- list(errorMatAll=errorMatAll, errorList=errorList)
    return(result)
}


