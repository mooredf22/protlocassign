
#' plot mixture of all two compartment profile combinations as panel;
#'    assumes eight compartments
#' @param Acup relative amount of a given cellular
#'    compartment protein that ends up in a given centrifugation fraction
#' @param totProt vector of amounts starting material in each fraction
#' @param NstartMaterialFractions number of fractions in starting material
#' @param eps small positive constant to add before taking a log transformation
#' @return color heat map of errors
#' @examples
#'   # See See Vignette 4 for a full explanation
#' @importFrom graphics par
#' @importFrom graphics text
#' @importFrom viridisLite magma
#' @import plot.matrix
#' @import grid
#' @import gridExtra
#' @export


mixtureHeatMap <- function(Acup, totProt, NstartMaterialFractions = 6,
    eps = 0.001) {

    numDataCols <- ncol(Acup)
    nCompart <- nrow(Acup)
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
                errorMat <- matrix(c(ae11, ae12, ae13,
                    ae21, ae22, ae23), nrow = 2, byrow = TRUE)
                # add up all errors
                errorMatAll <- errorMatAll + errorMat



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
    par(op)  # restore par
    return(errorMatAll)
}


