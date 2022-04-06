#' Plot set of simulated proteins distributed between two compartments.
#' 
#' For a set of simulated proteins distributed between two compartments, 
#'       plot the CPA estimates from a given type of profile
#'       versus true mixture proportions.
#'   
#' @param mixProtiProtjCPA data frame of CPA estimates for each 
#'       set of two compartment mixtures
#' @param NstartMaterialFractions Number of fractions that reconstitute 
#'       the starting material, e.g., a complete set of differential 
#'       centrifugation fractions.  For experiment AT5, it is 6 
#'       ( N, M, L1, L2, P, and S).
#' @param Loc1  row number of subcellular location 1 of mixture
#' @param Loc2  row number of subcellular location 2 of mixture
#' @param increment.prop  increments in proportion residing in 
#'     (from 0 to 1);
#'     Default is 0.1
#' @param errorReturn  return area of error region if TRUE
#' @param subTitle  subtitle for plot if present; NULL if not (default)
#' @param xaxisLab  plot label for x-axis if TRUE
#' @param yaxisLab  plot label for y-axis if TRUE
#' @return Plot CPA estimates versus true mixture proportions.


#' @examples
#' data(protNSA_test)
#' data(markerListJadot)
#' data(totProtAT5)
#' data(refLocProfAcup)
#' data(refLocProfNSA)
#' 
#' mixCytoLysoAcup <- proteinMix(AcupRef= refLocProfAcup,
#'                               increment.prop=0.1,
#'                               Loc1=1, Loc2=4)
                              
#' # Convert mixture Acup profiles to NSA profiles
#' mixCytoLysoNSA <- NSAfromAcup(Acup=mixCytoLysoAcup,
#'                       NstartMaterialFractions=6, totProt=totProtAT5)
#'                       
#' # Find constrained proportional values
#' mixCytoLysoCPAfromNSA <- fitCPA(profile=mixCytoLysoNSA,
#'                                 refLocationProfiles=refLocProfNSA,
#'                                 numDataCols=9)
#' 
#' # Plot the mixtures and fitted values and print the error
#' library(pracma)
#' mixturePlot(mixProtiProtjCPA=mixCytoLysoCPAfromNSA,
#'             NstartMaterialFractions=6, Loc1=1, Loc2=4,
#'             increment.prop=0.1, xaxisLab=TRUE, yaxisLab=TRUE,
#'             errorReturn = TRUE)
#'             
#' @importFrom graphics axis
#' @importFrom graphics points
#' @importFrom graphics segments
#' @importFrom graphics title
#' @importFrom pracma trapz
#' @export

mixturePlot <- function(mixProtiProtjCPA, NstartMaterialFractions = 6,
    Loc1, Loc2, increment.prop = 0.1, errorReturn = FALSE,
    subTitle = NULL, xaxisLab = TRUE, yaxisLab = TRUE) {
    # mixtures must be a list of equally spaced
    # proportions this program assumes exactly
    # eight subcellular compartments set up color
    # and point lists
    Loc1 <- as.integer(Loc1)
    Loc2 <- as.integer(Loc2)
    loc.list <- names(mixProtiProtjCPA)
    n.loc <- length(loc.list)
    col.list <- c("red", "blue", "orange", "darkgreen",
        "orange", "lightblue", "purple", "green")
    pch.list <- c(1, 2, 3, 4, 17, 6, 15, 8)
    col.list <- col.list[seq_len(n.loc)]
    pch.list <- pch.list[seq_len(n.loc)]
    plotLables <- data.frame(loc.list, col.list, pch.list)

    fracList <- seq(0, 1, increment.prop)  # for creating empty plot
    fracList2 <- fracList  #also need this for the empty plot

    # plot estimated proportion vs true
    # proportion
    plot(fracList ~ fracList2, type = "n", xlab = "",
        ylab = "", axes = FALSE)  # this is the empty plot
    if (yaxisLab)
        axis(2, labels = TRUE, at = c(0, 0.5, 1), las = 1)
    if (!yaxisLab)
        axis(2, labels = FALSE)
    if (xaxisLab)
        axis(1, labels = TRUE, at = c(0, 0.5, 1))
    if (!xaxisLab)
        axis(1, labels = FALSE)
    for (kk in seq_len(ncol(mixProtiProtjCPA))) {
        points(mixProtiProtjCPA[, kk] ~ fracList, col = col.list[kk],
            pch = pch.list[kk], cex = 0.85)
        # calculate sum of squares of errors

    }
    # make matrix input.prop listing mixtures
    prop.vec <- seq(0, 1, increment.prop)
    qrop.vec <- 1 - prop.vec
    nrow.out <- length(prop.vec)
    # mixMat <- matrix(0, nrow=nrow.out,
    # ncol=nrow(Acup)) # matrix of mixtures
    mixMat <- matrix(0, nrow = nrow.out, ncol = nrow.out)  # matrix of mixtures
    # each row is a 'protein' with mixed
    # residence each column is a subcellular
    # location, with the proportioal assignment
    mixMat[, Loc1] <- prop.vec
    mixMat[, Loc2] <- qrop.vec
    input.prop <- data.frame(mixMat)

    # each row of 'input.prop' is a mixture of p
    # of Loc1 1 and (1-p) of Loc2 Column Loc1 is
    # incremental increases in p from 0 to 1
    # Column Loc2 is incremental decreases in p
    # from 1 to 0 all other columns are 0 For
    # example, for Loc 1 = 1 (Cyto) and Loc2 = 4
    # (Lyso), here is the matrix.  row names are
    # not added here, since they are not needed.

    # Cyto Cyto1 Cyto2 Cyto3 Cyto4 Cyto5 Cyto6
    # Cyto7 Cyto8 Cyto9 Cyto10 0_Cyto:1_Lyso 0.0
    # 0 0 1.0 0 0 0 0 0 0 0 0.1_Cyto:0.9_Lyso 0.1
    # 0 0 0.9 0 0 0 0 0 0 0 0.2_Cyto:0.8_Lyso 0.2
    # 0 0 0.8 0 0 0 0 0 0 0 0.3_Cyto:0.7_Lyso 0.3
    # 0 0 0.7 0 0 0 0 0 0 0 0.4_Cyto:0.6_Lyso 0.4
    # 0 0 0.6 0 0 0 0 0 0 0 0.5_Cyto:0.5_Lyso 0.5
    # 0 0 0.5 0 0 0 0 0 0 0 0.6_Cyto:0.4_Lyso 0.6
    # 0 0 0.4 0 0 0 0 0 0 0 0.7_Cyto:0.3_Lyso 0.7
    # 0 0 0.3 0 0 0 0 0 0 0 0.8_Cyto:0.2_Lyso 0.8
    # 0 0 0.2 0 0 0 0 0 0 0 0.9_Cyto:0.1_Lyso 0.9
    # 0 0 0.1 0 0 0 0 0 0 0 1_Cyto:0_Lyso 1.0 0 0
    # 0.0 0 0 0 0 0 0 0



    areaErr <- 0
    for (kk in seq_len(ncol(mixProtiProtjCPA))) {
        # kk=1 use trapezoidal rule (trapz from
        # package pracma) areas are, first, area
        # under a triangle (assuming equally
        # spaced increments) and second, area
        # under the fitted CPA curve The
        # difference is the area between the
        # theoretical and observed mixture

        areaErr <- areaErr + abs(trapz(fracList, as.numeric(input.prop[,
            kk])) - trapz(fracList, as.numeric(mixProtiProtjCPA[,
            kk])))



    }

    segments(0, 0, 1, 1, col = col.list[Loc1])
    segments(0, 1, 1, 0, col = col.list[Loc2])
    segments(0, 0, 1, 0, col = "gray")
    titleText <- paste(loc.list[Loc1], "-", loc.list[Loc2])
    title(paste(titleText, "(", format(round(areaErr,
        3), nsmall = 3), ")\n", subTitle))  # guarantee 3 digits after decimal
    # plotLables loc.list col.list pch.list 1
    # Cyto red 1 open circle 2 ER blue 2 triangle
    # 3 Golgi orange 3 plus 4 Lyso darkgreen 4 X
    # 5 Mito orange 17 solid triangle 6 Nuc
    # lightblue 6 upside down triangle 7 Perox
    # purple 15 solid square 8 PM green 8
    # asterisk
    if (errorReturn) {
        areaErrOut <- data.frame(loc.list[Loc1], loc.list[Loc2],
            areaErr)
        names(areaErrOut) <- c("Loc1", "Loc2", "ErrorArea")
        return(areaErrOut)
    }
}
