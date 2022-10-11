#' Area-based error of mixtures
#' 
#' Compute area-based error for accuracy of CPA estimates for 
#'    sets of simulated proteins distributed between two compartments.
#'    
#' @param mixProtiProtjCPA data frame of CPA estimates for each set 
#'      of two compartment mixtures
#' @param NstartMaterialFractions  Number of fractions that reconstitute 
#'       the starting material, e.g., a complete set of differential 
#'       centrifugation fractions.  For experiment AT5, it is 6 
#'       ( N, M, L1, L2, P, and S).
#' @param Loc1  row number of subcellular location 1 of mixture
#' @param Loc2  row number of subcellular location 2 of mixture
#' @param increment.prop  increments in proportion residing in Loc1
#'         (from 0 to 1); Default is 0.1
#' @return The area between predicted and observed curves
#' @examples
#' data(refLocProfAcup)
#' data(refLocProfNSA)
#' data(totProtAT5)
#'
#' # Compute relative amount of each theoretical protein that resides 
#' #    in each fraction in a mixture set
#' mixCytoLysoAcup <- proteinMix(AcupRef=refLocProfAcup,
#'                               increment.prop=0.1,
#'                               Loc1=1, Loc2=4)
#'                               
#' # Convert theoretical protein profiles to NSA
#' mixCytoLysoNSA <- NSAfromAcup(Acup=mixCytoLysoAcup,
#'                            NstartMaterialFractions=6, totProt=totProtAT5)
#'                               
#' # Find constrained proportional values
#' mixCytoLysoCPAfromNSA <- fitCPA(profile=mixCytoLysoNSA,
#'                                 refLocationProfiles=refLocProfNSA,
#'                                 numDataCols=9)
#' # calculate the mixture error
#' mixtureAreaError(mixProtiProtjCPA=mixCytoLysoCPAfromNSA,
#'             NstartMaterialFractions=6, Loc1=1, Loc2=4,
#'             increment.prop=0.1)
#' @importFrom pracma trapz
#' @export


mixtureAreaError <- function(mixProtiProtjCPA, NstartMaterialFractions = 6,
    Loc1, Loc2, increment.prop = 0.1) {
    # mixtures must be a list of equally spaced
    # proportions this program assumes exactly
    # eight subcellular compartments set up color
    # and point lists
    Loc1 <- as.integer(Loc1)
    Loc2 <- as.integer(Loc2)
    loc.list <- names(mixProtiProtjCPA)

    fracList <- seq(0, 1, 0.1)

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
        # package pracma)
        areaErr <- areaErr + abs(trapz(fracList, as.numeric(input.prop[,
            kk])) - trapz(fracList, as.numeric(mixProtiProtjCPA[,
            kk])))
    }

    return(areaErr)

}
