#' compute a mixture of proteins in two compartments
#' @param AcupRef amount of given protein in fraction /
#'                amount of given protein in starting material
#' @param Loc1  row number of one compartment
#' @param Loc2  row number of other compartment
#' @param increment.prop fraction increment from 0 to 1
#' @return mixAmount relative amounts of proteins in the fractions
#' @examples
#' data(protNSA_test)
#' data(markerListJadot)
#' data(totProtAT5)
#' refLocationProfilesNSA <- locationProfileSetup(profile=protNSA_test,
#'                                 markerList=markerListJadot, numDataCols=9)
#' round(refLocationProfilesNSA, digits=3)
#' # Convert NSA reference profiles to Acup to prepare for forming mixtures
#' refLocationProfilesAcup <-
#'      AcupFromNSA(NSA=refLocationProfilesNSA, NstartMaterialFractions=6,
#'      totProt=totProtAT5)
#' round(refLocationProfilesAcup, digits=4)
#' # Compute mixtures
#' mixCytoLysoAcup <- proteinMix(AcupRef=refLocationProfilesAcup,
#'                               increment.prop=0.1,
#'                               Loc1=1, Loc2=4)
#' @export
proteinMix <- function(AcupRef, Loc1, Loc2, increment.prop = 0.1) {
    # replace mix.df with input.prop throughout
    # replace relAmount with Acup througout
    nrowRef <- nrow(AcupRef)
    if (Loc1 > Loc2) {
        Loc1orig <- Loc1
        Loc2orig <- Loc2
        Loc1 <- Loc2orig
        Loc2 <- Loc1orig
    }
    if (Loc2 > nrowRef)
        warning("Error, not enough rows\n")
    LocNames <- row.names(AcupRef)
    prop.vec <- seq(0, 1, increment.prop)
    qrop.vec <- 1 - prop.vec
    nrow.out <- length(prop.vec)
    Acup <- NULL
    mixProtNames <- NULL

    for (i in seq_len(nrow.out)) {
        Acup.i <- prop.vec[i] * AcupRef[Loc1, ] + qrop.vec[i] *
            AcupRef[Loc2, ]
        Acup <- rbind(Acup, Acup.i)
        mixProtNames.i <- paste(prop.vec[i], "_", LocNames[Loc1],
            ":", qrop.vec[i], "_", LocNames[Loc2],
            sep = "")
        mixProtNames <- c(mixProtNames, mixProtNames.i)

    }

    # mixMat <- matrix(0, nrow=nrow.out,
    # ncol=nrow(Acup)) # matrix of mixtures each
    # row is a 'protein' with mixed residence
    # each column is a subcellular location, with
    # the proportioal assignment mixMat[,Loc1] <-
    # prop.vec mixMat[,Loc2] <- qrop.vec
    # input.prop <- data.frame(mixMat)
    # names(input.prop) <- row.names(Acup)
    # rownames(input.prop) <- mixProtNames

    row.names(Acup) <- mixProtNames
    # result <- list(Acup=Acup,
    # input.prop=input.prop) result
    Acup
}


