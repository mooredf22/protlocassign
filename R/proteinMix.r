#' Simulate sets of protein profiles distributed between two compartments
#' 
#' Compute the Acup profiles of simulated proteins that are distributed 
#' between two compartments in specified proportions
#'  
#' @param AcupRef data frame containing Acup profiles for the 
#'      reference compartments
#' @param Loc1  row number of one compartment
#' @param Loc2  row number of other compartment
#' @param increment.prop increments in proportion residing in
#'        Loc1(from 0 to 1); Default is 0.1
#' @return A data frame containing Acup profiles for the simulated proteins
#' @examples
#' 
#' data(refLocProfAcup)
#' 
#' # Compute mixtures
#' mixCytoLysoAcup <- proteinMix(AcupRef=refLocProfAcup,
#'                               increment.prop=0.1,
#'                               Loc1=1, Loc2=4)
#' round(mixCytoLysoAcup, digits=4)
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
        warning("not enough rows\n")
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


    row.names(Acup) <- mixProtNames
 
    Acup
}


