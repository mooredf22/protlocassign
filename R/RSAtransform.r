#' Convert NSA to Acup
#' 
#' Convert NSA (normalized specific amount) profiles to Acup 
#' (relative amount) profiles
#'
#' @param NSA A matrix giving the specific amount or
#'         normalized specific amount
#'         profiles, either marker profiles from
#'         'cpaSetup' (which requires
#'         normalized specific amounts as input), or a list of
#'         many protein profiles
#' @param NstartMaterialFractions Number of fractions that reconstitute 
#'       the starting material, e.g., a complete set of differential 
#'       centrifugation fractions.  For experiment AT5, it is 6 
#'       ( N, M, L1, L2, P, and S).
#' @param totProt Total protein counts in each of the differential
#'         and nycodenz fractions; this is necessary to compute RSA's
#' @return Acup profiles
#'   amount of given protein in fraction /
#'   amount of that given protein in starting material
#' @examples
#' data(protNSA_test)
#' data(totProtAT5)
#' protAcup <- AcupFromNSA(protNSA_test[,seq_len(9)],
#'    NstartMaterialFractions = 6, totProt=totProtAT5)
#' protRSA <- RSAfromAcup(protAcup, NstartMaterialFractions = 6,
#'   totProt=totProtAT5)
#' protNSA <- NSAfromRSA(protRSA)
#' all.equal(protNSA, protNSA_test[,seq_len(9)] )
#' protRSA2 <- RSAfromNSA(protNSA_test[,seq_len(9)],
#'    NstartMaterialFractions = 6, totProt=totProtAT5)
#' all.equal(protRSA, protRSA2)
#' @export

AcupFromNSA <- function(NSA, NstartMaterialFractions,
    totProt = NULL) {
    if (ncol(NSA) != length(totProt)) {
        warning("Error from relAmtTransform:
            no. of rows of NSA must match length of totProt\n")
    }
  # just differential fractions
    startMaterialFractions <- NSA[, seq_len(NstartMaterialFractions)]
    nTotFractions <- length(totProt)  # number of all fractions
    NSAfractions <- NSA[, seq_len(nTotFractions)]
         # columns with NSA amounts; excludes additional columns

    AA <- data.frame(sweep(NSAfractions, 2, totProt,
        "*"))
    names(AA) <- colnames(NSAfractions)

    # Compute amount of given protein in fraction
    # / amount of given protein in starting
    # material use these values to create
    # mixtures these are in red in Excel
    # spreadsheet
    Acup <- data.frame(t(apply(AA, 1,
        function(x) x/sum(x[seq_len(NstartMaterialFractions)]))))
    names(Acup) <- colnames(AA)

    Acup
}



#' Convert Acup to RSA
#' 
#' Convert Acup (relative amount) profiles to RSA 
#' (relative specific amount) profiles
#' @param Acup amount of given protein in fraction /
#'             amount of that given protein in starting material
#' @param NstartMaterialFractions Number of fractions that reconstitute 
#'       the starting material, e.g., a complete set of differential 
#'       centrifugation fractions.  For experiment AT5, it is 6 
#'       ( N, M, L1, L2, P, and S).
#' @param totProt Total protein counts in each of the differential
#'          and nycodenz fractions; this is necessary to compute RSA's
#' @return rsa: relative specific amount
#' @examples
#' data(protNSA_test)
#' data(totProtAT5)
#' protAcup <- AcupFromNSA(protNSA_test[,seq_len(9)],
#'      NstartMaterialFractions = 6, totProt=totProtAT5)
#' protRSA <- RSAfromAcup(protAcup, NstartMaterialFractions = 6,
#'      totProt=totProtAT5)
#' protNSA <- NSAfromRSA(protRSA)
#' all.equal(protNSA, protNSA_test[,seq_len(9)] )
#' @export
RSAfromAcup <- function(Acup, NstartMaterialFractions,
    totProt = NULL) {

    if (ncol(Acup) != length(totProt)) {
        warning("Error from RSAfromAcup: no. of rows of Acup
              must match length of totProt\n")
    }
    # # # # # # # # rename variables # # # # # #
    # #

    # Difp -> t.h propFrac -> t.cup

  # total protein in the differential fractions
    t.h <- sum(totProt[seq_len(NstartMaterialFractions)])
  # proportion of protein in the differential fractions (9 component vector)
    t.cup <- totProt/t.h

    rsa <- sweep(Acup, 2, 1/t.cup, "*")
    names(rsa) <- colnames(Acup)
    rsa <- as.data.frame(rsa)
    # The following, which standardizes rsa rows
    # to sum to one, returns the original
    # markerProfiles !!
    #NSAfromRSA <- t(apply(rsa, 1, function(x) x/sum(x)))
    # this is deprecated; no need

    rsa
}


#' Convert NSA to RSA
#' 
#' Convert NSA (normalized specific amount) profiles to RSA
#' (relative specific amount) profiles
#' @param NSA A matrix giving the specific amount or
#'           normalized specific amount
#'           profiles, either marker profiles from 'cpaSetup'
#'           (which requires
#'           normalized specific amounts as input),
#'           or a list of many protein profiles
#' @param NstartMaterialFractions Number of fractions that reconstitute 
#'       the starting material, e.g., a complete set of differential 
#'       centrifugation fractions.  For experiment AT5, it is 6 
#'       ( N, M, L1, L2, P, and S).
#' @param totProt Total protein counts in each of the fractions;
#'           this is necessary to compute RSA's
#'
#' @return rsa: relative specific activity
#' @examples
#' data(protNSA_test)
#' data(totProtAT5)
#' protAcup <- AcupFromNSA(protNSA_test[,seq_len(9)],
#'              NstartMaterialFractions = 6, totProt=totProtAT5)
#' protRSA <- RSAfromAcup(protAcup, NstartMaterialFractions = 6,
#'              totProt=totProtAT5)
#' protNSA <- NSAfromRSA(protRSA)
#' all.equal(protNSA, protNSA_test[,seq_len(9)] )
#' @importFrom stats complete.cases
#' @export

RSAfromNSA <- function(NSA, NstartMaterialFractions,
    totProt = NULL) {

    missing.rows <- NSA[!complete.cases(NSA), ]

    if (ncol(NSA) != length(totProt)) {
        warning("Error from rsaDirect: no. of rows of NSA
            must match length of totProt\n")
    }
    # first columns of NSA
    startMaterialFractions <- NSA[, seq_len(NstartMaterialFractions)]
    nTotFractions <- length(totProt)  # number of all fractions
    # columns with NSA amounts; excludes additional columns
    NSAfractions <- NSA[, seq_len(nTotFractions)]


    protAbund <- AcupFromNSA(NSAfractions,
                      NstartMaterialFractions = NstartMaterialFractions,
        totProt = totProt)
    Acup <- protAbund
    rsa <- RSAfromAcup(Acup,
        NstartMaterialFractions = NstartMaterialFractions,
        totProt = totProt)

    rsa
}

#' Convert RSA to NSA
#' 
#' convert RSA (relative specific amount) profiles to NSA
#' (normalized specific amount) profiles
#' @param RSA relative specific activity
#' @return NSA A matrix giving the specific amount
#'        or normalized specific amount
#'        profiles, either marker profiles from 'cpaSetup' (which requires
#'        normalized specific amounts as input),
#'            or a list of many protein profiles
#' @examples
#' data(protNSA_test)
#' data(totProtAT5)
#' protAcup <- AcupFromNSA(protNSA_test[,seq_len(9)],
#'            NstartMaterialFractions = 6, totProt=totProtAT5)
#' protRSA <- RSAfromAcup(protAcup, NstartMaterialFractions = 6,
#'           totProt=totProtAT5)
#' protNSA <- NSAfromRSA(protRSA)
#' all.equal(protNSA, protNSA_test[,seq_len(9)] )
#' @export
NSAfromRSA <- function(RSA) {
    NSA <- data.frame(t(apply(RSA, 1, function(x) x/sum(x))))
    NSA
}

#' Convert Acup to NSA
#' 
#' Convert Acup (relative amount) profiles to NSA
#'  (normalized specific amount) profiles
#' @param Acup amount of given protein in fraction /
#'             amount of that given protein in starting material
#' @param NstartMaterialFractions Number of fractions that reconstitute 
#'       the starting material, e.g., a complete set of differential 
#'       centrifugation fractions.  For experiment AT5, it is 6 
#'       ( N, M, L1, L2, P, and S).
#' @param totProt Total protein counts in each of the differential
#'         and nycodenz fractions; this is necessary to compute RSA's
#' @return NSA: normalized specific amount
#' @examples
#' data(protNSA_test)
#' data(totProtAT5)
#' protAcup <- AcupFromNSA(protNSA_test[,seq_len(9)],
#'             NstartMaterialFractions = 6, totProt=totProtAT5)
#' protNSA <- NSAfromAcup(protAcup, NstartMaterialFractions = 6,
#'             totProt=totProtAT5)
#' all.equal(protNSA, protNSA_test[,seq_len(9)])
#' @export
NSAfromAcup <- function(Acup, NstartMaterialFractions,
    totProt = NULL) {
    RSA <- RSAfromAcup(Acup,
        NstartMaterialFractions = NstartMaterialFractions,
        totProt = totProt)
    NSA <- NSAfromRSA(RSA)
    NSA
}

#' Convert RSA to Acup
#' 
#' Convert RSA (relative specific amount) profiles to Acup
#' (relative amount) profiles
#' @param RSA amount of given protein in fraction /
#'            amount of that given protein in starting material
#' @param NstartMaterialFractions Number of fractions that reconstitute 
#'       the starting material, e.g., a complete set of differential 
#'       centrifugation fractions.  For experiment AT5, it is 6 
#'       ( N, M, L1, L2, P, and S).
#' @param totProt Total protein counts in each of the differential
#'            and nycodenz fractions; this is necessary to compute RSA's
#' @return Acup: relative amount
#' @examples
#' data(protNSA_test)
#' data(totProtAT5)
#' protAcup <- AcupFromNSA(protNSA_test[,seq_len(9)],
#'             NstartMaterialFractions = 6, totProt=totProtAT5)
#' protRSA <- RSAfromAcup(protAcup, NstartMaterialFractions = 6,
#'             totProt=totProtAT5)
#' protNSA <- NSAfromRSA(protRSA)
#' all.equal(protNSA, protNSA_test[,seq_len(9)] )
#' @export
AcupFromRSA <- function(RSA, NstartMaterialFractions,
    totProt = NULL) {
    NSA <- NSAfromRSA(RSA)
    Acup <- AcupFromNSA(NSA, NstartMaterialFractions = NstartMaterialFractions,
        totProt = totProt)
    Acup
}
