#' Convert NSA to Acup
#' 
#' Convert NSA (normalized specific amount) profiles to Acup 
#' (relative amount) profiles
#'
#' @param NSA matrix of normalized specific amount profiles 
#'    (see help file for protNSA_test)
#' @param NstartMaterialFractions number of fractions that reconstitute 
#'       the starting material, e.g., a complete set of differential 
#'       centrifugation fractions.  For experiment AT5, it is 6 
#'       ( N, M, L1, L2, P, and S).
#' @param totProt vector of total protein amounts (derived from a given amount 
#'   of starting material) in each of the fractions comprising the profile
#' @return Acup profiles;
#'   amount of given protein in fraction /
#'   amount of that given protein in starting material
#'   (see help file for protAcup_test).
#' @examples
#' data(protNSA_test)
#' data(totProtAT5)
#' protAcup_out1 <- AcupFromNSA(protNSA_test[,seq_len(9)],
#'    NstartMaterialFractions = 6, totProt=totProtAT5)
#' round(head(protAcup_out1), digits=4)
#' str(protAcup_out1)
#' @export

AcupFromNSA <- function(NSA, NstartMaterialFractions,
    totProt = NULL) {
    if (ncol(NSA) != length(totProt)) {
        warning("from relAmtTransform:
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
#' @param Acup data frame of specified protein (row name) Acup profiles;
#'   amount of given protein in fraction /
#'   amount of that given protein in starting material
#'   (see help file for protAcup_test).
#' @param NstartMaterialFractions Number of fractions that reconstitute 
#'       the starting material, e.g., a complete set of differential 
#'       centrifugation fractions.  For experiment AT5, it is 6 
#'       ( N, M, L1, L2, P, and S).
#' @param totProt vector of total protein amounts (derived from a 
#'         given amount of starting material) in each of the 
#'         fractions comprising the profile 
#' @return RSA profiles.  RSA is the ratio of two ratios: 
#'    the numerator is the amount of a given protein in a particular 
#'    fraction divided by the amount of that given protein in the 
#'    starting material while the denominator is amount of total protein 
#'    in a particular fraction divided by the amount of total protein 
#'    in the starting material. The RSA describes the fold-enrichment 
#'    (RSA>1) or depletion (RSA<1) of a protein during the 
#'    fractionation process, and is analogous to the 
#'    relative specific activity term used in classical 
#'    analytical subcellular fractionation
#' @examples
#' data(protAcup_test)
#' data(totProtAT5)
#' protRSA_out1 <- RSAfromAcup(protAcup_test[,seq_len(9)], 
#'      NstartMaterialFractions = 6,
#'      totProt=totProtAT5)
#' round(head(protRSA_out1), digits=4)
#' str(protRSA_out1)
#' @export
RSAfromAcup <- function(Acup, NstartMaterialFractions,
    totProt = NULL) {

    if (ncol(Acup) != length(totProt)) {
        warning("from RSAfromAcup: no. of rows of Acup
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
#' @param NSA data frame of specified protein(row name) NSA profiles
#'       (see help file for protNSA_test)
#' @param NstartMaterialFractions Number of fractions that reconstitute 
#'       the starting material, e.g., a complete set of differential 
#'       centrifugation fractions.  For experiment AT5, it is 6 
#'       ( N, M, L1, L2, P, and S).
#' @param totProt vector of total protein amounts 
#'           (derived from a given amount of starting material) 
#'           in each of the fractions comprising the profile 
#'
#' @return RSA profiles.  RSA is the ratio of two ratios: 
#'    the numerator is the amount of a given protein in a particular 
#'    fraction divided by the amount of that given protein in the 
#'    starting material while the denominator is amount of total protein 
#'    in a particular fraction divided by the amount of total protein 
#'    in the starting material. The RSA describes the fold-enrichment 
#'    (RSA>1) or depletion (RSA<1) of a protein during the 
#'    fractionation process, and is analogous to the 
#'    relative specific activity term used in classical 
#'    analytical subcellular fractionation
#' @examples
#' data(protNSA_test)
#' data(totProtAT5)
#' protRSA_out2 <- RSAfromNSA(protNSA_test[,seq_len(9)], 
#'              NstartMaterialFractions = 6,
#'              totProt=totProtAT5)
#' round(head(protRSA_out2), digits=4)
#' str(protRSA_out2)
#' @importFrom stats complete.cases
#' @export

RSAfromNSA <- function(NSA, NstartMaterialFractions,
    totProt = NULL) {

    missing.rows <- NSA[!complete.cases(NSA), ]

    if (ncol(NSA) != length(totProt)) {
        warning("from rsaDirect: no. of rows of NSA
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
#' Convert RSA (relative specific amount) profiles to NSA
#' (normalized specific amount) profiles
#' @param RSA data frame containing RSA profiles.  
#'    RSA is the ratio of two ratios: 
#'    the numerator is the amount of a given protein in a particular 
#'    fraction divided by the amount of that given protein in the 
#'    starting material while the denominator is amount of total protein 
#'    in a particular fraction divided by the amount of total protein 
#'    in the starting material. The RSA describes the fold-enrichment 
#'    (RSA>1) or depletion (RSA<1) of a protein during the 
#'    fractionation process, and is analogous to the 
#'    relative specific activity term used in classical 
#'    analytical subcellular fractionation.
#' @return A data frame of NSA profiles 
#'    (see help file for protNSA_test)
#' @examples
#' data(protRSA_test)
#' protNSA_out2 <- NSAfromRSA(protRSA_test[,seq_len(9)])
#' round(head(protNSA_out2), digits=4)
#' str(protNSA_out2)
#' @export
NSAfromRSA <- function(RSA) {
    NSA <- data.frame(t(apply(RSA, 1, function(x) x/sum(x))))
    NSA
}

#' Convert Acup to NSA
#' 
#' Convert Acup (relative amount) profiles to NSA
#'  (normalized specific amount) profiles
#' @param Acup data frame containing Acup profiles;
#'   amount of given protein in fraction /
#'   amount of that given protein in starting material
#'   (see help file for protAcup_test).
#' @param NstartMaterialFractions number of fractions that reconstitute 
#'       the starting material, e.g., a complete set of differential 
#'       centrifugation fractions.  For experiment AT5, it is 6 
#'       ( N, M, L1, L2, P, and S).
#' @param totProt vector of total protein amounts (derived from a given 
#'        amount of starting material) in each of the fractions 
#'        comprising the profile
#' @return A data frame of NSA profiles 
#'    (see help file for protNSA_test)
#' @examples
#' data(protAcup_test)
#' data(totProtAT5)
#' protNSA_out1 <- NSAfromAcup(protAcup_test[,seq_len(9)], 
#'             NstartMaterialFractions = 6,
#'             totProt=totProtAT5)
#' round(head(protNSA_out1), digits=4)
#' str(protNSA_out1)
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
#' @param RSA matrix of relative specific amount profiles 
#'      (see help file for protRSA_test)  
#' @param NstartMaterialFractions number of fractions that reconstitute 
#'       the starting material, e.g., a complete set of differential 
#'       centrifugation fractions.  For experiment AT5, it is 6 
#'       ( N, M, L1, L2, P, and S).
#' @param totProt vector of total protein amounts (derived from 
#'         a given amount of starting material) in each of the 
#'         fractions comprising the profile 
#' @return Acup profiles: amount of given protein in fraction / 
#'         amount of that given protein in starting material 
#'         (see help file protAcup_test).
#' @examples
#' data(protRSA_test)
#' data(totProtAT5)
#' protAcup_out2 <- AcupFromRSA(protRSA_test[,seq_len(9)],
#'             NstartMaterialFractions = 6, totProt=totProtAT5)
#' round(head(protAcup_out2), digits=4)
#' str(protAcup_out2)
#' @export
AcupFromRSA <- function(RSA, NstartMaterialFractions,
    totProt = NULL) {
    NSA <- NSAfromRSA(RSA)
    Acup <- AcupFromNSA(NSA, NstartMaterialFractions = NstartMaterialFractions,
        totProt = totProt)
    Acup
}
