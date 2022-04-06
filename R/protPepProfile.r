#' Combine a protein profile with its component peptide profiles
#' 
#' Service function for protPepProfile. Creates a data frame containing 
#'    the profile of a single protein followed by 
#'    its component peptide profiles.
#' 
#' @param i Protein i
#' @param flagPeps data frame containing profiles for all spectra, 
#'         information mapping these to proteins and peptides, 
#'         and an indicator for outliers
#' @param numRefCols  number of reference columns
#' @param numDataCols number of fractions in each profile
#' @param protProfileData data frame containing protein profiles
#' @return Data frame containing the profiles for a single protein 
#'         and its component peptides
#' @examples
#' set.seed(7336)
#' eps <- 0.029885209
#' data(spectraNSA_test)
#' flagSpectraBox <- outlierFind(protClass=spectraNSA_test,
#'                        outlierLevel='peptide', numRefCols=5, numDataCols=9,
#'                        outlierMeth='boxplot', range=3, eps=eps,
#'                        randomError=TRUE)
#'                         
#'  #examine numbers of spectra that are outliers
#' table(flagSpectraBox$outlier.num.spectra)
#' 
#' # find peptide-level profiles
#' pepProfiles <- profileSummarize(protsCombineCnew=flagSpectraBox,
#'       numRefCols=6, numDataCols=9, refColsKeep=c(1,2,4),eps=eps,
#'       GroupBy="peptideId", outlierExclude="spectra")
#'  # find oulier peptides
#' flagPepsBox <- outlierFind(protClass=pepProfiles,
#'       outlierLevel="protein",
#'       numRefCols=3, numDataCols=9, eps=eps)
#'       str(flagPepsBox, strict.width="cut", width=65)
#'       
#'  # combine the two types of outliers into one
#' flagSpectraPeps <- merge(x=flagSpectraBox,
#'       y=flagPepsBox[,c(4,17)], by="pepId")
#'       
#' # Create protein profile from component peptide and spectra
#' protNSA_1 <-
#'       profileSummarize(protsCombineCnew=flagSpectraPeps,
#'       numRefCols=7, numDataCols=9, eps=eps, GroupBy="protId",
#'       outlierExclude="spectraAndPeptide")
#'       
#' # combine protein and peptide profiles 
#' protPepProfileNSA <- protPepProfile_i(i=1, flagPeps=flagPepsBox,
#'       numRefCols=4, numDataCols=9,
#'       protProfileData=protNSA_1)
#' str(protPepProfileNSA, strict.width="cut", width=65)
#' 
#' @export
protPepProfile_i <- function(i, flagPeps, numRefCols,
    numDataCols, protProfileData) {
    # byProtein, eps) { flagPeps: one row per
    # peptide col 1: prot col 2: peptide col
    # numRefCols + 1: N (followed by M, L1, L2,
    # ... for numDataCols columns)
    # protProfileData: one row per protein; must
    # match proteins in flagPeps col 1: prot col
    # 2: N, M, L1, L2, ...  This program will
    # insert a copy of the protein in 'peptide'
    # after prot

    # testing: flagPeps <- tannousflagPeps
    # protProfileData <- tannousprotProfileData
    # i=1

    prot.i <- protProfileData$prot[i]

    Npep <- 1  # always only one peptide for the peptide-level data
    protDataNarrow.i <- protProfileData[protProfileData$prot ==
        prot.i, ]  # original form for protein data
    # blankInsert <- cbind(rep(0, Nspectra),
    # matrix(NA, nrow=nSpectra,
    # numDataCols=(numRefCols-3))) Create a
    # one-row matrix of protein data; put protein
    # name first (for 'prot' slot) and second
    # (for 'peptide' slot) This will allow for
    # using the second column for row names later
    # protData.i is formatted to combine with
    # pepData.i
    protData.i <- data.frame(rep(protDataNarrow.i[1],
        2), matrix(NA, ncol = (numRefCols - 2), nrow = 1),
        protDataNarrow.i[2:(1 + numDataCols + 2)])



    # drop protId and pepId columns
    pepData.i <- flagPeps[flagPeps$prot == prot.i,
        -((numRefCols + numDataCols + 3):(numRefCols +
            numDataCols + 4))]
    Nspectra <- nrow(pepData.i)
    # pepMatInsert <- data.frame(matrix(NA,
    # nrow=Nspectra, numDataCols=2))
    # names(pepMatInsert) <- c('Nspectra.Pep',
    # 'Npep.Pep') pepData.i <-
    # data.frame(pepDataT.i, pepMatInsert)

    names(protData.i) <- names(pepData.i)

    protPepData.i <- rbind(protData.i, pepData.i)
    protPepData.i
}

#'Interlace protein profiles and component peptide profiles 
#' 
#' Interlace profiles of proteins and their component peptides into a 
#' single data frame.
#' @param flagPeps data frame containing profiles for all spectra, 
#'       information mapping these to proteins and peptides, 
#'       and an indicator for outliers
#' @param numRefCols  number of reference columns
#' @param numDataCols number of data columns
#' @param protProfileData data frame containing protein profiles
#' @return Data frame containing protein and peptide profiles
#' @examples
#' set.seed(17356)
#' eps <- 0.029885209
#' data(spectraNSA_test)
#' flagSpectraBox <- outlierFind(protClass=spectraNSA_test,
#'                     outlierLevel='peptide', numRefCols=5, numDataCols=9,
#'                     outlierMeth='boxplot', range=3, eps=eps,
#'                     randomError=TRUE)
#'                     
#'  #examine numbers of spectra that are outliers
#' table(flagSpectraBox$outlier.num.spectra)
#' # find peptide-level profiles
#' pepProfiles <- profileSummarize(protsCombineCnew=flagSpectraBox,
#'       numRefCols=6, numDataCols=9, refColsKeep=c(1,2,4),eps=eps,
#'       GroupBy="peptideId", outlierExclude="spectra")
#'       
#'  # find oulier peptides
#' flagPepsBox <- outlierFind(protClass=pepProfiles,
#'       outlierLevel="protein",
#'       numRefCols=3, numDataCols=9, eps=eps)
#'       str(flagPepsBox, strict.width="cut", width=65)
#'       
#'  # combine the two types of outliers into one
#' flagSpectraPeps <- merge(x=flagSpectraBox,
#'       y=flagPepsBox[,c(4,17)], by="pepId")
#'       
#'  # create protein profile from component peptide and spectra
#' protNSA_1 <-
#'       profileSummarize(protsCombineCnew=flagSpectraPeps,
#'       numRefCols=7, numDataCols=9, eps=eps, GroupBy="protId",
#'       outlierExclude="spectraAndPeptide")
#'       
#' # combine protein and peptide profiles
#' protPepProfileNSA <- protPepProfile(flagPeps=flagPepsBox,
#'       numRefCols=4, numDataCols=9,
#'       protProfileData=protNSA_1)
#' str(protPepProfileNSA, strict.width="cut", width=65)
#'   # See Vignette 6 for a full explanation
#' @export
#'
protPepProfile <- function(flagPeps, numRefCols, numDataCols,
    protProfileData) {

    # insert a copy of the protein in 'peptide'
    # after prot
    n.prots <- nrow(protProfileData)
    protPepData <- NULL
    for (i in seq_len(n.prots)) {
        protPepData.i <- protPepProfile_i(i, flagPeps,
            numRefCols, numDataCols, protProfileData)
        # byProtein, eps)
        protPepData <- rbind(protPepData, protPepData.i)

    }
    protPepData
}
