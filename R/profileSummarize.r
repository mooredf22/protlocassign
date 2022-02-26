#' service function for profileSummarize
#'   find mean profiles for proteins, or for
#'   protein-peptide combinations for specific peptide or protein
#' @param i protID or pepID (only one)
#' @param uniqueLabel vector of either
#'   protId or pepId as specified by 'GroupBy' in
#'   profileSummarize
#' @param protsCombineCnew
#'   dataframe of protein profiles
#' @param numRefCols
#'   number of columns preceding the data
#' @param numDataCols number of numerical values in a
#'    profile
#' @param GroupBy 'protId' if average by
#'    protein; 'peptideId' if average by peptide
#' @param eps small value to add so that log
#'    argument is greater than zero
#' @param outlierExclude none, spectra, or
#'      spectraAndPeptide @return result_i coefficient
#'      estimates, numbers of peptides and spectra, and
#'      id for protein i
#' @return estimated means for all channels, and peptide or protein info
#' @examples
#' set.seed(17356)
#' eps <- 0.029885209
#' data(TLN1_test)
#' flagSpectraBox <- outlierFind(protClass=TLN1_test,
#'                               outlierLevel='peptide', numRefCols=5,
#'                               numDataCols=9,
#'                               outlierMeth='boxplot', range=3, eps=eps,
#'                               randomError=TRUE)
#' # examine numbers of spectra that are outliers
#' table(flagSpectraBox$outlier.num.spectra)
#' uniqueLabel <- flagSpectraBox$pepID
#' pepProfiles <- meansByProteins(i=uniqueLabel[1], uniqueLabel=uniqueLabel,
#'           protsCombineCnew=flagSpectraBox,
#'           numRefCols=6, numDataCols=9, GroupBy="peptideId",eps=eps,
#'           outlierExclude='spectra')
#' str(pepProfiles, strict.width='cut', width=65)
#' @importFrom stats var
#' @importFrom lme4 lmer
#' @importFrom lme4 fixef
#' @importFrom lme4 isSingular
#' @importFrom stats vcov
#' @importFrom stats sd
#' @importFrom stats median
#' @importFrom methods is
#
#' @export
meansByProteins <- function(i, uniqueLabel, protsCombineCnew,
    numRefCols, numDataCols, GroupBy, eps, outlierExclude) {
    # outlierExclude: none: don't exclude any
    # outliers spectra: exclude only spectra
    # within peptides spectraAndPeptide: exclude
    # spectra-within-peptide outliers and peptide
    # outliers i=217 i=2


    protData_i <- protsCombineCnew[uniqueLabel == i,
        ]
    prot_i <- protData_i$prot[1]
    protId_i <- protData_i$protId[1]
    pep_i <- protData_i$peptide[1]
    pepId_i <- protData_i$pepId[1]
    pepId_vec_i <- protData_i$pepId

    profileAll_i_x <- cbind(protData_i[, numRefCols +
        c(seq_len(numDataCols))], protData_i$protId,
        protData_i$pepId)
    names(profileAll_i_x)[numDataCols + 1] <- "protId"
    names(profileAll_i_x)[numDataCols + 2] <- "pepId"
    profileAll_i <- profileAll_i_x
    profileAll_i[, seq_len(numDataCols)] <- log2(profileAll_i_x[,
        seq_len(numDataCols)] + eps)  # transform to log2 scale



    # identify and eliminate outliers use_i <-
    # {protData_i$outCountBoxPlot == 0}
    if (outlierExclude == "none") {
        use_i <- rep(TRUE, nrow(protData_i))
    }
    if (outlierExclude == "spectra") {
        use_i <- {
            {
                protData_i$outlier.num.spectra == 0
            } & {
                !is.na(protData_i$outlier.num.spectra)
            }
        }
    }
    if (outlierExclude == "spectraAndPeptide") {
        use_i <- {
            {
                protData_i$outlier.num.spectra == 0
            } & {
                !is.na(protData_i$outlier.num.spectra)
            } & {
                protData_i$outlier.num.peptides ==
                  0
            } & {
                !is.na(protData_i$outlier.num.peptides)
            }
        }
    }
    profileXX_i <- profileAll_i[use_i, ]


    Nspectra <- nrow(profileXX_i)

    if (Nspectra == 0) {
        coef_est <- rep(NA, numDataCols)
        secoef_est <- rep(NA, numDataCols)
        Npep <- 0
        result_i <- c(coef_est, secoef_est, Nspectra,
            Npep, protId_i, pepId_i)
        return(result_i)
    }
    Npep <- length(unique(profileXX_i$pepId))


    # need at least 4 spectra and 3 unique
    # peptides
    lmerSingular <- FALSE
    lmerGood <- {
        {
            Nspectra >= 8
        } & {
            Npep >= 4
        } & {
            Nspectra/Npep > 4
        }
    }
    if (Nspectra == Npep)
        lmerGood <- FALSE  # don't do if one seq per spectrum
    if (GroupBy == "peptideId")
        lmerGood <- FALSE  # makes no sense if by peptide

    if (lmerGood) {
        coef_est <- rep(NA, numDataCols)
        secoef_est <- rep(NA, numDataCols)
        for (k in seq_len(numDataCols)) {
            # do for each fraction k=1
            y_k <- profileXX_i[, k]  # k'th column (fraction)

            varOK <- {
                var(y_k) > 0.001
            }
            if (varOK) {
            #result_k <- try(suppressMessages(lmer(y_k ~ 1 + (1 | pepId),
            #    data = profileXX_i)))
            lmerSingular <- FALSE

            result_k <- try(lmer(y_k ~ 1 + (1 | pepId),
                       data = profileXX_i))

            # Check for singular lmer fit;
            #  if so, set lmerSingular to FALSE and recalculate ALL means below


            if (isSingular(result_k)) {
              lmerSingular <- TRUE
              ## The following is for diagnostics only
              ##  error handling is as planned
              #cat("lmer matrix singular; recalculating average profile\n")
            }
            try(coef_est[k] <- fixef(result_k))
            try(secoef_est[k] <- sqrt(as.numeric(vcov(result_k))))
             }
            # variance too small for lmer:
            if (!varOK | is.na(coef_est[k])) {
                coef_est[k] <- mean(y_k)

                secoef_est[k] <- 0.001
            }
        }
    }

    if ({
        !lmerGood | lmerSingular
    } & {
        Nspectra == 1
    }) {
        coef_est <- as.numeric(profileXX_i[1, seq_len(numDataCols)])
        secoef_est <- rep(NA, numDataCols)
    }
    if ({
        !lmerGood | lmerSingular
    } & {
        Nspectra == 2
    }) {
        coef_est <- as.numeric(apply(profileXX_i[,
            seq_len(numDataCols)], 2, mean))
        secoef_est <- rep(NA, numDataCols)
    }
    if ({
        !lmerGood | lmerSingular
    } & {
        Nspectra >= 3
    }) {
        coef_est <- as.numeric(apply(profileXX_i[,
            seq_len(numDataCols)], 2, mean))
        secoef_est <- as.numeric(apply(profileXX_i[,
            seq_len(numDataCols)], 2, sd))/sqrt(Nspectra)
        coef_median_est <- as.numeric(apply(profileXX_i[,
            seq_len(numDataCols)], 2, median))
    }

    # apply(profile99.i,2,mean) # test; should be
    # similar to coef_est
    raw_coef_orig <- 2^coef_est - eps
    facAdj <- sum(raw_coef_orig)  # adjustment factor to make things add to 1
    raw_coef <- raw_coef_orig/facAdj
    coef_est <- log2(raw_coef + eps)

    # if GroupBy = 'protId', prot_i is the
    # protein number, and i is the same if
    # GroupBy = 'pepId', prot_i is the protein
    # number, and i is the peptide number
    result_i <- c(coef_est, secoef_est, Nspectra, Npep,
        protId_i, pepId_i)
    return(result_i)
}



# ================================================================
#' protProfileSummarize calculates mean and SE for each channel in
#' each prot using random effect model or arithmetic mean
#' random effect model can avoid dominance of a sequence with
#' too many spectra
#'
#' @param  protsCombineCnew a matrix of protein information and normalized
#'          specific amounts and outlier information
#' @param  numRefCols number of columns before Mass Spectrometry
#'          data columns
#' @param  numDataCols how many columns in MS data
#' @param  refColsKeep  which reference (non-data) columns to keep
#' @param  eps epsilon to avoid log(0)
#' @param  GroupBy 'protId' if average by protein; 'peptideId'
#'         if average by peptide
#' @param  outlierExclude 'none', 'spectra', or
#'         'spectraAndpeptide' (default)
#'          according to exclusion level
#' @param     cpus 1 (default);
#'            if cpus > 1 use BiocParallel with SnowParm(cpus)
#' @importFrom lme4 lmer
#' @importFrom BiocParallel bplapply
#' @export
#' @return protProfileSummaryIdentity: mean or weighted mean normalized
#'          specific amount profiles
#' @examples
#' set.seed(17356)
#' eps <- 0.029885209
#' data(TLN1_test)
#' flagSpectraBox <- outlierFind(protClass=TLN1_test,
#'            outlierLevel='peptide', numRefCols=5, numDataCols=9,
#'            outlierMeth='boxplot', range=3, eps=eps,
#'            randomError=TRUE, cpus=1)
#' # examine numbers of spectra that are outliers
#' table(flagSpectraBox$outlier.num.spectra)

#' pepProfiles <- profileSummarize(protsCombineCnew=flagSpectraBox,
#'            numRefCols=6, numDataCols=9, refColsKeep=c(1,2,4),eps=eps,
#'            GroupBy='peptideId', outlierExclude='spectra', cpus=1)
#' str(pepProfiles, strict.width='cut', width=65)
#'
# ================================================================
profileSummarize <- function(protsCombineCnew, numRefCols,
    numDataCols, refColsKeep = c(1, 2), eps, GroupBy = "protId",
    outlierExclude = "spectraAndPeptide", cpus=1) {

    # outlierExclude: none: don't exclude any
    # outlers spectra: (default) exclude only
    # spectra within peptides spectraAndPeptide:
    # exclude spectra-within-peptide outliers and
    # peptide outliers

    # # # # # # # # # # # # # # # # # # # # # # #
    # # # # # # # # # # # # # # # # # # # # # # #
    # # # # # # # # #
    if (GroupBy == "protId")
        uniqueLabel <- protsCombineCnew$protId
    if (GroupBy == "peptideId")
        uniqueLabel <- protsCombineCnew$pepId
    # if (GroupBy == 'protId') indList <-
    # unique(protsCombineCnew$protId) if (GroupBy
    # == 'peptideId') indList <-
    # unique(protsCombineCnew$pepId)
    indList <- unique(uniqueLabel)  # 28 Jan 2022

    if (!(outlierExclude %in% c("none", "spectra",
        "spectraAndPeptide"))) {
        cat("outlierExclude must be one of none, spectra, or peptide\n")
        stop()
    }
    if ((GroupBy == "peptideId") & (outlierExclude ==
        "spectraAndPeptide")) {
        cat("if GroupBy == 'peptideId' then outlierExclude
            cannot equal 'spectraAndPeptide' \n")
        stop()
    }

    # indList <- unique(uniqueLabel)
    n_prots <- length(indList)  # either proteins or proteins/peptides


  if (cpus == 1) {
    system.time(result <- lapply(indList, meansByProteins,
                uniqueLabel = uniqueLabel, protsCombineCnew = protsCombineCnew,
                numRefCols = numRefCols, numDataCols = numDataCols,
                GroupBy = GroupBy, eps = eps, outlierExclude = outlierExclude))
  }
  if (cpus > 1) {
    system.time(result <- bplapply(indList, meansByProteins,
               uniqueLabel = uniqueLabel, protsCombineCnew = protsCombineCnew,
               numRefCols = numRefCols, numDataCols = numDataCols,
               GroupBy = GroupBy, eps = eps, outlierExclude = outlierExclude,
               BPPARAM = SnowParam(cpus)))
  }

 # }
    # convert list of matrices to one matrix
    temp <- do.call(what = "rbind", result)
    temp_df <- data.frame(temp)

    names_channels <- names(protsCombineCnew)[numRefCols +
        c(seq_len(numDataCols))]
    names(temp_df)[seq_len(numDataCols)] <- names_channels
    names(temp_df)[(numDataCols + 1):(numDataCols +
        numDataCols)] <- paste(names_channels, ".se",
        sep = "")
    names(temp_df)[2 * numDataCols + seq_len(4)] <- c("Nspectra",
        "Npep", "protId", "pepId")

    if (GroupBy == "protId")
        refColsKeep <- 2  # only keep the 'prot' column
    # indices of second data set:
    index <- match(temp_df$pepId, protsCombineCnew$pepId)
    protsMini <- protsCombineCnew[index, refColsKeep]

    if (GroupBy == "protId") {

        ncolTemp <- ncol(temp_df)
        # remove protId and pepId (last two columns):
        temp_df <- temp_df[, -c(ncolTemp - 1, ncolTemp)]
        protProfileSummaryRes <- data.frame(protsMini,
            temp_df)  # 'protsMini' is protein names
        names(protProfileSummaryRes)[1] <- "prot"
    }
    if (GroupBy == "peptideId") {

        protProfileSummaryRes <- data.frame(protsMini,
            temp_df)
    }



    # for 'Identity' (untransformed) results,
    # drop the standard error terms, which don't
    # apply
    if (GroupBy == "protId") {
        # note that first two columns are protein
        # and peptide names
        protProfileSummaryIdentity <- protProfileSummaryRes[,
            -((2 + numDataCols):(1 + 2 * numDataCols))]
        # now reverse the log2 transformation
        protProfileSummaryIdentity[, 2:(1 + numDataCols)] <-
            2^protProfileSummaryRes[,
            c(2:(1 + numDataCols))] - eps
        # now eliminate negative values due to
        # roundoff error
        negVals <- {
            {
                protProfileSummaryIdentity[, 2:(1 +
                  numDataCols)] < 0
            } & {
                protProfileSummaryIdentity[, 2:(1 +
                  numDataCols)] > -1e-10
            }
        }
        protProfileSummaryIdentity[, 2:(1 + numDataCols)][negVals] <- 0
    }
    if (GroupBy == "peptideId") {
        nKeep <- length(refColsKeep)
        # note that first column is protein name.
        # (There is no peptide name here)
        # protProfileSummaryIdentity <-
        # protProfileSummaryRes[,-((2 +
        # numDataCols):(1 + 2*numDataCols))]

        # drop the standard errors of the mean
        # profiles (9 columns)
        protProfileSummaryIdentity <- protProfileSummaryRes[,
            -((nKeep + 1 + numDataCols):(nKeep + 2 *
                numDataCols))]
        # now reverse the log2 transformation
        protProfileSummaryIdentity[, (nKeep + 1):(nKeep +
            numDataCols)] <- 2^protProfileSummaryRes[,
            (nKeep + 1):(nKeep + numDataCols)] - eps
        # now eliminate negative values due to
        # roundoff error
        negVals <- {
            {
                protProfileSummaryIdentity[, (nKeep +
                  1):(nKeep + numDataCols)] < 0
            } & {
                protProfileSummaryIdentity[, (nKeep +
                  1):(nKeep + numDataCols)] > -1e-10
            }
        }
        protProfileSummaryIdentity[, (nKeep + 1):(nKeep +
            numDataCols)][negVals] <- 0
    }

    protProfileSummaryIdentity
}
