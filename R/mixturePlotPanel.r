#' Panel of mixture plots of all pairwise simulated combinations
#' 
#' For all sets of simulated proteins distributed between two compartments, 
#'   plot the CPA estimates for a given type of profile versus true 
#'   mixture proportions.; aAssumes eight or fewer compartments 
#'   for proper formatting of the plots. 
#' 
#' @param refLocationProfilesAcup data frame containing the Acup 
#'        profiles for the reference compartments
#' @param totProt vector of total protein amounts 
#'        (derived from a given amount of starting material) 
#'        in each of the fractions comprising the profile 
#' @param NstartMaterialFractions Number of fractions that reconstitute 
#'       the starting material, e.g., a complete set of differential 
#'       centrifugation fractions.  For experiment AT5, it is 6 
#'       ( N, M, L1, L2, P, and S).
#' @param errorReturn return all area-based errors if true
#' @param fitType use RSA, NSA, or Acup
#' @param log2Transf use log2-transformed values. Default is "FALSE"
#' @param eps small value to add so that log argument is greater than zero
#' @importFrom graphics par
#' @importFrom graphics text
#' @importFrom graphics legend
#' @importFrom graphics layout
#' @export
#' @return Panel of plots and a list of errors for all mixtures
#' @examples
#' # See Tutorial 4 for details.

mixturePlotPanel <- function(refLocationProfilesAcup, totProt,
                             NstartMaterialFractions,
                             errorReturn=FALSE,
                             fitType, log2Transf=FALSE, eps=0.001) {
  #windows(width=8,height=10)
  numDataCols <- ncol(refLocationProfilesAcup)
  numCompart <- nrow(refLocationProfilesAcup)
  numPlots <- numCompart*(numCompart - 1)/2   # numCompart choose two
  refLocationProfilesRSA <- RSAfromAcup(Acup=refLocationProfilesAcup,
                               NstartMaterialFractions=NstartMaterialFractions,
                               totProt=totProt)
  refLocationProfilesNSA <- NSAfromRSA(refLocationProfilesRSA)

  if (numCompart > 8) {
    warning("Error: too many compartments to plot on one page\n")
    warning("Error: too many compartments to plot on one page\n")
  }

  # case where there are 7 or 8 compartments:
  if(numCompart == 8) {
  layout(rbind(c(3,1,1,1,1),
               c(3,4,5,6,7),
               c(3,8,9,10,11),
               c(3,12,13,14,15),
               c(3,16,17,18,19),
               c(3,20,21,22,23),
               c(3,20,21,22,23),
               c(3,24,25,26,27),
               c(3,24,25,26,27),
               c(3,28,29,30,31),
               c(3,2,2,2,2)),

         heights=c(1.10,2,2,2,2,2,2,2, 0.35),
         widths=c(0.4,2,2,2,2),respect=FALSE)
  }

  # case where there are 7 compartments
  if(numCompart == 7) {
    layout(rbind(c(3,1,1,1,1),
                 c(3,4,5,6,7),
                 c(3,8,9,10,11),
                 c(3,12,13,14,15),
                 c(3,16,17,18,19),
                 c(3,20,21,22,23),
                 c(3,24,25,26,27),
                 c(3,2,2,2,2)),

           heights=c(1.10,2,2,2,2,2,2, 0.35),
           widths=c(0.4,2,2,2,2),respect=FALSE)
  }

  # case where there are 6 compartments
  if(numCompart == 6) {
    layout(rbind(c(3,1,1,1,1),
                 c(3,4,5,6,7),
                 c(3,8,9,10,11),
                 c(3,12,13,14,15),
                 c(3,16,17,18,19),
                 c(3,2,2,2,2)),

           heights=c(1.10,2,2,2,2, .35),
           widths=c(0.4,2,2,2,2),respect=FALSE)
  }

  # case where there are 5 compartments (or fewer)
  if(numCompart == 5) {
    layout(rbind(c(3,1,1,1,1),
                 c(3,4,5,6,7),
                 c(3,8,9,10,11),
                 c(3,12,13,14,15),
                 c(3,2,2,2,2)),

           heights=c(1.10,2,2,2, .35),
           widths=c(0.4,2,2,2,2),respect=FALSE)
  }

  #layout.show(31)
  # this program assumes exactly eight subcellular compartments
  # set up color and point lists
  loc.list <- rownames(refLocationProfilesAcup)


  n.loc <- length(loc.list)
  col.list <- c("red", "blue", "orange", "darkgreen", "orange", "lightblue",
                "purple", "green")
  pch.list <- c(1, 2, 3, 4, 17, 6, 15, 8)
  col.list <- col.list[seq_len(n.loc)]
  pch.list <- pch.list[seq_len(n.loc)]

  # create header (region 1)
  x <- c(0,5)
  y <- c(0,1.1)
  par(mar=c(0,0,0,0))
  plot(y ~ x,type="n",axes=FALSE, cex=1.9, xlab="", ylab="")
  #aa <- 1.34
  #log2Transf <- F
  #log2Transf <- T
  #fitType <- "Acup"
  #fitType <- "original"

  #if (!log2Transf) transF <- "no transformation"
  #if (log2Transf) transF <- "log2 transformed"
  if (!log2Transf) transF <- ""
  if (log2Transf) transF <- "log2"

  text(x=2.5,y=0.8,paste("Synthetic Protein CPAs,", transF, fitType), cex=1.9)
  legend(x="bottom", legend=loc.list, col=col.list, pch=pch.list, ncol=8)

  # create x-label at bottom of plot (region 2)
  x <- c(0,5)
  y <- c(0,0.5)
  par(mar=c(0,0,0,0))
  plot(y ~ x,type="n",axes=FALSE, cex=1, xlab="", ylab="")
  #aa <- 1.34

  text(x=2.5,y=0.3,"True proportion", cex=1.7, xlab="", ylab="")

  # create y-label on left side of plot (region 3)
  x <- c(0,5)
  y <- c(0,0.5)
  par(mar=c(0,0,0,0))
  plot(y ~ x,type="n",axes=FALSE,cex=1.7, xlab="", ylab="")
  #aa <- 1.34

  text(x=2.5,y=0.3,"Estimated proportion", cex=1.5, srt=90)
  # mixtures must be a list of equally spaced proportions

  plotLables <- data.frame(loc.list, col.list, pch.list)
  par(mar=c(3,3,3,3))
  kk <- 0
  mixErrorMat <- NULL
  for (i in seq_len((numCompart-1))) {
    for (j in (i+1):numCompart) {
      #mixturePlot(i=i, j=j, type=fitType, log2Transf=log2Transf)
      mixProtiProtjAcup <- proteinMix(refLocationProfilesAcup, Loc1=i, Loc2=j)

      mixProtiProtjRSA <- RSAfromAcup(Acup=mixProtiProtjAcup,
                            NstartMaterialFractions=NstartMaterialFractions,
                            totProt=totProt)
      if (!log2Transf & {fitType == "RSA"}) {
       mixProtiProtjCPA <- fitCPA(profile=mixProtiProtjRSA,
                                 refLocationProfiles=refLocationProfilesRSA,
                                  numDataCols=numDataCols)
      }

      if (log2Transf& {fitType == "RSA"}) {
      # Take a log2 transformation
        log2MixProtiProtjRSA <- log2(mixProtiProtjRSA + eps)
        log2refLocationProfilesRSA <- log2(refLocationProfilesRSA + eps)
        mixProtiProtjCPA <- fitCPA(profile=log2MixProtiProtjRSA,
                               refLocationProfiles=log2refLocationProfilesRSA,
                               numDataCols=numDataCols)
      }
      if (!log2Transf & {fitType == "NSA"}) {
        mixProtiProtjNSA <- NSAfromRSA(mixProtiProtjRSA)
        mixProtiProtjCPA <- fitCPA(profile=mixProtiProtjNSA,
                                 refLocationProfiles=refLocationProfilesNSA,
                                 numDataCols=numDataCols)
      }

      if (log2Transf& {fitType == "NSA"}) {
        # Now convert to normalized specific amounts
        mixProtiProtjNSA <- NSAfromRSA(mixProtiProtjRSA)

        # Take a log2 transformation
        log2MixProtiProtjNSA <- log2(mixProtiProtjNSA + eps)
        log2refLocationProfilesNSA <- log2(refLocationProfilesNSA + eps)


        mixProtiProtjCPA <- fitCPA(profile=log2MixProtiProtjNSA,
                              refLocationProfiles=log2refLocationProfilesNSA,
                              numDataCols=numDataCols)
      }

      if (!log2Transf & {fitType == "Acup"}) {
        mixProtiProtjCPA <- fitCPA(profile=mixProtiProtjAcup,
                                 refLocationProfiles=refLocationProfilesAcup,
                                 numDataCols=numDataCols)
      }

      if (log2Transf& {fitType == "Acup"}) {


        mixProtiProtjCPA <- fitCPA(profile=log2(mixProtiProtjAcup + eps),
                      refLocationProfiles=log2(refLocationProfilesAcup + eps),
                      numDataCols=numDataCols)
      }

      mixResult <- mixturePlot(mixProtiProtjCPA=mixProtiProtjCPA,
                           NstartMaterialFractions=NstartMaterialFractions,
                           Loc1=i, Loc2=j,
                           increment.prop=0.1, errorReturn=errorReturn)
      if (errorReturn) {
        mixErrorMat <- rbind(mixErrorMat, mixResult)
        mixErrorMat
      }
    }
  }


  if (errorReturn) return(mixErrorMat)
}
