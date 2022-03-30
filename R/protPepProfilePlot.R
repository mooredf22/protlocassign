#' Plot profile and peptide profiles with reference profiles
#' 
#' This function plots the average profiles of any protein in the dataset,
#'   the peptide profiles, and also the reference profile for each compartment
#'
#' @param protName Name of the protein to plot
#' @param protProfile protein profile
#' @param Nspectra indicator for if there are columns in profile for Nspectra
#'         (number of spectra) and Npep (number of peptides)
#' @param pepProfile peptide profiles
#' @param numRefCols number of reference columns
#'        (preceding the data profile columns)
#' @param numDataCols  number of fractions per protein
#' @param n.compartments number of compartments (8 in Jadot data)
#' @param refLocationProfiles A matrix refLocationProfiles giving the
#'          abundance level profiles of the subcellular locations
#'        n.compartments = 8 columns are subcellular locations, and
#'        numDataCols rows are the fraction names
#' @param assignPropsMat A matrix of assignment proportions, from the
#'          constrained proportional assignment algorithm,
#'        and optionally upper and lower 95 percent confidence limits
#' @param propCI True if lower and upper confidence intervals
#'            are included in assignPros
#' @param transType label for transformation; default is ""
#' @param yAxisLabel label for y-axis
#' @return plot of average, peptide, and reference profiles
#' @examples
#'   # See Vignette 6 for a full explanation
#' @importFrom graphics layout
#' @export

protPepPlotfun <- function(protName, protProfile, Nspectra=TRUE,
                   pepProfile=NULL, numRefCols, numDataCols, n.compartments=8,
                   refLocationProfiles, assignPropsMat, propCI=FALSE,
                   transType="", yAxisLabel="") {


  oldpar <- par(no.readonly=TRUE)
  #on.exit(par(oldpar))
  protsOK <- {rownames(protProfile) == rownames(assignPropsMat)}

  temp <- protIndex(protName, protProfile, exactMatch=TRUE)
  protName.i <- protName
  # must be exact to avoid duplicate finds

  # this can be a vector, matrix, or vector, so handling is complicated
  #if (is.matrix(tempx)) temp <- tempx   # leave it alone
  # if temp is not a matrix, can then test for being NA with no error returned
  if (!is.matrix(temp)) {
    if (is.na(temp)[1]) {
      warning(paste(protName, " not found \n"))  # first element is Na
      return(temp)
    }
  }
  if (nrow(temp) > 1) {
    warning(paste("more than one protein matches pattern \n"))
    return(temp)
  }
  protPlot <- temp[1,1]   # works even if temp is a vector
  # protPlot is the index number of the protein

  # just the protProfiles
  meanProteinLevels <- protProfile[,seq_len(numDataCols)]

  protNames <- rownames(protProfile)
  protNameUnique <- make.unique(protNames)
  rownames(meanProteinLevels) <- protNameUnique
  subCellNames <- rownames(refLocationProfiles)

  fractions.list <- colnames(refLocationProfiles)  # column names[protPlot])

  # # # # # # # # # # # # # #
  #  Do the following if "finalList" (the full list of peptides and spectra)
  #   are available
  # # # # # # # # # # # # # #

  #if (!is.null(finalList)) {


  # # # # # # # # # # # # # #
  #  Do the following if "peptideList" (the list of mean peptide profiles)
  #   is available
  # # # # # # # # # # # # # #
  if (!is.null(pepProfile)) {

    pepProfile.i <- pepProfile[pepProfile$prot == protName.i,]
    means.peptides.i <- pepProfile.i[,numRefCols + seq_len(numDataCols)]
    n.spectra.i <- pepProfile.i$Nspectra
    outlierFlagVec.i <- {pepProfile.i$outlier.num.peptides > 0}
    n.unique.peptide.i <- nrow(pepProfile.i)
  }


  yy <- as.numeric(meanProteinLevels[protPlot,])

  xvals <- seq_len(numDataCols)
  #max.y <- max(channelsAll, na.rm=TRUE)
  #max.y <- max(meanProteinLevels, na.rm=TRUE)
  min.y <-0
  #if (dataType == "raw") min.y <- 0
  loc.list <- subCellNames
  #windows(width=10, height=8)
  #par(mfrow=c(3,3))
  layout(rbind(c(14,1,1,1,15),
               c(14,2,3,4,15),
               c(14,5,5,5,15),
               c(14,6,7,8,15),
               c(14,9,9,9,15),
               c(14,10,11,12,15),
               c(14,13,13,13,15)),
         #c(16,12,13,14,15,17),
         #c(16,18,18,18,18,17)),
         heights=c(0.75,2,0.25,2,0.25,2,0.25),
         widths=c(0.4,2,2,2,0.4),respect=FALSE)



  #layout.show(15)
  x <- c(0,5)
  y <- c(0,0.5)
  par(mar=c(0,0,0,0))
  plot(y ~ x,type="n",axes=FALSE,cex=1)
  #aa <- 1.34

  #if (length(indAssignProp.prot) > 0) {
    text(x=2.5,y=0.3,paste(protName.i), cex=2)

    NpeptidesPlot <- protProfile$Npep[protPlot]
    NspectraPlot <- protProfile$Nspectra[protPlot]
    NpeptidesPlotText <- " peptides and "

    #browser()

    if (NpeptidesPlot == 1) NpeptidesPlotText <- " peptide and "
    NspectraPlotText <- " spectra"
    if (NspectraPlot == 1) NspectraPlotText <- " spectrum"
    text(x=2.5,y=0.1, paste(NpeptidesPlot, NpeptidesPlotText,
                            NspectraPlot , NspectraPlotText), cex=2)
  #}
 # if (length(indAssignProp.prot) ==0) {
 #   text(x=2.5,y=0.3,protName.i, cex=2)
 # }

  # max.y <- max(c(max(means.peptides.i), max(refLocationProfiles[i,])))
  min.y <- 0
  par(mar=c(2,4,2,1.5))
  # The plots are in alphabetical order
  # re-arrange the plots so that they are in this order
  #   Mito (7)     Lyso  (4)   Perox (5)
  #   ER (2)       Golgi (1)   PM (8)
  #   Cyto (3 )    Nuc (6)
  #
  #loc.ord <- c(7, 4, 5, 2, 1, 8, 3, 6)
  loc.ord <- c(5, 4, 7, 2, 3, 8, 1, 6)
  for (i in loc.ord) {    # do all the subcellular locations
    # i=5
    if (TRUE) {
      #if ({loc.ord[i] == 4} | {loc.ord[i] == 7}) {
      if ({i == 2} | {i == 1}) {
        x <- c(0,5)
        y <- c(0,0.5)
        par(mar=c(0,1,0,0))
        plot(y ~ x,type="n",axes=FALSE,cex=1, ylab="")
        par(mar=c(2,3.5,2,1.5))
      }
    }

    assign.i <- names(meanProteinLevels)[i]  # channel i name

    assignLong.i <- subCellNames[i]

    #mean.i <- markerLoc[,i]
    mean.i <- as.numeric(refLocationProfiles[i,])
    max.y <- max(c(max(means.peptides.i, na.rm=TRUE),
                   max(refLocationProfiles[i,])))
    n.uniq.peptide.i <- nrow(means.peptides.i)
    par(mar=c(2,4.1,2,1.5))
    plot(mean.i ~ xvals,  axes=FALSE, type="l",
         ylim=c(min.y, max.y), ylab=transType)
    axis(1,at=xvals,labels=fractions.list, las=2)
    if (FALSE) {
    axis(side=1,at=xvals,labels=FALSE, cex.axis=0.6)
    text(x = seq_len(length(fractions.list)),
         ## Move labels to just below bottom of chart.
         y = par("usr")[3] - max(mean.i)/12,
         ## Use names from the data list.
         labels = fractions.list,
         ## Change the clipping region.
         xpd = NA,
         ## Rotate the labels by 45 degrees.
         srt = 45,
         ## Adjust the labels to almost 100% right-justified.
         adj = 0.965,
         ## label size.
         cex = 0.8)
    }
    axis(2)

  if (TRUE) {
    for (j in seq_len(n.uniq.peptide.i)) {
      lwdplot <- 1
      colplot <- "cyan"
      if (n.spectra.i[j] > 1) {
        lwdplot <- 2
        colplot <- "deepskyblue"
      }
      if (n.spectra.i[j] > 2) {
        lwdplot <- 3
        colplot <- "dodgerblue3"
      }
      if (n.spectra.i[j] > 5)   {
        lwdplot <- 4
        colplot <- "blue"
      }

      if (outlierFlagVec.i[j] != 1) {
        lines(as.numeric(means.peptides.i[j,]) ~ xvals,
              cex=0.5, lwd=lwdplot, col=colplot)
      }
      if (outlierFlagVec.i[j] == 1) {
        lines(as.numeric(means.peptides.i[j,]) ~ xvals,
                                          cex=0.5, lwd=2, col="orange")
      }

    }
  }




    lines(yy ~ xvals, col="red", lwd=2)
    lines(mean.i ~ xvals, lwd=4, col="black")
    lines(mean.i ~ xvals, lwd=2, col="yellow", lty=2)
    title(paste(assignLong.i, "\n p = ",
                round(assignPropsMat[protPlot,i], digits=2 )))
  }
  x <- c(0,5)
  y <- c(0,0.5)
  par(mar=c(0,0,0,0))
  plot(y ~ x, type="n", axes=FALSE)
  if (TRUE) {
    legend(x=1, y=0.4, legend=c("Reference profile",
                "Average profile", "1 spectrum", "2 spectra", "3-5 spectra",
                "6+ spectra", "Outliers"),
         col=c("black", "red", "cyan", "deepskyblue",
               "dodgerblue3", "blue", "orange"),
         lwd=c(4,2,1,2,3,4,2), lty=c(1,1,1,1,1,1,1))
     legend(x=1, y=0.4, legend=c("Reference profile", "Average profile",
        "1 spectrum", "2 spectra", "3-5 spectra", "6+ spectra", "Outliers"),
         col=c("yellow", "red", "cyan", "deepskyblue",
               "dodgerblue3", "blue", "orange"),
         lwd=c(2,2,1,2,3,4,2), lty=c(2,1,1,1,1,1,7))
  }
  if (FALSE) {
    legend(x=1, y=0.4, legend=c("Reference profile", "Average profile"),
           col=c("yellow", "red"), lwd=c(2,2), lty=c(1,1))
    legend(x=1, y=0.4, legend=c("Reference profile", "Average profile"),
           col=c("black", "red"), lwd=c(2,2), lty=c(2,1))
  }
  x <- c(0,5)
  y <- c(0,0.5)

  plot(y ~ x, type="n", axes=FALSE)
  par(mar=c(0,0,0,0))
  plot(y ~ x, type="n", axes=FALSE)


  text(x=2.5, y=0.25, labels=yAxisLabel, srt=90, cex=2 )
}
