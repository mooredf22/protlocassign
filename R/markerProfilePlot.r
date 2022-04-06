

#' Plot profiles of reference proteins
#'
#' This function plots profiles of reference proteins
#'    and also the average profile for each compartment
#'
#' @param refLoc  name of the reference subcellular
#'          compartment to plot
#' @param markerList list of reference proteins
#' @param profile data frame of specified protein (row name) profiles
#' @param ylab label for y-axis of plot, e.g., NSA or RSA or Acup
#' @param refLocationProfiles data frame of profiles for the 
#'       reference compartments
#' @param ylab label for y=axis, eg, 'NSA'
#' @param refProtPlot indices of reference proteins to plot; default is NULL
#' @importFrom graphics par
#' @importFrom graphics axis
#' @importFrom graphics text
#' @importFrom graphics lines
#' @importFrom graphics title
#' @export
#' @return Plot of profiles of reference (marker) proteins
#' @examples
#' data(protNSA_test)
#' data(markerListJadot)
#' data(refLocProfNSA)
#' markerProfilePlot(refLoc='PM', profile=protNSA_test,
#'    markerList=markerListJadot,
#'    refLocationProfiles=refLocProfNSA, ylab='NSA')


markerProfilePlot <- function(refLoc, profile, markerList,
    refLocationProfiles, ylab = "", refProtPlot = NULL) {
    n.channels <- ncol(refLocationProfiles)
    protNames <- rownames(profile)
    # add protName column for merge
    profileWithProteins <- data.frame(protNames, profile)
    names(profileWithProteins)[1] <- "protName"
    meanReferenceProts <- merge(x = markerList, y = profileWithProteins,
        by.x = "protName", by.y = "protName", all.x = FALSE,
        sort = FALSE)
    markerLoc <- t(refLocationProfiles)
    max.val <- max(markerLoc)

    fractions.list <- rownames(markerLoc)
    compartments.list <- rownames(refLocationProfiles)

    # stop if refLoc not found
    if (!(refLoc %in% compartments.list)) {
        warning("reference compartment not found\n")
        return(NULL)
    }

    location.list <- colnames(markerLoc)

    n.loc <- ncol(markerLoc)

    # i=1 loc.i <- location.list[i]
    loc.i <- refLoc

    # channels.i <-
    # meanReferenceProts[meanReferenceProts$referenceCompartment
    # == loc.i,2+1:n.channels]
    channels.i <- meanReferenceProts[meanReferenceProts$referenceCompartment ==
        loc.i, 2 + seq_len(n.channels)]
    refprotsvec.i <-
     as.character(meanReferenceProts[meanReferenceProts$referenceCompartment ==
        loc.i, 1])
    n.refprots.i <- length(refprotsvec.i)

    mean.i <- as.numeric(markerLoc[, {
        loc.i == location.list
    }])
    xvals <- seq_len(length(mean.i))

    max.val <- max(channels.i)
    # ylim=c(0,max.val),
    par(mar = c(6.1, 4.1, 4.1, 4.1))  # change the margins
    plot(mean.i ~ xvals, ylim = c(0, max.val), axes = "F",
        type = "n", ylab = ylab, xlab = "")

    axis(1, at = xvals, labels = fractions.list, cex.axis = 1,
        las = 2)

    axis(2, las = 1)
    lines(mean.i ~ xvals, lwd = 5, lty = 1, col = "black")
    lines(mean.i ~ xvals, lwd = 2, lty = 2, col = "yellow")
    for (j in seq_len(nrow(channels.i))) {

        means.j <- as.numeric(channels.i[j, ])

        if (is.null(refProtPlot)) {
            lines(as.numeric(means.j) ~ xvals, col = "red")
        }
        if (!is.null(refProtPlot)) {
            if (refProtPlot == j) {
                lines(as.numeric(means.j) ~ xvals,
                  col = "red")
                refprot.j <- refprotsvec.i[j]
            }
        }
    }


    if (is.null(refProtPlot)) {
        title(main = paste(loc.i, "profiles\n", as.character(n.refprots.i),
            " reference proteins"))
    }
    if (!is.null(refProtPlot)) {
        title(main = paste(loc.i, "profiles\n", "reference protein ",
            as.character(refprot.j)))
    }

    # }
    return(NULL)

}

