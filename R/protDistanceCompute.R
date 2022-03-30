
#' Compute distances between proteins
#' 
#' Compute distances between a particular protein or profile
#'  and all other proteins,
#'  and list the nearest ones
#'
#' @param protName  Name of protein to which distances are to be computed
#' @param n.nearest Number of nearest proteins to list
#' @param distProts distance matrix created by, for example
#' @param protNames A list of all proteins in a dataset
#' @param profile dataframe of profiles for proteins
#' @return List of the proteins in protName closest to protName or to profile
#' @export
#' @examples
#' data(protNSA_test)
#' distUse <- dist(protNSA_test[,seq_len(9)], method='euclidean')
#' protsUse <- names(protNSA_test)
#' nearestProts(protName='CTSD', n.nearest=10,  distProts=distUse,
#'   protNames=protsUse, profile=protNSA_test[,seq_len(9)])
nearestProts <- function(protName, n.nearest = 5, distProts,
    protNames, profile) {
    distProtsMat <- as.matrix(distProts)

    ref <- protIndex(protName, profile)
    if (is.character(ref)) {
        return(ref)
    }
    if (nrow(ref) > 1) {
        warning("More than one protein matches protName\n")
        return(ref)
    }

    ind.ref <- ref[1, 1]
    # vector of distances to the reference protein
    vect.dist <- distProtsMat[ind.ref, ]



    nearest.list <- sort(vect.dist)

    resultAll <- data.frame(names(nearest.list), as.numeric(nearest.list))
    names(resultAll) <- c("protName", "euclidean distance")
    result <- resultAll[seq_len(n.nearest), ]
    result

}

# nearestProts('AADAC', n.nearest=10,
# distProts=distUse, protNames=protsUse)

# nearestProts('AAD', n.nearest=10,
# distProts=distUse, protNames=genesUse)

