
#' Identifies proteins with profiles similar to that of a specified protein.
#' 
#' Uses a distance matrix to identify proteins with profiles nearest 
#'    to that of a specified protein; can be used to identify proteins 
#'    with similar subcellular locations.
#'
#' @param protName  name of protein to which distances are to be computed
#' @param n.nearest number of nearest proteins to list
#' @param distProts distance matrix created by, for example, using the 'dist'
#'     function in R
#' @param protNames list of all proteins in a dataset
#' @param profile data frame of profiles for proteins
#' @return List of the proteins in data set with profiles closest to 
#'       that of specified protein (protName)
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


