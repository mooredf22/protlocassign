\name{totProtAT5}
\alias{totProtAT5}
\docType{data}
\title{Amount of total protein in each fraction 
}
\value{amount of total protein in each fraction}
\description{
This is a vector that gives the amount of total protein derived 
from a gram of tissue homogenate in the nine fractions of 
experiment AT5 from Tannous et al. 2000. 
The values are in the same order as those in protNSA_test: 
there is a value for each of the six differential centrifugation 
fractions (N, M, L1, L2, P and S) and the three Nycodenz fractions. 
Note that the sum of the six differential centrifugation fractions 
is used to estimate the total protein in the tissue homogenate.  
}
\usage{data(totProtAT5)}
\format{
 A vector containing total protein values for fractions 
      “N” “M” “L1” “L2” “P” “S” “Nyc.1” “Nyc.2” “Nyc.3””
}

\references{
Tannous, A.; Boonen, M.; Zheng, H.; Zhao, C.; Germain, C. J.; Moore, D. F.; 
Sleat, D. E.; Jadot, M.; Lobel, P., Comparative Analysis 
of Quantitative Mass Spectrometric Methods for Subcellular Proteomics. 
J Proteome Res 2020, 19, (4), 1718-1730.
}
\examples{
data(totProtAT5)
totProtAT5
}
\keyword{datasets}
