\name{totProtAT5}
\alias{totProtAT5}
\docType{data}
\title{totProtAT5
%%   ~~ data name/kind ... ~~
}
\description{
Adjustment factors for the tmtMS2 data, to adjust the relative protein
amounts as measured
by the mass spectrometer to the amounts of protein in each channel in
the original sample.
This is a vector with nine values: there is one value for each of the
six fraction channels and the three Nycodenz channels.
Note that the sum of the three Nycodenz channels
(Nyc.1, Nyc.2, Nyc.3) equals the amount in L1.
These values are used to adjust the levels in the data set
"geneProfileSummaryTMTms2" to re-create the levels
in the original samples.
}
\usage{data("totProtAT5")}
\format{
  The format is:
 names: "N" "M" "L1" "L2" "P" "S" "Nyc.1" "Nyc.2" "Nyc.3""
}

\references{
A Tannous et al. (2020)
%%  ~~ possibly secondary sources and usages ~~
}
\examples{
data(totProtAT5)
totProtAT5
}
\keyword{datasets}
