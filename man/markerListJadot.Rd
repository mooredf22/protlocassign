\name{markerListJadot}
\alias{markerListJadot}
\docType{data}
\title{Reference proteins from Jadot et al. (2016)

}
\description{
These are the reference proteins curated by Jadot et al. (2016)
for each of the subcellular locations Cytosol, ER, Golgi, Lysosome,
Mitochondria, Nucleus, Peroxisomes, and Plasma Membrane.
The protein "TPP1" was originally used as a reference for Lysosomes,
but has been removed here.

}
\format{
  A data frame with 39 observations on the following 2 variables.
  \describe{
    \item{\code{protName}}{a factor with levels \code{ACP2}
     \code{ATP1A1} \code{B4GALT1} \code{CANX}
    \code{CAT} \code{CD38} \code{CPS1} \code{CTSD} \code{DPP3}
    \code{ENPP1} \code{GANAB} \code{GLB1} \code{GLUD1} \code{GPT}
    \code{HAO1} \code{HEXA} \code{HIST1H1D} \code{ITGB1} \code{LMNA}
    \code{LMNB1} \code{LMNB2} \code{MAN1A2} \code{MAN2A1} \code{MDH2}
    \code{MGAT1}  \code{NT5E} \code{PC} \code{PCK1}
    \code{PGM1} \code{PHYH} \code{POR} \code{RPN2} \code{SDHA}
    \code{ST6GAL1} \code{UGGT1} \code{UGT2B37} \code{UOX}}
    \item{\code{referenceCompartment}}{a factor with levels
    \code{Cytosol} \code{ER} \code{Golgi} \code{Lysosome}
    \code{Mitochondria} \code{Nucleus} \code{Peroxisomes}
    \code{Plasma membrane}}
  }
}

\source{
Jadot M, Boonen M, Thirion J, Wang N, Xing J, Zhao C, Tannous A,
Qian M, Zheng H, Everett JK, Moore DF, Sleat DE, and Lobel P (2016)
Accounting for protein subcellular localization: a compartmental
map of the rat liver proteome.
Molecular and Cellular Proteomics 16, 194-212.
doi:10.1074/mcp.M116.064527 PMCID: PMC5294208

}

\examples{
data(markerListJadot)
str(markerListJadot)
}
\keyword{datasets}
