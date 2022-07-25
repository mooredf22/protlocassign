\docType{package}
\name{protlocassign}
\alias{protlocassign}
\title{Constrained proportional assignment of proteins}
\description{
The `protlocassign` package implements constrained proportional assignment
to assign proteins proportionately to their subcellular organelle 
residences. It uses quantitative mass spectrometry data from 
subcellular fractionation experiments to create profiles 
describing the distribution of proteins and peptides across 
centrifugation fractions. The package identifies outlier 
spectra and computes weighted means across peptides using 
random effects models. Then, using reference proteins that 
reside in single compartments, it estimates the distribution 
of each protein across these compartments by obtaining optimal 
mixtures of the reference profiles.
}

\references{
Moore DF, Sleat DE, Lobel P A method to estimate the distribution 
  of proteins across
  multiple compartments using data from quantitative proteomics subcellular 
  fractionation experiments.
  Journal of Proteome Research, 2022, 21, 6, 1371-1381 
  https://doi.org/10.1021/acs.jproteome.1c00781
  
Tannous, A.; Boonen, M.; Zheng, H.; Zhao, C.;
 Germain, C. J.; Moore, D. F.;
 Sleat, D. E.; Jadot, M.; Lobel, P., Comparative Analysis of
 Quantitative Mass Spectrometric Methods for Subcellular Proteomics.
 J Proteome Res 2020, 19, (4), 1718-1730.
}

