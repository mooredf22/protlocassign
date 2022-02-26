# protlocassign: A method for proportionally assigning protein residence among subcellular organelles

## Introduction
This package implements constrained proportional assignment to assign proteins proportionately to their subcellular organelle residences. It uses quantitative mass spectrometry data from subcellular fractionation experiments to create profiles describing the distribution of proteins and peptides across centrifugation fractions. The package identifies outlier spectra and computes weighted means across peptides using random effects models. Then, using reference proteins that reside in single compartments, it estimates the distribution of each protein across these compartments by obtaining optimal mixtures of the reference profiles. 

## Installation
Start by installing the devtools package from CRAN, by typing

```
install.packages("devtools")
```

Then install the protlocassign package by typing

```
devtools::install_github("mooredf22/protlocassign")
library(protlocassign)
```

This will make the programs and datasets available. Instructions for use of the package are in the vignette.

## References

Jadot, M.; Boonen, M.; Thirion, J.; Wang, N.; Xing, J.; Zhao, C.; Tannous, A.; Qian, M.; Zheng, H.; Everett, J. K., Accounting for protein subcellular localization: A compartmental map of the rat liver proteome. Molecular & Cellular Proteomics 2017, 16, (2), 194-212.

Tannous, A.; Boonen, M.; Zheng, H.; Zhao, C.; Germain, C. J.; Moore, D. F.; Sleat, D. E.; Jadot, M.; Lobel, P., Comparative Analysis of Quantitative Mass Spectrometric Methods for Subcellular Proteomics. J Proteome Res 2020, 19, (4), 1718-1730.


