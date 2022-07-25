# protlocassign: A method for proportionally assigning protein residence among subcellular organelles

## Introduction
This package implements constrained proportional assignment to 
assign proteins proportionately to their subcellular organelle residences. 
It uses quantitative mass spectrometry data from subcellular fractionation 
experiments to create profiles describing the distribution of proteins 
and peptides across centrifugation fractions. The package identifies 
outlier spectra and computes weighted means across peptides using 
random effects models. Then, using reference proteins that reside in 
single compartments, it estimates the distribution of each protein across 
these compartments by obtaining optimal mixtures of the reference profiles. 

## Installation
Start by installing the package from Bioconductor. For this you will
need the Bioconductor installation facility `BiocManager`, which may
be installed from CRAN as follows:

```
install.packages("BiocManager")
```

Then obtain the `protlocassign` package as follows:

```
BiocManager::install("protlocassign")
```

Alternatively, one may install a development version directly from
Github by first installing the devtools package from CRAN:

```
install.packages("devtools")
```

and then installing the github version of `protlocassign`:

```
devtools::install_github("mooredf22/protlocassign")
```

Next, install other required packages:

```
BiocManager::install(c("BB", "pracma", "lme4", "outliers", "BiocParallel"))
```

Finally, load the package:

```
library(protlocassign)
```

This will make the programs and datasets available. 
Instructions for use of the package are in the vignette.

## References

Jadot, M.; Boonen, M.; Thirion, J.; Wang, N.; Xing, J.; Zhao, C.; 
Tannous, A.; Qian, M.; Zheng, H.; Everett, J. K., 
Accounting for protein subcellular localization: A compartmental map of 
the rat liver proteome. 
Molecular & Cellular Proteomics 2017, 16, (2), 194-212.

Moore DF, Sleat D, Lobel P. A method to estimate the distribution of 
proteins across multiple compartments using data from quantitative 
proteomics subcellular fractionation experiments. 
Journal of Proteome Research, 2022, 21, 6, 1371-1381 
https://doi.org/10.1021/acs.jproteome.1c00781

Tannous, A.; Boonen, M.; Zheng, H.; Zhao, C.; Germain, C. J.; 
Moore, D. F.; Sleat, D. E.; Jadot, M.; Lobel, P., 
Comparative analysis of quantitative mass spectrometric methods 
for subcellular proteomics. 
J Proteome Res 2020, 19, (4), 1718-1730.


