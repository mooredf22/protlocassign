# test outlier search function

test_that("Test outlier search function", {
  set.seed(17356)
  eps <- 0.029885209
  data(spectraNSA_test)
  flagSpectraBox <- outlierFind(protClass=spectraNSA_test,
                                outlierLevel='peptide', numRefCols=5,
                                numDataCols=9,
                                outlierMeth='boxplot', range=3, eps=eps,
                                randomError=TRUE, multiprocess=FALSE)
  
  # examine breakdown of spectral according to the number of fractions 
  #  in their profiles that are outliers
  outlierTable <- as.numeric(table(flagSpectraBox$outlier.num.spectra))
  outlierCorrect <- c(6244, 345, 74, 24, 6, 2, 1)
  expect_equal(outlierTable, outlierCorrect)
  
})