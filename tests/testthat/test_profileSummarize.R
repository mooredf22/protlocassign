

test_that("Test profile summarize function on ACP2::CPLQDFLR", {
  set.seed(17356)
  eps <- 0.029885209
  data(spectraNSA_test)
  flagSpectraBox <- outlierFind(protClass=spectraNSA_test,
                       outlierLevel='peptide', numRefCols=5, numDataCols=9,
                       outlierMeth='boxplot', range=3, eps=eps,
                       randomError=TRUE, multiprocess=FALSE)
  pepProfiles <- profileSummarize(protsCombineCnew=flagSpectraBox,
                       numRefCols=6, numDataCols=9, 
                       refColsKeep=c(1,2,4),eps=eps,
                       GroupBy='peptideId', outlierExclude='spectra', 
                       multiprocess = FALSE)
  profActual <- as.numeric(round(pepProfiles[1,3:14], digits=3)) # first row
  profExpected <- c(349, 0.041, 0.053, 0.123, 0.088, 0.03, 0.022, 0.125, 
                  0.451, 0.067, 3, 1)
  expect_equal(profActual, profExpected)
})