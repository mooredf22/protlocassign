

test_that("test 50-50 mixture cpa fit and area error", {
  data(refLocProfAcup)
  data(refLocProfNSA)
  data(totProtAT5)

  # Compute relative amount of each theoretical protein that resides 
  #    in each fraction in a mixture set
  mixCytoLysoAcup <- proteinMix(AcupRef=refLocProfAcup,
                              increment.prop=0.1,
                              Loc1=1, Loc2=4)

  # Convert theoretical protein profiles to NSA
  mixCytoLysoNSA <- NSAfromAcup(Acup=mixCytoLysoAcup,
                              NstartMaterialFractions=6, totProt=totProtAT5)

  # Find constrained proportional values
  mixCytoLysoCPAfromNSA <- fitCPA(profile=mixCytoLysoNSA,
                                refLocationProfiles=refLocProfNSA,
                                numDataCols=9)
  # select the row with a 50/50 mixture
  mix50_50actual <- round(as.numeric(mixCytoLysoCPAfromNSA[6,]), digits=4)
  mix50_50expect <- c(0.2122, 0, 0, 0.7878, 0, 0, 0, 0)
  expect_equal(mix50_50actual, mix50_50expect)
  
  # calculate the mixture error
  errorActual <- mixtureAreaError(mixProtiProtjCPA=mixCytoLysoCPAfromNSA,
                 NstartMaterialFractions=6, Loc1=1, Loc2=4,
                 increment.prop=0.1)
  expect_equal(round(errorActual, digits=4), 0.4078)
  
})


