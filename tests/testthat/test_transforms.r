library(protlocassign)

data(protNSA_test)
data(totProtAT5)
protAcup <- AcupFromNSA(protNSA_test[,seq_len(9)], 
                        NstartMaterialFractions = 6, totProt=totProtAT5)
protNSA <- NSAfromAcup(protAcup, 
                       NstartMaterialFractions = 6, totProt=totProtAT5)
expect_equal(protNSA, protNSA_test[,seq_len(9)])


data(protNSA_test)
data(totProtAT5)
protAcup <- AcupFromNSA(protNSA_test[,seq_len(9)], 
                        NstartMaterialFractions = 6, totProt=totProtAT5)
protRSA <- RSAfromAcup(protAcup, 
                       NstartMaterialFractions = 6, totProt=totProtAT5)
protNSA <- NSAfromRSA(protRSA)
expect_equal(protNSA, protNSA_test[,seq_len(9)] )



data(protNSA_test)
data(markerListJadot)
refLocationProfilesNSA <- locationProfileSetup(profile=protNSA_test,
                              markerList=markerListJadot, numDataCols=9)
protCPAfromNSA_test <- fitCPA(profile=protNSA_test["TLN1",],
                              refLocationProfiles=refLocationProfilesNSA,
                              numDataCols=9)
test_cpa <- as.numeric(protCPAfromNSA_test)
expect_equal(test_cpa[1], 0.35747217, tolerance=0.001)
expect_equal(test_cpa[2], 0)
expect_equal(test_cpa[3], 0)
expect_equal(test_cpa[4], 0)
expect_equal(test_cpa[5], 0)
expect_equal(test_cpa[6], 0.06179167, tolerance=0.001)
expect_equal(test_cpa[7], 0)
expect_equal(test_cpa[8], 0.58073616, tolerance=0.001)



