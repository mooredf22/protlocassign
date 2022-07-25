library(protlocassign)

test_that("sequence of transformations return original", {
data(protNSA_test)
data(totProtAT5)
# inputting NSA data into "AcupFromNSA" and then into "NSAfromAcup"
#  should yield the original NSA data.
protAcup <- AcupFromNSA(protNSA_test[,seq_len(9)], 
                        NstartMaterialFractions = 6, totProt=totProtAT5)
protNSA <- NSAfromAcup(protAcup, 
                       NstartMaterialFractions = 6, totProt=totProtAT5)
expect_equal(protNSA, protNSA_test[,seq_len(9)])


#data(protNSA_test)
#data(totProtAT5)
# inputting NSA data into "AcupFromNSA", then into "RSAfromAcup",
# and finally into "NSAfromRSA"
#  should yield the original NSA data.
protAcup <- AcupFromNSA(protNSA_test[,seq_len(9)], 
                        NstartMaterialFractions = 6, totProt=totProtAT5)
protRSA <- RSAfromAcup(protAcup, 
                       NstartMaterialFractions = 6, totProt=totProtAT5)
protNSA <- NSAfromRSA(protRSA)
expect_equal(protNSA, protNSA_test[,seq_len(9)] )
})

test_that("reverse sequence of transformations return original", {
# inputting NSA data into "RSAFromNSA", then into "AcupfromRSA",
# and finally into "NSAfromAcup"
#  should yield the original NSA data.
protRSA <- RSAfromNSA(protNSA_test[,seq_len(9)], 
                        NstartMaterialFractions = 6, totProt=totProtAT5)
protAcup <- AcupFromRSA(protRSA, 
                       NstartMaterialFractions = 6, totProt=totProtAT5)
protNSA <- NSAfromAcup(protAcup, 
                       NstartMaterialFractions = 6, totProt=totProtAT5)
expect_equal(protNSA, protNSA_test[,seq_len(9)] )
})




