library(protlocassign)

# fitCPA, a key component of the package, computes
#  assignment proportions for proteins based on
#  a set of pre-assigned marker proteins
test_that("CPA routine makes correct assignments for test protein", {
data(protNSA_test)
data(markerListJadot)
refLocationProfilesNSA <- locationProfileSetup(profile=protNSA_test,
                                               markerList=markerListJadot, 
                                               numDataCols=9)
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

})


# used to assess the goodness-of-fit of a weighted mixture of 
# reference profiles to a specified profile but with some 
# proportions constrained to zero. 
# There is a subset of varying parameters while other parameters 
# are fixed at a particular value (typically zero).

test_that("CPA routine makes correct assignments for test protein
          with only specified components allowed to vary", {
data(protRSA_test)
data(refLocProfRSA)


protCPAfromRSA_i_out1b <- fCPAone(profile=protRSA_test[1,],
                                  refLocationProfiles=refLocProfRSA,
                                  numDataCols=9, startProps=NULL,
                                  maxit=10000,
                                  ind.vary=c(2,4), minVal=FALSE)
round(protCPAfromRSA_i_out1b, digits=4)

test_cpaC <- round(protCPAfromRSA_i_out1b, digits=4)
testOutC <- c(0.0000,  0.1415,  0.0000,  0.8585,  0.0000,  0.0000,  0.0000,  
              0.0000, 41.0000,  8.0000)
for (j in 1:10) {
  expect_equal(test_cpaC[j], testOutC[j])
}
})



test_that("CPA routine makes correct assignments for test protein
          with all combinations of a specified 
          number of components allowed to vary", {

            data(protRSA_test)
            data(refLocProfRSA)
            # all combinations of pairs of compartments
            protCPAfromRSA_out2 <- fCPAsubsets(profile=protRSA_test[1,],
                            refLocationProfiles=refLocProfRSA,
                            numDataCols=9, startProps=NULL, nCPAcomparts=2)
            testOut2 <- round(head(protCPAfromRSA_out2), digits=4)
            
            # correct values for first row:
            ERLyso <- c(0.0000, 0.1415, 0.0000, 0.8585, 0.000, 0.0000,    
                        0, 0.0000, 0.0920)
            for (j in 1:9) {
              expect_equal(testOut2[1,j], ERLyso[j])
            }
          })

test_that("projSimplex projects vector to simplex with 
          constraint on sum", {
            y2 <- c(0.1, 0.3)
            p2 <- projSimplex(y2)
            expect_equal(p2, c(0.4, 0.6))
})
            