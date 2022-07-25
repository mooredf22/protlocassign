# Qfun4 is an internal function used to assess the goodness-of-fit
# of a weighted mixture of reference profiles to a specified profile
# This tests that certain candidate CPA values 
# return the pre-computed goodness-of-fit values

test_that("Qfun4 returns correct pre-computed goodness-of-fit values", {
A <-c(1, 0, 0, 0)
B <- c(0, .5, 0.5,0)
C <-c(0, 0, 0, 1)
gmat<-cbind(A, B,C)
# make vector for a specified profile to be evaluated
yy <- c(0.2, 0.15, 0.15, 0.5)
# make vector for candidate CPA value to be evaluated
pvec1 <- c(1, 0, 0)  # 100% in compartment A, 0% in compartment B2,
#  0% in compartment C
pvec2<-c(0, 1, 0)  # 0% in compartment A, 100% in compartment B2, 
#  0% in compartment C
pvec3<-c(0,0,1)  # 0% in compartment A, 0% in compartment B2, 
#  100% in compartment C
pvec4<-c(0.2, 0.3, 0.5) #20% in compartment A, 30% in compartment B2, 
#  50% in compartment C

expect_equal(Qfun4(pvec1, yy, gmat), 0.935)
expect_equal(Qfun4(pvec2, yy, gmat), 0.535)
expect_equal(Qfun4(pvec3, yy, gmat), 0.335)
expect_equal(Qfun4(pvec4, yy, gmat), 0)

})


test_that("Qfun4subset returns correct pre-computed goodness-of-fit 
          values with some parameters fixed", {
# Suppose a profile consists of four fractions and there are three
#  reference compartment, designated A, B, and C

A <- c(1, 0, 0, 0)
B <- c(0, .5, 0.5,0)
C <- c(0, 0, 0, 1)
gmat<-cbind(A, B,C)
yy <- c(0.2, 0.15, 0.15, 0.5)
# Make vector for a candidate CPA value to be evaluated
pvec.vary <- c(0.3, 0.7)  # 30%  in compartment B, 70% in compartment C
ind.vary <- c(2,3)  # compartments 2 and 3 may vary
ind.fixed <- 1   # fix value of compartment 1
par.fixed <- 0   # fix that value at 0

expect_equal(Qfun4subset(pvec.vary, yy, gmat, ind.vary=ind.vary,
            ind.fixed=ind.fixed, par.fixed=par.fixed),
            0.08)
})