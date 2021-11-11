# This code gives an example of how the maringal
# distribution could be written as a user-defined distribution.

# We want the complementary (right-tail) CDF of the Weibull.

dccdf_weib <- nimbleFunction(
  run = function(x = double(), shape = double(), scale = double(),
                 log = integer(0, default = 0)) {
     # One could use pweibull or write out the equation.
    ans <- pweibull(x, shape, scale, lower.tail = FALSE, log.p = log)
    return(ans)
    returnType(double())
  }
)

C_dccdf_weib <- compileNimble(dccdf_weib)

testModel <- nimbleModel(
  nimbleCode({y ~ dccdf_weib(shape = 2, scale = 3)}),
  data = list(y = 1.5))

testModel$calculate()
pweibull(1.5, shape = 2, scale = 3, lower.tail = FALSE, log.p = TRUE)
# That shows that the log probability calculation for y is correct.
CtestModel <- compileNimble(testModel)
# That shows that the model compiles.
CtestModel$calculate()
