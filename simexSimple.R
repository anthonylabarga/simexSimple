## This function performs simulation-extrapolation to estimate regression coefficients
## under conditions of covariate noise.
##
## We suppose that we observe a response vector Y and a corrupted covariate vector X_w = X_u + nu
## where nu is a noise vector with mean 0 and variance sigma^2_u. We suppose that Y and X_u
## are related through the simple linear model Y = Beta_0 + Beta_1*X_u+epsilon.
##
## The procedure is to simulate B datasets corrupted by a measurement error (1+theta)sigma^2_u,
## where theta are equidistant values in some set thetaRange. For each theta, B regressions are
## fit with simulated data, and the average coefficient estimates stored as BetaHat. For each
## coefficient, a quadratic curve is fit between theta and BetaHat_i. A final coefficient
## estimate is produced by extrapolating this curve to theta = -1.

require(compiler)

simexSimple <- function( response, covariate, noiseVar, numDatasets = 100, thetaSpace = seq( .1, 2, .1))
{    
  NoiseVec <- sqrt(thetaSpace*noiseVar)
  
  ExtrapolationMatrix <- matrix(0,length(thetaSpace),(NCOL(covariate)+1))
  
  for( i in 1:length(thetaSpace))
  {
    
    for( j in 1:numDatasets)
    {
      simCov <- covariate + rnorm(NROW(covariate),0,NoiseVec[i])
      ExtrapolationMatrix[i,] <- ((j-1)/j)*ExtrapolationMatrix[i,]+(1/j)*coef( lm( response ~ simCov))
    }
  }
  ExtrapolationMatrix <- cbind(thetaSpace, ExtrapolationMatrix)
  betaHat <- c(sum(coef(lm(ExtrapolationMatrix[,2] ~ ExtrapolationMatrix[,1] + I(ExtrapolationMatrix[,1]^2)))*c(1, -1, 1)),
               sum(coef(lm(ExtrapolationMatrix[,3] ~ ExtrapolationMatrix[,1] + I(ExtrapolationMatrix[,1]^2)))*c(1, -1, 1)))
  return(betaHat)
}

simexSimpleC <- cmpfun(simexSimple)
