# Semiparametric-Estimation-of-Fractal-Index
MATLAB code accompanying the paper Bennedsen (2020): "Semiparametric estimation and inference on the fractal index of Gaussian and conditionally Gaussian time series data”, 2020. Econometric Reviews, Volume 39, Issue 9, p. 875-903.

The code includes 3 files (additional comments are commented into the files):

(1) illustration_run.m:
This script illustrates the use of the two functions below, i.e., the estimators of the fractal index. As an initial example, try running this script with sig2_u = 0 (no noise) and sig2_u = 0.05 (noise) to see the difference this makes to the estimators and the test for the presence of noise.

(2) estimateFractalIndex.m:
The Matlab command "[a_hat,stddev] = estimateFractalIndex(X,p,m,B)” will calculate an estimate (and associated standard deviation) of the fractal index of the process underlying the data X, as explained in the paper, cf. Theorem 3.1. The (optional) parameters are: p, the power parameter used in the calculation of the variogram; m, the bandwidth; B, the number of simulations to use to estimate the asymptotic variance of the estimate (cf. Appendix B).

(3) estimateFractalIndexRobust.m:
The Matlab command "[a_star,stddev_star,A] = estimateFractalIndexRobust(X,p,m,B,kap)” will calculate an estimate which is robust to noise in the data (and associated standard deviation, which is only valid in the case of no noise) of the fractal index of the process underlying the data X, as explained in the paper, cf. Theorem 3.2. The value “A” is the test statistic for the test of no noise in the data (cf. Theorem 3.3 and Corollary 3.2). Under the null of no noise A is asymptotically standard normally distributed. The (optional) parameters are: p, the power parameter used in the calculation of the variogram; m, the bandwidth; B, the number of simulations to use to estimate the asymptotic variance of the estimate (cf. Appendix B); kap, the gap-parameter used for the robust estimator.

PLEASE NOTE. An additional file is needed for the above to run, which has to be supplied by the user: namely, a function "simFBM(H,n)" that simulates n+1 observations of an fBm with Hurst index H. Many such function can be found on the web (cf., e.g., MathWorks); alternatively, it can be coded up using the methods in the book “Stochastic Simulation” by Asmussen and Glynn (2007, Springer NY), Section XI.3 and XI.6)
