%%% -----------------------------------------------------------------------
%%% This script illustrates the use of two estimators of the fractal index
%%% of a continuous-time stochastic process.
%%%
%%% See the README file, comments below, and the comments in the 
%%% subfunctions for more information.
%%%
%%% As an illustration: Try running the code with sig2_u = 0 (H0, no noise) 
%%% and sig2_u = 0.05 (H1, with noise). This will illustrate the effect of
%%% noise on the estimators.
%%%
%%% -----------------------------------------------------------------------
%%%
%%% REQUIRED: Function "simFBM(H,n)" that simulates n+1 observations of an 
%%%           fBm with Hurst index H.
%%%
%%% -----------------------------------------------------------------------
%%%
%%% Copyright: Mikkel Bennedsen, June 12, 2019.
%%% 
%%% Please cite: Bennedsen (2019): "Semiparametric estimation and inference
%%% on the fractal index of Gaussian and conditionally Gaussian time series
%%% data"
%%% -----------------------------------------------------------------------

clear; close all; clc;
%% Initialization
rng(42);

alpha = -0.20; % True fractal/roughness index for the simulated data.
sig2_u = 0.00; % Variance of noise. Set equal to 0 for no noise (H0) and >0 for noise in the data (H1).

n = 2000;      % Number of observations.

p = 2;         % Power parameter for variogram
m = 5;         % Bandwidth
kap = 10;      % "Gap"-kappa-parameter for robust estimator

B = 1000;      % Number of fBm-simulations for approximation of AVAR (cf. Appendix B).


%% Simulate data (requires "simFBM"-function)
X = simFBM(alpha+0.5,n-1); % Simulate n observations

if sig2_u > 0
    X = X + sqrt(sig2_u)*randn(n,1); % Add noise.
end

%% Estimate alpha using the methods of the paper
[a_hat,stddev] = estimateFractalIndex(X,p,m,B); % Estimate fractal index of X using usual OLS estimator (Theorem 3.1)
[a_star,stddev_star,A] = estimateFractalIndexRobust(X,p,m,B,kap); % Estimate fractal index of X using robust OLS estimator (Theorem 3.2). "A" is the test-statistic for the test of noise in the data (Theorem 3.3 and Corollary 3.2)

%% Display
disp(' ');

disp(['True alpha          = ',num2str(alpha)]);

disp(' ');
disp(' ---------- Estimation results  for OLS estimator (Theorem 3.1) ---------- ');
disp(' ');
disp(['Estimate (OLS)      = ',num2str(a_hat)]);
disp(['Bias     (OLS)      = ',num2str(a_hat-alpha)]);
disp(['Std. dev (OLS)      = (',num2str(stddev),')']);

disp(' ');
disp(' ---------- Estimation results for robust OLS estimator (Theorem 3.2) ---------- ');
disp(' ');
disp(['Estimate (robust)   = ',num2str(a_star)]);
disp(['Bias     (robust)   = ',num2str(a_star-alpha)]);
disp(['Std. dev (robust)   = (',num2str(stddev_star),')']); % Note: This estimate of the std. dev. is only valid when no noise in present (cf. Theorem 3.2(ii)).

disp(' ');
disp(' ---------- Test for the presence of noise ---------- ');
disp(' ');
disp(['Test stat           = ', num2str(A)]); % N(0,1) under null of no noise
disp(['p-value (two-sided) = ',num2str(2*normcdf(-abs(A)))])

disp(' ');
if 2*normcdf(-abs(A)) < 0.01
     disp('Reject H0 of no noise at 1 percent level.');
elseif 2*normcdf(-abs(A)) < 0.05
     disp('Reject H0 of no noise at 5 percent level.');
elseif 2*normcdf(-abs(A)) < 0.10  
    disp('Reject H0 of no noise at 10 percent level.');
else
    disp('Cannot reject H0 of no noise at 5 percent level.');
end
