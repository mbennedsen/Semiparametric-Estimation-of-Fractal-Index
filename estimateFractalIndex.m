function [a_hat,stddev] =  estimateFractalIndex(X,p,m,B)
%%% Function that estimates the fractal index of a process + the standard
%%% deviation of the estimate.
%%%
%%% Copyright: Mikkel Bennedsen, June 12, 2019.
%%% 
%%% Please cite: Bennedsen (2019): "Semiparametric estimation and inference
%%% on the fractal index of Gaussian and conditionally Gaussian time series
%%% data"
%%%
%%% INPUT:
%%% X: data
%%% p: power parameter (default = 2)
%%% m: bandwidth (default = 5)
%%% B: Number of simulations used for estimating AVAR (default = 1000).
%%%
%%% OUTPUT:
%%% a_hat : Estimate of fractal index using Proposition 2.1
%%% stddev: Estimate of the standard deviation of a_hat using Theorem 3.1
%%%
%%% REQUIRED:
%%% Function "simFBM(H,n)" that simulates n+1 observations of an fBm with
%%% Hurst index H.

%% Set defaults
if nargin < 1
    error('Function needs data input');
elseif nargin < 2
    p = 2; 
    m = 5;
    B = 1000;
elseif nargin < 3
    m = 5;
    B = 1000;
elseif nargin < 4
    B = 1000;
end

n = length(X);

%% Some calculations
x_m = (log(1:m) - mean(log(1:m)))';

mp = 2^(p/2)/sqrt(pi)*gamma((p+1)/2);
m2p = 2^(2*p/2)/sqrt(pi)*gamma((2*p+1)/2);

%% Estimate alpha (Proposition 3.1)
gam_hat   = nan(m,1);
for i = 1:m
    gam_hat(i)  = mean( abs(X(1+i:end) - X(1:end-i)).^p ); % Variogram of X.
end

a_hat = x_m'*log(gam_hat)/(x_m'*x_m)/p - 0.5; % Estimate of fractal index.

%% alpha only allowed to be in (-0.5,0.5). Check this:
if a_hat < -0.4995
    a_hat = -0.4995;
    warning('Estimate of alpha < -0.4995. Set to -0.4995.');
elseif a_hat > 0.4995
    a_hat = 0.4995;
    warning('Estimate of alpha > 0.4995. Set to 0.4995.');    
end

%% Calculate std. dev. (Theorem 3.1) if required
if nargout > 1
    
    H_star = a_hat + 0.5; % Value of Hurst index used in approximating AVAR (if testing H0: alpha = alpha0 for some alpha0, then replace "a_hat" with alpha0, cf. Remark 3.2)
    gam2_hat = mean( abs(X(1+1:end) - X(1:end-1)).^(2*p) );
    Sp_hat = sqrt(gam2_hat/m2p)/(gam_hat(1)/mp); % Correction for possible heteroskedasticity / non-Gaussianity through volatility modulation (cf. Proposition 3.2).

    %%% Approximate LAMBDA matrix by Monte Carlo simulation (Appendix B)
    gam_hat_star = nan(B,m);
    for b = 1:B
        Xstar = simFBM(H_star,n-1);

        for i = 1:m
            gam_hat_star(b,i)  = mean( abs(Xstar(1+i:end) - Xstar(1:end-i)).^p );
        end
    end
    LAM = nan(m,m);
    tmp_gam = cov(gam_hat_star);
    for i = 1:m
        for j = 1:m
            LAM(i,j)     = n*tmp_gam(i,j)/((i/n)^(p*H_star)*mp)/((j/n)^(p*H_star)*mp);
        end
    end

    sig2_hat = x_m'*LAM*x_m/(x_m'*x_m)^2/p^2; % Estimate of AVAR from Theorem 3.1

    stddev = sqrt(sig2_hat)*Sp_hat/sqrt(n); % Std. dev. of estimate of alpha (cf. Corollary 3.1)
end