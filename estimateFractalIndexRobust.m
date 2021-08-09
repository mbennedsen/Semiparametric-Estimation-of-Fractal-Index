function [a_star,stddev,A] =  estimateFractalIndexRobust(X,p,m,B,kap)
%%% Function that estimates the fractal index of a process + the standard
%%% deviation of the estimate using an estimator robust to noise. Also
%%% calculated is a test statistic for H0: "no noise in the observations".
%%%
%%% Copyright: Mikkel Bennedsen, June 12, 2019.
%%% 
%%% Please cite: Bennedsen (2019): "Semiparametric estimation and inference
%%% on the fractal index of Gaussian and conditionally Gaussian time series
%%% data"
%%%
%%% INPUT:
%%% X  : data
%%% p  : power parameter (default = 2)
%%% m  : bandwidth (default = 5)
%%% B  : Number of simulations used for estimating AVAR (default = 1000)
%%% kap: "Gap"-parameter for use in estimator (default = max(floor(n^0.25),2) )
%%%
%%% OUTPUT:
%%% a_star: Robust estimate of fractal index using Proposition 3.4 (NOTE: Consistent with and without noise in the data!)
%%% stddev: Estimate of the standard deviation of a_star using Theorem 3.2. (NOTE: only consistent under the null of no noise!)
%%% A     : Test statistic for the presence of noise in the data using Corollary 3.2 (asymptotically N(0,1) distributed under the null of no noise).
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

if nargin < 5
    kap = max(floor(n^0.25),2);
end

%% Some calculations
x_m = (log(1:m) - mean(log(1:m)))';

mp = 2^(p/2)/sqrt(pi)*gamma((p+1)/2);
m2p = 2^(2*p/2)/sqrt(pi)*gamma((2*p+1)/2);

%% Estimate alpha (Proposition 3.1)
gam_hat   = nan(m,1);
f_hat   = nan(m,1);
for i = 1:m
    gam_hat(i)  = mean( abs(X(1+i:end) - X(1:end-i)).^p ); % Variogram of X
    gam_hat2    = mean( abs(X(1+kap*i:end) - X(1:end-kap*i)).^p );
                    
    f_hat(i) = abs( gam_hat2^(2/p) - gam_hat(i)^(2/p) );
end

a_hat  = x_m'*log(gam_hat)/(x_m'*x_m)/p - 0.5; % Usual estimate of alpha.
a_star = x_m'*log(f_hat)/(x_m'*x_m)/2 - 0.5;   % Robust estimate of alpha.

%% alpha only allowed to be in (-0.5,0.5). Check this:
if a_hat < -0.4995
    a_hat = -0.4995;
    warning('Estimate of alpha < -0.4995. Set to -0.4995.');
elseif a_hat > 0.4995
    a_hat = 0.4995;
    warning('Estimate of alpha > 0.4995. Set to 0.4995.');    
end

if a_star < -0.4995
    a_star = -0.4995;
    warning('Robust estimate of alpha < -0.4995. Set to -0.4995.');
elseif a_star > 0.4995
    a_star = 0.4995;
    warning('Robust estimate of alpha > 0.4995. Set to 0.4995.');    
end


%% Calculate SIGMA_1 matrix (Appendix C.1) if required
if nargout > 1
    
    k_star = floor(m/kap)+1;
    m_vec = [1:m,k_star*kap:kap:m*kap];
    
    SIG1mat = zeros(m,length(m_vec));    
    if k_star == 1
        SIG1mat(1:m,1:m) = -(2/p)/(kap^(2*a_star+1)-1) * eye(m);
        SIG1mat(k_star:m , m+1:end) = (2/p)*kap^(2*a_star+1)/(kap^(2*a_star+1)-1) * eye(m);
    else
        SIG1mat(1:m,1:m) = -(2/p)/(kap^(2*a_star+1)-1)*eye(m);
        for i = 1:(k_star-1)
            SIG1mat(i,i*kap) = (2/p)*kap^(2*a_star+1)/(kap^(2*a_star+1)-1);
        end
    
        SIG1mat(k_star:m , m+1:end) = (2/p)*kap^(2*a_star+1)/(kap^(2*a_star+1)-1) * eye(m-k_star+1);
    end
    %% Calculate std. dev. (Theorem 3.2)
    H_star = a_star + 0.5; % Value of Hurst index used in approximating AVAR (if testing H0: alpha = alpha0 for some alpha0, then replace "a_hat" with alpha0, cf. Remark 3.2)
    gam2_hat = mean( abs(X(1+1:end) - X(1:end-1)).^(2*p) );
    Sp_hat = sqrt(gam2_hat/m2p)/(gam_hat(1)/mp); % Correction for possible heteroskedasticity / non-Gaussianity through volatility modulation (cf. Proposition 3.2).

    %%% Approximate LAMBDA matrix by Monte Carlo simulation (Appendix B)
    gam_hat_star = nan(B,length(m_vec));
    for b = 1:B
        Xstar = simFBM(H_star,n-1);

        for i = 1:length(m_vec)
            gam_hat_star(b,i)  = mean( abs(Xstar(1+m_vec(i):end) - Xstar(1:end-m_vec(i))).^p );
        end
    end
    LAM = nan(length(m_vec),length(m_vec));
    tmp_gam = cov(gam_hat_star);
    for i = 1:length(m_vec)
        for j = 1:length(m_vec)
            LAM(i,j)     = n*tmp_gam(i,j)/((m_vec(i)/n)^(p*H_star)*mp)/((m_vec(j)/n)^(p*H_star)*mp);
        end
    end

    sig2_star = x_m'*(SIG1mat*LAM*SIG1mat')*x_m/(x_m'*x_m)^2/4; % Estimate of AVAR from Theorem 3.2

    stddev = sqrt(sig2_star)*Sp_hat/sqrt(n); % Std. dev. of estimate of alpha
end

%% Calculate SIGMA_2 matrix (Appendix C.2) if required
if nargout > 2
    
    SIG2mat = zeros(m,length(m_vec));
    if k_star == 1
        SIG2mat(1:m,1:m) = -(2/p)*kap^(2*a_star+1)/(kap^(2*a_star+1)-1) * eye(m);
        SIG2mat(k_star:m , m+1:end) = (2/p)*kap^(2*a_star+1)/(kap^(2*a_star+1)-1) * eye(m);
    else
        SIG2mat(1:m,1:m) = -(2/p)*kap^(2*a_star+1)/(kap^(2*a_star+1)-1)*eye(m);
        for i = 1:(k_star-1)
            SIG2mat(i,i*kap) = (2/p)*kap^(2*a_star+1)/(kap^(2*a_star+1)-1);
        end
        SIG2mat(k_star:m , m+1:end) = (2/p)*kap^(2*a_star+1)/(kap^(2*a_star+1)-1) * eye(m-k_star+1);
    end
    
    sig2_star2 = x_m'*(SIG2mat*LAM*SIG2mat')*x_m/(x_m'*x_m)^2/4; % Estimate of AVAR from Theorem 3.3
    
    A = sqrt(n)*(a_star - a_hat)/Sp_hat/sqrt(sig2_star2); % Test statistic for H0: "no noise in observations". Asymptotically N(0,1) under H0. (Theorem 3.3 and Corollary 3.2)
    
end