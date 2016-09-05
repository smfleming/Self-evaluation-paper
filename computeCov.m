function bigSigma = computeCov(sigma1, sigma2, rho)
% function bigSigma = computeCov(sigma1, sigma2, rho)
% Compute covariance matrix

bigSigma =  [sigma1^2 rho*sigma1*sigma2;...
             rho*sigma1.*sigma2 sigma2^2];

