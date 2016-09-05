function [nR_S1 nR_S2 data] = metaModelFit_generate_data(ntrials, nratings, sigma_act, sigma_conf, rho, thresh, inv_sigma_act, inv_sigma_conf, inv_rho)
% Generate data from meta inference model
%
% Output coded as for meta-d model
%
% If only generative sigmas and rho are provided, assume that inverse
% parameters are the same as generative
%
% SF 2014

if nargin < 1
    ntrials = 500;
    nratings = 4;
    sigma_act = 1;
    sigma_conf = 1;
    rho = 0.5;
    thresh = [0.5 0.7 0.9];
end
if nargin < 7
    inv_sigma_act = sigma_act;
    inv_sigma_conf = sigma_conf;
    inv_rho = rho;
end

bigSigma = computeCov(sigma_act, sigma_conf, rho);

% Initialise response arrays
nC_rS1 = zeros(1, nratings);
nI_rS1 = zeros(1, nratings);
nC_rS2 = zeros(1, nratings);
nI_rS2 = zeros(1, nratings);

discrete_ind = [-Inf thresh Inf];
for i = 1:ntrials
    
    flipd = round(rand);
    if flipd == 0
        d = -1;
    else
        d = 1;
    end
    
    r = mvnrnd([d d], bigSigma, 1);
    xa = r(1);
    xp = r(2);
    
    if xa > 0
        a = 1;
    else
        a = 0;
    end
    
    % Get confidence
    conf = computeMetaConf(xp, a, inv_sigma_act, inv_sigma_conf, inv_rho);
    
    % Get discrete confidence
    level = conf > discrete_ind;
    rating = sum(level);
    
    if d == -1 & a == 0
        nC_rS1(rating) = nC_rS1(rating) + 1;
        data.correct(i) = 1;
    elseif d == -1 & a == 1
        nI_rS1(rating) = nI_rS1(rating) + 1;
        data.correct(i) = 0;
    elseif d == 1 & a == 0
        nI_rS2(rating) = nI_rS2(rating) + 1;
        data.correct(i) = 0;
    elseif d == 1 & a == 1
        nC_rS2(rating) = nC_rS2(rating) + 1;
        data.correct(i) = 1;
    end
    
    % Store data as vectors
    data.conf(i) = conf;
    data.c(i) = rating;
    data.d(i) = d;
    data.a(i) = a;
    data.theta(i) = 1;
end

% Sort data into nR_S1 and nR_S2
nR_S1 = [nC_rS1(end:-1:1) nI_rS2];
nR_S2 = [nI_rS1(end:-1:1) nC_rS2];

