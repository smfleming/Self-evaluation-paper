function conf = sampleMetaConf_singleTrial(xp, action, sigma, sigma_a, sigma_p)
% function conf = selfSelf_singleTrial(xp, action, sigma, sigma_a, sigma_p)
% 
% Runs single-trial inference (via JAGS) of meta-inference model
% xp - perceptual sample
% aciton - subject's action
%
% sigma = shared noise (SD)
% sigma_a = action noise (SD)
% sigma_p = perception noise (SD)
% 
% SF 2014

if nargin < 1
    action = 1;
    xp = 0.5;
    sigma = 0.7;
    sigma_a = 1;
    sigma_p = 1;
elseif nargin > 0 & nargin < 5
    error('Data and precisions not specified')
end

%% Sampling
% MCMC Parameters
mcmc_params.nchains = 1; % How Many Chains?
mcmc_params.nburnin = 10; % How Many Burn-in Samples?
mcmc_params.nsamples = 1000;  %How Many Recorded Samples?
mcmc_params.nthin = 1; % How Often is a Sample Recorded?
mcmc_params.doparallel = 0; % Parallel Option
mcmc_params.dic = 1;
% Initialize Unobserved Variables
for i=1:mcmc_params.nchains
    S.x = 0;
    mcmc_params.init0(i) = S;
end

% Assign variables to the observed nodes
datastruct = struct('x2', xp, 'sigma', sigma, 'sigma1', sigma_a, 'sigma2', sigma_p, 'a', action);

% Select model file and parameters to monitor
model_file = 'sampleMetaConf.txt';
monitorparams = {'d','a'};

% Use JAGS to Sample
[samples, stats] = matjags( ...
    datastruct, ...
    fullfile(pwd, model_file), ...
    mcmc_params.init0, ...
    'doparallel' , mcmc_params.doparallel, ...
    'nchains', mcmc_params.nchains,...
    'nburnin', mcmc_params.nburnin,...
    'nsamples', mcmc_params.nsamples, ...
    'thin', mcmc_params.nthin, ...
    'dic', mcmc_params.dic,...
    'monitorparams', monitorparams, ...
    'savejagsoutput' , 0 , ...
    'verbosity' , 0 , ...
    'cleanup' , 1 , ...
    'workingdir' , 'tmpjags' );
act = samples.a(:);
act(act == 0) = -1;

conf = sum(samples.d(:) == act(:))./length(samples.d(:));


