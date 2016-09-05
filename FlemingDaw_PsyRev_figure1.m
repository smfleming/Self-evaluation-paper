% Script to produce demo plots of first-order, post-decisional and
% second-order model confidence and graphical intuitions
%
% For second-order model simulation, the script shows the effect of Xconf sample on inference on Xact
% calls JAGS model to sample conditional distribution when d is unknown (note that this uses alternative parameterisation of covariance matrix in terms of
% sigma, sigma_act and sigma_conf)
% This requires JAGS and matjags to be installed, see http://psiexp.ss.uci.edu/research/programs_data/jags/
%
% Is savePlots is switched on, will write out figure to file
%
% Steve Fleming 2016

clear all
close all

addpath('~/Dropbox/Utils/matjags/')
figDir = '~/Dropbox/Research/Metacognition/BN_model/selfSelf/figures/';

savePlots = 0;

action = 1;
xp = -0.6;
sigma = 1;
sigma_a = 1;
sigma_p = 1;


%% FIGURE 1A - first-order model example
d1 = 1;
d0 = -1;
x = linspace(-3,3,400);
x_samp = 0.4;

d1_dist = normpdf(x, d1, sigma);
d0_dist = normpdf(x, d0, sigma);

% Likelihood
h1 = figure;
set(gcf, 'Position', [400 400 400 300]);
plot(x, d0_dist, 'k--', 'LineWidth', 2);
hold on
plot(x, d1_dist, 'k', 'LineWidth', 2); 
p=patch([0 3 3 0],[0 0 1 1],[0.8 0.8 0.8], 'LineStyle', 'None');
set(p,'FaceAlpha',0.5);
set(gca, 'FontSize', 18, 'YLim', [0 1], 'XLim', [-3 3], 'XTick', [-1 0 1])
xlabel('Decision variable');
ylabel('P(X_{act}|d)');
axis square

ind = dsearchn(x', x_samp);
line([x_samp x_samp], [0 0.7], 'LineStyle', '--', 'Color', [0.5 0.5 0.5], 'LineWidth', 2);
conf = computeFirstOrderConf(x_samp, sign(x_samp), sigma);
text(0.2, 0.9, 'a = R', 'FontSize', 16)
text(-2, 0.8, ['Confidence = ' num2str(round(conf,2))], 'FontSize', 16)
box off
if savePlots
    saveas(h1, [figDir 'firstOrderModel_plot.png'], 'png');
end

%% FIGURE 1B - post-decisional model example
d1 = 1;
d0 = -1;
sigma = 1;
x = linspace(-3,3,400);
x_act = 0.4;
x_conf = -0.6;

d1_dist = normpdf(x, d1, sigma);
d0_dist = normpdf(x, d0, sigma);

% Likelihood
h2 = figure;
set(gcf, 'Position', [400 400 400 300]);
plot(x, d0_dist, 'k--', 'LineWidth', 2);
hold on
plot(x, d1_dist, 'k', 'LineWidth', 2);
p=patch([0 3 3 0],[0 0 1 1],[0.8 0.8 0.8], 'LineStyle', 'None');
set(p,'FaceAlpha',0.5);
set(gca, 'FontSize', 18, 'YLim', [0 1], 'XLim', [-3 3], 'XTick', [-1 0 1])
xlabel('Decision variable');
ylabel('P(X_{conf}|d)');
axis square

ind = dsearchn(x', x_samp);
line([x_act x_act], [0 0.7], 'LineStyle', '--', 'Color', [0.5 0.5 0.5], 'LineWidth', 2);
line([x_conf x_conf], [0 0.7], 'LineStyle', '-', 'Color', [0.5 0.5 0.5], 'LineWidth', 2);
conf = computeFirstOrderConf(x_conf, sign(x_samp), sigma);
text(0.2, 0.9, 'a = R', 'FontSize', 16)
text(-2, 0.8, ['Confidence = ' num2str(round(conf,2))], 'FontSize', 16)
box off
if savePlots
    saveas(h2, [figDir 'postDecisionModel_plot.png'], 'png');
end

%% FIGURE 1C - second-order model example
% MCMC Parameters
mcmc_params.nchains = 1; % How Many Chains?
mcmc_params.nburnin = 1000; % How Many Burn-in Samples?
mcmc_params.nsamples = 100000;  %How Many Recorded Samples?
mcmc_params.nthin = 1; % How Often is a Sample Recorded?
mcmc_params.doparallel = 0; % Parallel Option
mcmc_params.dic = 0;
% Initialize Unobserved Variables
for i=1:mcmc_params.nchains
    S.x = 0;
    mcmc_params.init0(i) = S;
end

%% Plot bivariate PDF
% Assign variables to the observed nodes
datastruct = struct('sigma', sigma, 'sigma1', sigma_a, 'sigma2', sigma_p);

% Select model file and parameters to monitor
model_file = 'sampleMetaConf.txt';
monitorparams = {'x1','x2'};

% Use JAGS to Sample while clamping x2
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
    'verbosity' , 1 , ...
    'cleanup' , 1 , ...
    'workingdir' , 'tmpjags' );

% Approximate bivariate pdf of samples
x1 = samples.x1;
x2 = samples.x2;
d = [mean(x1) mean(x2)];
bigSigma = cov(x1,x2);

xact = linspace(-6, 6, 500);
xconf = xact;
[X1, X2] = meshgrid(xact, xconf);
F = mvnpdf([X1(:) X2(:)],d, bigSigma);
F = reshape(F,length(xconf),length(xact));

h3 = figure;
set(gcf, 'Position', [200 200 400 350], 'colormap', gray)
a1 = axes;
contour(xact, xconf, F)
p=patch([0 6 6 0],[-6 -6 6 6],[0.8 0.8 0.8], 'LineStyle', 'None');
set(p,'FaceAlpha',0.5);
set(gca, 'FontSize', 20, 'XLim', [-6 6]);
line([-6 2], [xp xp], 'LineStyle', '--', 'Color', [0.5 0.5 0.5], 'LineWidth', 2);
xlabel('Decision variable');
ylabel('Confidence variable');
text(3.5, -3, 'a = R', 'FontSize', 16)
box off

%% Use JAGS to sample while clamping x2
% Assign variables to the observed nodes
datastruct = struct('x2', xp, 'sigma', sigma, 'sigma1', sigma_a, 'sigma2', sigma_p, 'a', action);

% Select model file and parameters to monitor
model_file = 'sampleMetaConf.txt';
monitorparams = {'d','a','x1'};

% Use JAGS to Sample while clamping x2
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
    'verbosity' , 1 , ...
    'cleanup' , 1 , ...
    'workingdir' , 'tmpjags' );
act = samples.a(:);
act(act == 0) = -1;
conf = mean(act' == samples.d);

axes('Position', [0.25 0.2 0.65 0.1], 'Color', 'None')
f = ksdensity(samples.x1, xact);
f(xact <= 0) = 0;
plot(xact, f, 'k', 'LineWidth', 2)
box off
set(gca, 'XLim', [-3.5 5], 'YLim', [0 1], 'XTickLabel', [], 'FontSize', 11);
ylabel('P(X_{act}|X_{conf}, a)')

if savePlots
    saveas(h3, [figDir 'secondOrderModel_plot.png'], 'png');
end