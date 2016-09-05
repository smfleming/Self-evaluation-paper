% Wrapper to generate confidence distributions
clear all
close all

print_fig = 0;
addpath('~/Dropbox/Utils/graphics/export_fig/')
figDir = '~/Dropbox/Research/Metacognition/BN_model/selfSelf/model fits/';

Ntrials = 1000;
thresholds = linspace(0, 1, 5);
nratings = length(thresholds)-1;
model = 'sigma+rho';

%% Simulate changing sigmaConf
gen_sigmaConf = linspace(0.5, 3, 5);
h1 = figure(1);
colors(:,3) = ones(length(gen_sigmaConf),1);
colors(:,[1 2]) = repmat(linspace(0, 0.8, length(gen_sigmaConf))', 1, 2);
Nsamp = 100000; % how many samples to use for calculating confidence distributions
for i = 1:length(gen_sigmaConf)
    
    % simulate data
    out = genConf(Ntrials, Nsamp, thresholds, 2, gen_sigmaConf(i), 0.7);
    
    nI = sum(out.allCounts(1:nratings));
    nC = sum(out.allCounts(nratings+1:end));
    nTot = nI + nC;
    d1 = norminv(nC./nTot) - norminv(nI./nTot);
    sigma(i) = 2/d1;
    
    pArray = [log(2) log(0.7./(1-0.7))]; % initialise parameters
    [params1(i,:) fit] = fitConf_fmin(out.allCounts+1, Nsamp, sigma(i), model, pArray);
    
    figure(1);
    plot(1:nratings, fit.simProbs(1:nratings), 'Color', colors(i,:), 'LineStyle', '--', 'LineWidth', 2);
    hold on
    plot(nratings+1:nratings*2, fit.simProbs(nratings+1:nratings*2), 'Color', colors(i,:), 'LineWidth', 2);
    plot(out.allCounts./sum(out.allCounts), 'o ', 'Color', colors(i,:), 'MarkerSize', 8);
    hold on
    
end
figure(1);
set(gca, 'XTickLabel', 1:nratings, 'XTick', 1:(nratings*2), 'FontSize', 14)
ylabel('P(conf = y)')
xlabel('Confidence')
title('Sigma_{conf}')
box off
if print_fig
    export_fig([figDir 'paramRecovery_sigmaConf.png'], '-transparent', '-painters', h1)
end

%% Simulate changing rho
gen_rho = [0.2 0.5 0.8];
h2 = figure;
colors = ones(length(gen_sigmaConf),1);
colors(:,[2 3]) = repmat(linspace(0, 0.8, length(gen_sigmaConf))', 1, 2);
for i = 1:length(gen_rho)
    
    % simulate data
    out = genConf(Ntrials, Nsamp, thresholds, 2, 1.5, gen_rho(i));
    
    nI = sum(out.allCounts(1:nratings));
    nC = sum(out.allCounts(nratings+1:end));
    nTot = nI + nC;
    d1 = norminv(nC./nTot) - norminv(nI./nTot);
    sigma = 2/d1;
    
    pArray = [log(2) log(0.7./(1-0.7))]; % initialise parameters
    [params2(i,:) fit] = fitConf_fmin(out.allCounts+1, Nsamp, sigma, model, pArray);
    
    plot(1:nratings, fit.simProbs(1:nratings), 'Color', colors(i,:), 'LineStyle', '--', 'LineWidth', 2);
    hold on
    plot(nratings+1:nratings*2, fit.simProbs(nratings+1:nratings*2), 'Color', colors(i,:), 'LineWidth', 2);
    plot(out.allCounts./sum(out.allCounts), 'o ', 'Color', colors(i,:), 'MarkerSize', 8);
    hold on
end
set(gca, 'XTickLabel', 1:nratings, 'XTick', 1:(nratings*2), 'FontSize', 14)
ylabel('P(conf = y)')
xlabel('Confidence')
title('\rho')
box off
if print_fig
    export_fig([figDir 'paramRecovery_rho.png'], '-transparent', '-painters', h2)
end