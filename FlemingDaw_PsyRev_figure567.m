% Simulated data for effects of action on metacognitive inference
%
% Steve Fleming 2014

clear all
close all

savePlots = 0;
addpath('~/Dropbox/Utils/graphics/export_fig/');
figDir = '~/Dropbox/Research/Metacognition/BN_model/Self-evaluation-paper/figures/';
N = 10000; % trials per coherence level

%% Effect of action on p(dhat=1|Xconf,a) as function of Xconf, second-order
xp_vec = linspace(-2,2,200);
rho = 0.4;
sigma_act = 1;
sigma_conf = 1;

for xpi = 1:length(xp_vec)
    xp = xp_vec(xpi);
    pD_pre(xpi) = computeFirstOrderConf(xp, 1, sigma_act);
    pD_R(xpi) = computeMetaConf(xp, 1, sigma_act, sigma_conf, rho);
    pD_L(xpi) = 1-computeMetaConf(xp, -1, sigma_act, sigma_conf, rho);
end

h1 = figure;
set(gcf, 'Position', [200 200 400 400])
plot(xp_vec, pD_pre, 'Color', [0.5 0.5 0.5], 'LineWidth', 2);
hold on
plot(xp_vec, pD_L, 'k--', 'LineWidth', 2);
plot(xp_vec, pD_R, 'k', 'LineWidth', 2);
set(gca, 'YLim', [0 1], 'FontSize', 20);
xlabel('X_{conf}', 'FontSize', 20);
ylabel('p(d = 1)','FontSize', 20);
legend({'a=?','a=L','a=R'}, 'Location', 'NorthWest')
legend boxoff
box off
if savePlots
    export_fig([figDir 'actionSim_A_psychometric.png'], '-transparent', '-painters', h1)
end

%% Effect of sigma_conf on p(dhat) and confidence
clear pD_pre pD_R pD_L
rho = 0.4;
sigma_act = 1;
sigma_conf_vec = linspace(0.5,2,100);
xp = 0;
for sigmai = 1:length(sigma_conf_vec)
    
    sigma_conf = sigma_conf_vec(sigmai);
    
    % Pre-action confidence is equivalent to max of first-order conf
    mean_cor_pre(sigmai) = max([computeFirstOrderConf(xp, -1, sigma_act), computeFirstOrderConf(xp, 1, sigma_act)]);
    pD_pre(sigmai) = computeFirstOrderConf(xp, 1, sigma_act);
    
    % Post-action confidence
    mean_cor_L(sigmai) = computeMetaConf(xp, -1, sigma_act, sigma_conf, rho);
    mean_cor_R(sigmai) = computeMetaConf(xp, 1, sigma_act, sigma_conf, rho);
    pD_R(sigmai) = mean_cor_R(sigmai);
    pD_L(sigmai) = 1-mean_cor_L(sigmai);
    
end

% Plot post-action confidence as function of sigma_conf
h2 = figure;
set(gcf, 'Position', [200 200 750 350])
subplot(1,2,1)
pos = get(gca, 'Position');
pos(2) = 0.2;
pos(4) = 0.7;
set(gca, 'Position', pos)
plot(sigma_conf_vec, pD_pre, 'Color', [0.5 0.5 0.5], 'LineWidth', 2);
hold on
plot(sigma_conf_vec, pD_L, 'k--', 'LineWidth', 2);
plot(sigma_conf_vec, pD_R, 'k', 'LineWidth', 2);
set(gca, 'YLim', [0 1], 'FontSize', 20);
xlabel('\sigma_{conf}','FontSize', 20);
ylabel('p(d = 1)','FontSize', 20);
legend({'a=?','a=L','a=R'}, 'Location', 'NorthWest')
legend boxoff
box off

subplot(1,2,2)
pos = get(gca, 'Position');
pos(2) = 0.2;
pos(4) = 0.7;
set(gca, 'Position', pos)
plot(sigma_conf_vec, mean_cor_pre, 'Color', [0.5 0.5 0.5], 'LineWidth', 2);
hold on
plot(sigma_conf_vec, mean_cor_L, 'k', 'LineWidth', 2);
set(gca, 'YLim', [0 1], 'FontSize', 20);
xlabel('\sigma_{conf}','FontSize', 20);
ylabel('Confidence','FontSize', 20);
legend({'a=?','a=L or R'}, 'Location', 'NorthWest')
legend boxoff
box off

if savePlots
    export_fig([figDir 'actionSim_B_sigmaConf.png'], '-transparent', '-painters', h2)
end

%% Effect of rho on p(dhat) and confidence
clear pD_pre pD_R pD_L
sigma_act = 1;
sigma_conf = 1;
rho_vec = linspace(0.1,0.9,100);
xp = 0;
for sigmai = 1:length(rho_vec)
    
    rho = rho_vec(sigmai);
    
    % Pre-action confidence is max of first-order conf
    mean_cor_pre(sigmai) = max([computeFirstOrderConf(xp, -1, sigma_act), computeFirstOrderConf(xp, 1, sigma_act)]);
    pD_pre(sigmai) = computeFirstOrderConf(xp, 1, sigma_act);
    
    % Post-action confidence
    mean_cor_L(sigmai) = computeMetaConf(xp, -1, sigma_act, sigma_conf, rho);
    mean_cor_R(sigmai) = computeMetaConf(xp, 1, sigma_act, sigma_conf, rho);
    pD_R(sigmai) = mean_cor_R(sigmai);
    pD_L(sigmai) = 1-mean_cor_L(sigmai);
    
end

% Plot post-action confidence as function of sigma_a
h3 = figure;
set(gcf, 'Position', [200 200 750 350])
subplot(1,2,1)
pos = get(gca, 'Position');
pos(2) = 0.2;
pos(4) = 0.7;
set(gca, 'Position', pos)
plot(rho_vec, pD_pre, 'Color', [0.5 0.5 0.5], 'LineWidth', 2);
hold on
plot(rho_vec, pD_L, 'k--', 'LineWidth', 2);
plot(rho_vec, pD_R, 'k', 'LineWidth', 2);
set(gca, 'YLim', [0 1], 'FontSize', 20);
xlabel('\rho','FontSize', 20);
ylabel('p(d = 1)','FontSize', 20);
legend({'a=?','a=L','a=R'}, 'Location', 'NorthEast')
legend boxoff
box off

subplot(1,2,2)
pos = get(gca, 'Position');
pos(2) = 0.2;
pos(4) = 0.7;
set(gca, 'Position', pos)
plot(rho_vec, mean_cor_pre, 'Color', [0.5 0.5 0.5], 'LineWidth', 2);
hold on
plot(rho_vec, mean_cor_L, 'k', 'LineWidth', 2);
set(gca, 'YLim', [0 1], 'FontSize', 20);
xlabel('\rho','FontSize', 20);
ylabel('Confidence','FontSize', 20);
legend({'a=?','a=L or R'}, 'Location', 'NorthEast')
legend boxoff
box off

if savePlots
    export_fig([figDir 'actionSim_C_rho.png'], '-transparent', '-painters', h3)
end

%% Simulate rate-choice experimental predictions
theta = [0 0.032 0.064 0.128 0.256 0.512 1];
rho = 0.6;
sigma_act = 1;
sigma_conf = 1;
bigSigma = computeCov(sigma_act, sigma_conf, rho);

% Simulate x's, decisions and confidence for first order + second order
meta_cor_rateChoose = nan(length(theta), N);
meta_cor_chooseRate = nan(length(theta), N);
first_cor_rateChoose = nan(length(theta), N);
first_cor_chooseRate = nan(length(theta), N);
d = nan(length(theta), N);
a = nan(length(theta), N);

for s = 1:length(theta)
    
    for i = 1:N
        
        if rand > 0.5
            d(s,i) = 1;
        else
            d(s,i) = -1;
        end
        
        r = mvnrnd([d(s,i).*theta(s) d(s,i).*theta(s)], bigSigma, 1);
        xa = r(1);
        xp = r(2);
        
        if xa > 0
            a(s,i) = 1;
        else
            a(s,i) = -1;
        end
        
        % Compute rate-choose confidence, second-order
        meta_cor_rateChoose(s,i) = max([computeMetaConf(xp, -1, sigma_act, sigma_conf, rho) computeMetaConf(xp, 1, sigma_act, sigma_conf, rho)]);
        
        % Compute choose-rate confidence, second-order
        meta_cor_chooseRate(s,i) = computeMetaConf(xp, a(s,i), sigma_act, sigma_conf, rho);
        
        % Compute rate-choose confidence, first-order
        first_cor_rateChoose(s,i) = max([computeFirstOrderConf(xa, -1, sigma_act) computeFirstOrderConf(xa, 1, sigma_act)]);
        
        % Compute choose-rate confidence, first-order
        first_cor_chooseRate(s,i) = computeFirstOrderConf(xa, a(s,i), sigma_act);
    end
end

for s = 1:length(theta)
    meta_conf_rateChoose_err(s) = mean(meta_cor_rateChoose(s, d(s,:) ~= a(s,:)));
    meta_conf_rateChoose_cor(s) = mean(meta_cor_rateChoose(s, d(s,:) == a(s,:)));
    meta_conf_chooseRate_err(s) = mean(meta_cor_chooseRate(s, d(s,:) ~= a(s,:)));
    meta_conf_chooseRate_cor(s) = mean(meta_cor_chooseRate(s, d(s,:) == a(s,:)));
    first_conf_rateChoose_err(s) = mean(first_cor_rateChoose(s, d(s,:) ~= a(s,:)));
    first_conf_rateChoose_cor(s) = mean(first_cor_rateChoose(s, d(s,:) == a(s,:)));
    first_conf_chooseRate_err(s) = mean(first_cor_chooseRate(s, d(s,:) ~= a(s,:)));
    first_conf_chooseRate_cor(s) = mean(first_cor_chooseRate(s, d(s,:) == a(s,:)));
end

h4a = figure;
set(gcf, 'Position', [200 200 400 400])
plot(theta, meta_conf_chooseRate_cor, 'g', 'LineWidth', 2);
hold on
plot(theta, meta_conf_chooseRate_err, 'r', 'LineWidth', 2);
plot(theta, meta_conf_rateChoose_cor, 'g--', 'LineWidth', 2);
plot(theta, meta_conf_rateChoose_err, 'r--', 'LineWidth', 2);
xlabel('Stimulus strength (\theta)', 'FontSize', 20);
ylabel('Confidence','FontSize', 20);
set(gca, 'YLim', [0.5 1], 'FontSize', 20);
legend({'Choose-rate cor', 'Choose-rate err', 'Rate-choose cor', 'Rate-choose err'}, 'Location', 'SouthWest');
legend boxoff
box off
if savePlots
    export_fig([figDir 'actionSim_D_rateChoose_meta.png'], '-transparent', '-painters', h4a)
end

% Just show one theta level similar to Wierzchon data
h4b = figure;
set(gcf, 'Position', [200 200 400 400])
plot([0.95 1.95], [meta_conf_chooseRate_cor(end) meta_conf_chooseRate_err(end)], 'k^-', 'LineWidth', 2, 'MarkerSize', 12)
hold on
plot([1.05 2.05], [meta_conf_rateChoose_cor(end) meta_conf_rateChoose_err(end)], 'ko--', 'LineWidth', 2, 'MarkerSize', 12)
legend('Choose-rate', 'Rate-choose', 'FontSize', 20)
ylabel('Confidence (a.u.)', 'FontSize', 20)
set(gca, 'YLim', [0.5 1], 'XTick', [1 2], 'XTickLabel', {'Correct', 'Error'}, 'FontSize', 20)
box off
legend boxoff
if savePlots
    export_fig([figDir 'actionSim_D_rateChoose_meta_snapshot.png'], '-transparent', '-painters', h4b)
end

h4c = figure;
set(gcf, 'Position', [200 200 400 400])
plot(theta, first_conf_chooseRate_cor, 'g', 'LineWidth', 2);
hold on
plot(theta, first_conf_chooseRate_err, 'r', 'LineWidth', 2);
plot(theta, first_conf_rateChoose_cor, 'g--', 'LineWidth', 2);
plot(theta, first_conf_rateChoose_err, 'r--', 'LineWidth', 2);
xlabel('Stimulus strength (\theta)', 'FontSize', 20);
ylabel('Confidence','FontSize', 20);
set(gca, 'YLim', [0.5 1], 'FontSize', 20);
legend({'Choose-rate cor', 'Choose-rate err', 'Rate-choose cor', 'Rate-choose err'}, 'Location', 'SouthWest');
legend boxoff
box off
if savePlots
    export_fig([figDir 'actionSim_D_rateChoose_firstOrder.png'], '-transparent', '-painters', h4c)
end

%% Compute dependence of rate-choose and choose-rate bias/resolution on parameter settings for single theta value
theta = 1;
rho = 0.6;
sigma_act = 1;
sigma_conf = 1;

% variable sigma_conf
sigma_conf_vec = linspace(0.5,1.5,10);
for sigmai = 1:length(sigma_conf_vec)
    meta_cor_rateChoose = nan(length(theta), N);
    meta_cor_chooseRate = nan(length(theta), N);
    first_cor_rateChoose = nan(length(theta), N);
    first_cor_chooseRate = nan(length(theta), N);
    d = nan(length(theta), N);
    a = nan(length(theta), N);
    sigma_conf = sigma_conf_vec(sigmai);
    bigSigma = computeCov(sigma_act, sigma_conf, rho);
    
    for i = 1:N
        
        if rand > 0.5
            d(i) = 1;
        else
            d(i) = -1;
        end
        
        r = mvnrnd([d(i) d(i)], bigSigma, 1);
        xa = r(1);
        xp = r(2);
        
        if xa > 0
            a(i) = 1;
        else
            a(i) = -1;
        end
        
        % Compute rate-choose confidence, second-order
        meta_cor_rateChoose(i) = max([computeMetaConf(xp, -1, sigma_act, sigma_conf, rho) computeMetaConf(xp, 1, sigma_act, sigma_conf, rho)]);
        
        % Compute choose-rate confidence, second-order
        meta_cor_chooseRate(i) = computeMetaConf(xp, a(i), sigma_act, sigma_conf, rho);
        
    end
    % Compute difference in bias and resolution for this particular setting
    % of parameters
    bias_rateChoose(sigmai) = mean(meta_cor_rateChoose);
    bias_chooseRate(sigmai) = mean(meta_cor_chooseRate);
    resolution_rateChoose(sigmai) = mean(meta_cor_rateChoose(d == a)) - mean(meta_cor_rateChoose(d ~= a));
    resolution_chooseRate(sigmai) = mean(meta_cor_chooseRate(d == a)) - mean(meta_cor_chooseRate(d ~= a));
    
end
h5 = figure;
set(gcf, 'Position', [200 200 400 400])
plot(sigma_conf_vec, bias_chooseRate, 'k', 'LineWidth', 2);
hold on
plot(sigma_conf_vec, bias_rateChoose, 'k--', 'LineWidth', 2);
xlabel('\sigma_{conf}', 'FontSize', 20);
ylabel('Bias','FontSize', 20);
set(gca, 'YLim', [0.5 1], 'FontSize', 20);
legend({'Choose-rate', 'Rate-choose'}, 'Location', 'SouthWest');
legend boxoff
box off
if savePlots
    export_fig([figDir 'actionSim_sigmaConf_bias.png'], '-transparent', '-painters', h5)
end

h6 = figure;
set(gcf, 'Position', [200 200 400 400])
plot(sigma_conf_vec, resolution_chooseRate, 'k', 'LineWidth', 2);
hold on
plot(sigma_conf_vec, resolution_rateChoose, 'k--', 'LineWidth', 2);
xlabel('\sigma_{conf}', 'FontSize', 20);
ylabel('Correct conf - error conf','FontSize', 20);
set(gca, 'FontSize', 20);
legend({'Choose-rate', 'Rate-choose'}, 'Location', 'NorthEast');
legend boxoff
box off
if savePlots
    export_fig([figDir 'actionSim_sigmaConf_resol.png'], '-transparent', '-painters', h6)
end

% variable rho
rho_vec = linspace(0.1,0.9,10);
theta = 1;
rho = 0.6;
sigma_act = 1;
sigma_conf = 1;
clear bias_rateChoose bias_chooseRate resolution_rateChoose resolution_chooseRate
for sigmai = 1:length(rho_vec)
    meta_cor_rateChoose = nan(length(theta), N);
    meta_cor_chooseRate = nan(length(theta), N);
    first_cor_rateChoose = nan(length(theta), N);
    first_cor_chooseRate = nan(length(theta), N);
    d = nan(length(theta), N);
    a = nan(length(theta), N);
    rho = rho_vec(sigmai);
    bigSigma = computeCov(sigma_act, sigma_conf, rho);
    
    for i = 1:N
        
        if rand > 0.5
            d(i) = 1;
        else
            d(i) = -1;
        end
        
        r = mvnrnd([d(i) d(i)], bigSigma, 1);
        xa = r(1);
        xp = r(2);
        
        if xa > 0
            a(i) = 1;
        else
            a(i) = -1;
        end
        
        % Compute rate-choose confidence, second-order
        meta_cor_rateChoose(i) = max([computeMetaConf(xp, -1, sigma_act, sigma_conf, rho) computeMetaConf(xp, 1, sigma_act, sigma_conf, rho)]);
        
        % Compute choose-rate confidence, second-order
        meta_cor_chooseRate(i) = computeMetaConf(xp, a(i), sigma_act, sigma_conf, rho);
        
    end
    % Compute difference in bias and resolution for this particular setting
    % of parameters
    bias_rateChoose(sigmai) = mean(meta_cor_rateChoose);
    bias_chooseRate(sigmai) = mean(meta_cor_chooseRate);
    resolution_rateChoose(sigmai) = mean(meta_cor_rateChoose(d == a)) - mean(meta_cor_rateChoose(d ~= a));
    resolution_chooseRate(sigmai) = mean(meta_cor_chooseRate(d == a)) - mean(meta_cor_chooseRate(d ~= a));
    
end
h7 = figure;
set(gcf, 'Position', [200 200 400 400])
plot(rho_vec, bias_chooseRate, 'k', 'LineWidth', 2);
hold on
plot(rho_vec, bias_rateChoose, 'k--', 'LineWidth', 2);
xlabel('\rho', 'FontSize', 20);
ylabel('Bias','FontSize', 20);
set(gca, 'YLim', [0.5 1], 'FontSize', 20);
legend({'Choose-rate', 'Rate-choose'}, 'Location', 'SouthWest');
legend boxoff
box off
if savePlots
    export_fig([figDir 'actionSim_rho_bias.png'], '-transparent', '-painters', h7)
end

h8 = figure;
set(gcf, 'Position', [200 200 400 400])
plot(rho_vec, resolution_chooseRate, 'k', 'LineWidth', 2);
hold on
plot(rho_vec, resolution_rateChoose, 'k--', 'LineWidth', 2);
xlabel('\rho', 'FontSize', 20);
ylabel('Correct conf - error conf','FontSize', 20);
set(gca, 'FontSize', 20);
legend({'Choose-rate', 'Rate-choose'}, 'Location', 'NorthEast');
legend boxoff
box off
if savePlots
    export_fig([figDir 'actionSim_rho_resol.png'], '-transparent', '-painters', h8)
end
