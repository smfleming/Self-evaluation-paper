% Check correspondence between rate-choose max and sampling

clear all
close all

addpath(genpath('~/Dropbox/Utils/matjags'));

% Set d
d = 1;
% Evaluate comparision for each parameter
sigma = 1;
sigma_a = 1.5;
sigma_p = 1;
xp_space = linspace(-2,2,20);

for xpi = 1:length(xp_space)
    
    % JAGS
    jags_conf(xpi) = sampleMetaConf_singleTrial_noAction(xp_space(xpi), sigma, sigma_a, sigma_p);
    
    % Analytic
    ana_conf(xpi) = max([computeMetaConf_old(xp_space(xpi), -1, sigma, sigma_a, sigma_p), computeMetaConf_old(xp_space(xpi), 1, sigma, sigma_a, sigma_p)]);
    
end

h = figure;
set(gcf, 'Position', [200 200 800 300]);
subplot(1,3,1)
plot(xp_space, jags_conf', 'LineWidth', 2);
xlabel('xp', 'FontSize', 14);
ylabel('confidence', 'FontSize', 14);
axis square
set(gca, 'FontSize', 12);

subplot(1,3,2)
plot(xp_space, ana_conf', 'LineWidth', 2);
xlabel('xp', 'FontSize', 14);
ylabel('confidence', 'FontSize', 14);
axis square
set(gca, 'FontSize', 12);

subplot(1,3,3)
plot(jags_conf(:), ana_conf(:), 'bx ', 'MarkerSize', 7, 'LineWidth', 1.5);
line([0 1], [0 1], 'LineStyle', '--', 'Color', 'k');
xlabel('Sampling', 'FontSize', 14);
ylabel('Max', 'FontSize', 14);
axis square
set(gca, 'FontSize', 12);

%% Simulate rate-choose experimental predictions under alternative definition without conditioning on action
theta = [0 0.032 0.064 0.128 0.256 0.512 1];
sigma = 1;
sigma_act = 1;
sigma_conf = 1;
N = 1000;

% Simulate x's, decisions and confidence for first order + second order
meta_cor_rateChoose = nan(length(theta), N);
meta_cor_chooseRate = nan(length(theta), N);
d = nan(length(theta), N);
a = nan(length(theta), N);

for s = 1:length(theta)
    
    for i = 1:N
        
        if rand > 0.5
            d(s,i) = 1;
        else
            d(s,i) = -1;
        end
        
        x = normrnd(d(s,i)*theta(s), sigma);
        xa = normrnd(x, sigma_act);
        xp = normrnd(x, sigma_conf);
        
        if xa > 0
            a(s,i) = 1;
        else
            a(s,i) = -1;
        end
        
        % Compute rate-choose confidence, second-order
        meta_cor_rateChoose(s,i) = sampleMetaConf_singleTrial_noAction(xp, sigma, sigma_act, sigma_conf);
        
        % Compute choose-rate confidence, second-order
        meta_cor_chooseRate(s,i) = computeMetaConf_old(xp, a(s,i), sigma, sigma_act, sigma_conf);
        
    end
end

for s = 1:length(theta)
    meta_conf_rateChoose_err(s) = mean(meta_cor_rateChoose(s, d(s,:) ~= a(s,:)));
    meta_conf_rateChoose_cor(s) = mean(meta_cor_rateChoose(s, d(s,:) == a(s,:)));
    meta_conf_chooseRate_err(s) = mean(meta_cor_chooseRate(s, d(s,:) ~= a(s,:)));
    meta_conf_chooseRate_cor(s) = mean(meta_cor_chooseRate(s, d(s,:) == a(s,:)));
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