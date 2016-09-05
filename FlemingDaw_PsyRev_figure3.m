% Simulations of divergent pattern of confidence ratings as a function of stimulus
% strength and decision accuracy for both the first- and second-order models
%
% Steve Fleming 2015

clear all
close all

addpath('~/Dropbox/Utils/graphics/export_fig/');
figDir = '~/Dropbox/Research/Metacognition/BN_model/selfSelf/figures/';
savePlots = 1;

%% Simulate X-pattern from first- and second-order models

theta = [0 0.032 0.064 0.128 0.256 0.512 1];
sigmaAct = 1;
sigmaConf = 1;
rho = 0.5;
bigSigma = computeCov(sigmaAct, sigmaConf, rho);
N = 10000;    % trials per level of coherence

% Simulate x's, decisions and confidence for meta model
mean_cor = nan(length(theta), N);
xa = nan(length(theta), N);
xp = nan(length(theta), N);
d = nan(length(theta), N);
a = nan(length(theta), N);

for s = 1:length(theta)
    for i = 1:N
        
        if rand > 0.5
            d(s,i) = 1;
        else
            d(s,i) = -1;
        end
        
        % Sample bivariate xact and xconf for second-order computation, use
        % xact for first-order models
        r = mvnrnd([d(s,i).*theta(s) d(s,i).*theta(s)], bigSigma, 1);
        xa(s,i) = r(1);
        xp(s,i) = r(2);
        
        if xa(s,i) > 0
            a(s,i) = 1;
            flip_a = 1;
        else
            a(s,i) = -1;
            flip_a = 0;
        end
        
        % Sample new xnew for postdecisional model
        xnew(s,i) = normrnd(d(s,i).*theta(s), sigmaAct);
        
        % Compute confidence
        firstOrder_mean_cor(s,i) = computeFirstOrderConf(xa(s,i), flip_a, sigmaAct);
        postDecisional_mean_cor(s,i) = computeFirstOrderConf(xa(s,i) + xnew(s,i), flip_a, sqrt(2*(sigmaAct^2)));
        secondOrder_mean_cor(s,i) = computeMetaConf(xp(s,i), flip_a, sigmaAct, sigmaConf, rho);
        
    end
    fprintf('Simulation completed for level %d... \n', s);
end

%% PLOTS
%% First-order

for s = 1:length(theta)
    acc(s) = mean(d(s,:) == a(s,:));
    conf_err(s) = mean(firstOrder_mean_cor(s, d(s,:) ~= a(s,:)));
    conf_cor(s) = mean(firstOrder_mean_cor(s, d(s,:) == a(s,:)));
end

% Figure for divergent confidence as function of theta
h1 = figure;
set(gcf, 'Position', [400 400 400 300]);
plot(theta, conf_cor, 'g', 'LineWidth', 3);
hold on
plot(theta, conf_err, 'r', 'LineWidth', 3);
set(gca, 'FontSize', 18, 'XLim', [0 max(theta)+0.1], 'YLim', [0 1]);
xlabel('Stimulus strength', 'FontSize', 20);
ylabel('Confidence', 'FontSize', 20);
box off
if savePlots
    export_fig([figDir 'firstOrder_conf_behav.png'], '-transparent', '-painters', h1)
end

% Figure for confidence as function of |x|
unzip_d = reshape(d,1,length(theta).*N);
unzip_a = reshape(a,1,length(theta).*N);
unzip_acc = unzip_d == unzip_a;
unzip_abs_x = abs(reshape(xa,1,length(theta).*N));
unzip_conf = reshape(firstOrder_mean_cor,1,length(theta).*N);

% bin into quantiles
q_x = quantile(unzip_abs_x, [0.33 0.66]);
q_x = [0 q_x Inf];
for b = 1:3
    mean_cor_conf(b) = mean(unzip_conf(unzip_abs_x > q_x(b) & unzip_abs_x < q_x(b+1) & unzip_acc));
    mean_err_conf(b) = mean(unzip_conf(unzip_abs_x > q_x(b) & unzip_abs_x < q_x(b+1) & ~unzip_acc));
end

h2 = figure;
set(gcf, 'Position', [400 400 400 300]);
plot(mean_cor_conf, 'g', 'LineWidth', 3);
hold on
plot(mean_err_conf, 'r', 'LineWidth', 3);
set(gca, 'FontSize', 18, 'XTick', [1 2 3], 'XTickLabel', {'Low', 'Med', 'High'}, 'XLim', [0.5 3.5], 'YLim', [0.5 1]);
xlabel('X_{act}', 'FontSize', 20);
ylabel('Confidence', 'FontSize', 20);
box off
if savePlots
    export_fig([figDir 'firstOrder_conf_neural.png'], '-transparent', '-painters', h2)
end

%% Post-decisional

for s = 1:length(theta)
    acc(s) = mean(d(s,:) == a(s,:));
    conf_err(s) = mean(postDecisional_mean_cor(s, d(s,:) ~= a(s,:)));
    conf_cor(s) = mean(postDecisional_mean_cor(s, d(s,:) == a(s,:)));
end

% Figure for divergent confidence as function of theta
h1 = figure;
set(gcf, 'Position', [400 400 400 300]);
plot(theta, conf_cor, 'g', 'LineWidth', 3);
hold on
plot(theta, conf_err, 'r', 'LineWidth', 3);
set(gca, 'FontSize', 18, 'XLim', [0 max(theta)+0.1], 'YLim', [0 1]);
xlabel('Stimulus strength', 'FontSize', 20);
ylabel('Confidence', 'FontSize', 20);
box off
if savePlots
    export_fig([figDir 'postDecisional_conf_behav.png'], '-transparent', '-painters', h1)
end

% Figure for confidence as function of |x|
unzip_d = reshape(d,1,length(theta).*N);
unzip_a = reshape(a,1,length(theta).*N);
unzip_acc = unzip_d == unzip_a;
xconf_post = xa + xnew;
unzip_abs_x = abs(reshape(xconf_post,1,length(theta).*N));
unzip_conf = reshape(postDecisional_mean_cor,1,length(theta).*N);

% bin into quantiles
q_x = quantile(unzip_abs_x, [0.33 0.66]);
q_x = [0 q_x Inf];
for b = 1:3
    mean_cor_conf(b) = mean(unzip_conf(unzip_abs_x > q_x(b) & unzip_abs_x < q_x(b+1) & unzip_acc));
    mean_err_conf(b) = mean(unzip_conf(unzip_abs_x > q_x(b) & unzip_abs_x < q_x(b+1) & ~unzip_acc));
end

h2 = figure;
set(gcf, 'Position', [400 400 400 300]);
plot(mean_cor_conf, 'g', 'LineWidth', 3);
hold on
plot(mean_err_conf, 'r', 'LineWidth', 3);
set(gca, 'FontSize', 18, 'XTick', [1 2 3], 'XTickLabel', {'Low', 'Med', 'High'}, 'XLim', [0.5 3.5], 'YLim', [0.5 1]);
xlabel('X_{conf}', 'FontSize', 20);
ylabel('Confidence', 'FontSize', 20);
box off
if savePlots
    export_fig([figDir 'postDecisional_conf_neural.png'], '-transparent', '-painters', h2)
end

%% Second order
for s = 1:length(theta)
    acc(s) = mean(d(s,:) == a(s,:));
    conf_err(s) = mean(secondOrder_mean_cor(s, d(s,:) ~= a(s,:)));
    conf_cor(s) = mean(secondOrder_mean_cor(s, d(s,:) == a(s,:)));
end

% Figure for divergent confidence as function of theta
h3 = figure;
set(gcf, 'Position', [400 400 400 300]);
plot(theta, conf_cor, 'g', 'LineWidth', 3);
hold on
plot(theta, conf_err, 'r', 'LineWidth', 3);
set(gca, 'FontSize', 18, 'YLim', [0 1], 'XLim', [0 max(theta)+0.1]);
xlabel('Stimulus strength', 'FontSize', 20);
ylabel('Confidence', 'FontSize', 20);
box off
if savePlots
    export_fig([figDir 'secondOrder_conf_behav.png'], '-transparent', '-painters', h3)
end

% Figure for confidence as function of |xconf|
unzip_d = reshape(d,1,length(theta).*N);
unzip_a = reshape(a,1,length(theta).*N);
unzip_acc = unzip_d == unzip_a;
unzip_abs_xp = abs(reshape(xp,1,length(theta).*N));
unzip_conf = reshape(secondOrder_mean_cor,1,length(theta).*N);

% bin into quantiles
q_xp = quantile(unzip_abs_xp, [0.33 0.66]);
q_xp = [0 q_xp Inf];
for b = 1:3
    mean_cor_conf(b) = mean(unzip_conf(unzip_abs_xp > q_xp(b) & unzip_abs_xp < q_xp(b+1) & unzip_acc));
    mean_err_conf(b) = mean(unzip_conf(unzip_abs_xp > q_xp(b) & unzip_abs_xp < q_xp(b+1) & ~unzip_acc));
end

h4 = figure;
set(gcf, 'Position', [400 400 400 300]);
plot(mean_cor_conf, 'g', 'LineWidth', 3);
hold on
plot(mean_err_conf, 'r', 'LineWidth', 3);
set(gca, 'FontSize', 18, 'YLim', [0 1], 'XTick', [1 2 3], 'XTickLabel', {'Low', 'Med', 'High'}, 'XLim', [0.5 3.5]);
xlabel('X_{conf}','FontSize',18);
ylabel('Confidence','FontSize',18);
box off
if savePlots
    export_fig([figDir 'secondOrder_conf_xconf.png'], '-transparent', '-painters', h4)
end

% Figure for confidence as function of |xact|
unzip_abs_xa = abs(reshape(xa,1,length(theta).*N));
% bin into quantiles
q_xa = quantile(unzip_abs_xa, [0.33 0.66]);
q_xa = [0 q_xa Inf];
for b = 1:3
    mean_cor_conf(b) = mean(unzip_conf(unzip_abs_xa > q_xa(b) & unzip_abs_xa < q_xa(b+1) & unzip_acc));
    mean_err_conf(b) = mean(unzip_conf(unzip_abs_xa > q_xa(b) & unzip_abs_xa < q_xa(b+1) & ~unzip_acc));
end

h5 = figure;
set(gcf, 'Position', [400 400 400 300]);
plot(mean_cor_conf, 'g', 'LineWidth', 3);
hold on
plot(mean_err_conf, 'r', 'LineWidth', 3);
set(gca, 'FontSize', 18, 'YLim', [0 1], 'XTick', [1 2 3], 'XTickLabel', {'Low', 'Med', 'High'}, 'XLim', [0.5 3.5]);
xlabel('X_{act}','FontSize',18);
ylabel('Confidence','FontSize',18);
box off
if savePlots
    export_fig([figDir 'secondOrder_conf_xact.png'], '-transparent', '-painters', h5)
end

%% Show X-pattern as function of parameters rho and sigma_conf
theta = [0 0.032 0.064 0.128 0.256 0.512 1];
N = 10000;    % trials per level of coherence

%% 1) Sigma_conf
sigmaAct = 1;
sigmaConfVec = linspace(0.5, 1.5, 10);
rho = 0.4;

h6 = figure;
cols = [0.2 0.2 0.2; 0.5 0.5 0.5; 0.8 0.8 0.8];
set(gcf, 'Position', [400 400 800 300]);
displayIndex = 1;
for sigmai = 1:length(sigmaConfVec)
    sigmaConf = sigmaConfVec(sigmai);
    bigSigma = computeCov(sigmaAct, sigmaConf, rho);
    
    % Simulate x's, decisions and confidence for meta model
    mean_cor = nan(length(theta), N);
    xa = nan(length(theta), N);
    xp = nan(length(theta), N);
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
            xa(s,i) = r(1);
            xp(s,i) = r(2);
            
            if xa(s,i) > 0
                a(s,i) = 1;
                flip_a = 1;
            else
                a(s,i) = -1;
                flip_a = 0;
            end
            
            % Compute confidence
            mean_cor(s,i) = computeMetaConf(xp(s,i), flip_a, sigmaAct, sigmaConf, rho);
            
        end
        fprintf('Simulation completed for level %d... \n', s);
    end
    
    for s = 1:length(theta)
        acc(s) = mean(d(s,:) == a(s,:));
        conf_err(sigmai, s) = mean(mean_cor(s, d(s,:) ~= a(s,:)));
        conf_cor(sigmai, s) = mean(mean_cor(s, d(s,:) == a(s,:)));
    end
    slope_cor(sigmai) = conf_cor(sigmai,end) - conf_cor(sigmai,1);
    slope_err(sigmai) = conf_err(sigmai,end) - conf_err(sigmai,1);
    
    if sigmai == 1 | sigmai == length(sigmaConfVec) | sigmai == round(length(sigmaConfVec)./2)
        subplot(1,2,1);
        plot(theta, conf_cor(sigmai,:), 'LineWidth', 3, 'Color', cols(displayIndex,:));
        hold on
        plot(theta, conf_err(sigmai,:), 'LineWidth', 3, 'LineStyle', '--', 'Color', cols(displayIndex,:));
        displayIndex = displayIndex + 1;
    end
end

set(gca, 'FontSize', 18, 'YLim', [0 1], 'XLim', [0 max(theta)+0.1]);
xlabel('Stimulus strength', 'FontSize', 20);
ylabel('Confidence', 'FontSize', 20);
box off

subplot(1,2,2);
plot(sigmaConfVec, slope_cor, 'go ', 'MarkerSize', 8)
hold on
plot(sigmaConfVec, slope_err, 'ro ', 'MarkerSize', 8)
set(gca, 'FontSize', 18)
xlabel('\sigma_{conf}', 'FontSize', 20);
ylabel('Slope', 'FontSize', 20);
legend('Correct', 'Error', 'Location', 'SouthEast')
legend boxoff
box off

if savePlots
    export_fig([figDir 'secondOrder_Xpattern_sigmaConf.png'], '-transparent', '-painters', h6)
end

%% 2) Rho
sigmaAct = 1;
sigmaConf = 1;
rhoVec = linspace(0.1,0.9,10);

h7 = figure;
cols = [0.2 0.2 0.2; 0.5 0.5 0.5; 0.8 0.8 0.8];
set(gcf, 'Position', [400 400 800 300]);
displayIndex = 1;
for sigmai = 1:length(rhoVec)
    rho = rhoVec(sigmai);
    bigSigma = computeCov(sigmaAct, sigmaConf, rho);
    
    % Simulate x's, decisions and confidence for meta model
    mean_cor = nan(length(theta), N);
    xa = nan(length(theta), N);
    xp = nan(length(theta), N);
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
            xa(s,i) = r(1);
            xp(s,i) = r(2);
            
            if xa(s,i) > 0
                a(s,i) = 1;
                flip_a = 1;
            else
                a(s,i) = -1;
                flip_a = 0;
            end
            
            % Compute confidence
            mean_cor(s,i) = computeMetaConf(xp(s,i), flip_a, sigmaAct, sigmaConf, rho);
            
        end
        fprintf('Simulation completed for level %d... \n', s);
    end
    
    for s = 1:length(theta)
        acc(s) = mean(d(s,:) == a(s,:));
        conf_err(sigmai, s) = mean(mean_cor(s, d(s,:) ~= a(s,:)));
        conf_cor(sigmai, s) = mean(mean_cor(s, d(s,:) == a(s,:)));
    end
    slope_cor(sigmai) = conf_cor(sigmai,end) - conf_cor(sigmai,1);
    slope_err(sigmai) = conf_err(sigmai,end) - conf_err(sigmai,1);
    
    if sigmai == 1 | sigmai == length(rhoVec) | sigmai == round(length(rhoVec)./2)
        subplot(1,2,1);
        plot(theta, conf_cor(sigmai,:), 'LineWidth', 3, 'Color', cols(displayIndex,:));
        hold on
        plot(theta, conf_err(sigmai,:), 'LineWidth', 3, 'LineStyle', '--', 'Color', cols(displayIndex,:));
        displayIndex = displayIndex + 1;
    end
end

set(gca, 'FontSize', 18, 'YLim', [0 1], 'XLim', [0 max(theta)+0.1]);
xlabel('Stimulus strength', 'FontSize', 20);
ylabel('Confidence', 'FontSize', 20);
box off

subplot(1,2,2);
plot(rhoVec, slope_cor, 'go ', 'MarkerSize', 8)
hold on
plot(rhoVec, slope_err, 'ro ', 'MarkerSize', 8)
set(gca, 'FontSize', 18)
xlabel('\rho', 'FontSize', 20);
ylabel('Slope', 'FontSize', 20);
legend('Correct', 'Error', 'Location', 'SouthEast')
legend boxoff
box off

if savePlots
    export_fig([figDir 'secondOrder_Xpattern_rho.png'], '-transparent', '-painters', h7)
end