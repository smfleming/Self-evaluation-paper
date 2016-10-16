% Plots of how confidence relates to internal state in first-order and
% second-order models
%
% SF 2015

clear all
close all

savePlots = 1;
addpath(genpath('~/Dropbox/Utils/graphics/'));
figDir = '~/Dropbox/Research/Metacognition/BN_model/Self-evaluation-paper/figures/';

%% First-order model confidence
base = linspace(-1, 1, 500);
sigma_act = [0.5 0.75];
colors = [0.2 0.2 0.2; 0.8 0.8 0.8];

for sigmai = 1:length(sigma_act)
    for i = 1:length(base)
        
        conf(sigmai, i) = computeFirstOrderConf(base(i), sign(base(i)), sigma_act(sigmai));
        pd(sigmai, i) = normpdf(base(i),1,sigma_act(sigmai))./(normpdf(base(i),-1,sigma_act(sigmai)) + normpdf(base(i),1,sigma_act(sigmai)));
        
    end
end

h1 = figure;
set(gcf, 'Position', [500 500 400 300])
for sigmai = 1:length(sigma_act)    % run loop twice so that legend colours work out
    hold on
    plot(base(251:end), conf(sigmai,251:end), 'LineWidth', 2, 'Color', colors(sigmai,:))
    
end
for sigmai = 1:length(sigma_act)
    plot(base(1:250), conf(sigmai,1:250), 'LineWidth', 2, 'LineStyle', '--', 'Color', colors(sigmai,:))
end
set(gca, 'YLim',[0 1], 'XTick',[-1 0 1], 'YTick', [0 0.5 1], 'FontSize', 14, 'Position', [0.1 0.2 0.75 0.75]);
xlabel('X_{conf}', 'FontSize', 18)
ylabel('Confidence', 'FontSize', 18)
line([-1 1],[0.5 0.5], 'LineStyle', '--', 'Color', 'k')
l1 = legend({num2str(sigma_act(1)), num2str(sigma_act(2))}, 'Location', 'SouthEast')
pos = get(l1, 'Position');
pos(1) = 0.65; pos(2) = 0.6;
set(l1, 'Position', pos);
text(0.5, 0.6, '\sigma', 'FontSize', 16);
legend boxoff
% Text for additional legend
line([0.7 0.95],[0.4 0.4], 'LineStyle', '--', 'Color', [0 0 0])
line([0.7 0.95],[0.32 0.32], 'LineStyle', '-', 'Color', [0 0 0])
text(0.6, 0.4, 'L', 'FontSize', 14);
text(0.6, 0.32, 'R', 'FontSize', 14);
box off
axis square

if savePlots
    export_fig([figDir 'firstOrder_conf_x.png'], '-transparent', '-painters', h1)
end

%% Post-decisional model confidence
sigma_conf = [0.5 0.75];
for sigmai = 1:length(sigma_conf)
    for i = 1:length(base)
        
        conf_R(sigmai, i) = computeFirstOrderConf(base(i), 1, sqrt(2.*sigma_conf(sigmai).^2));
        conf_L(sigmai, i) = computeFirstOrderConf(base(i), -1, sqrt(2.*sigma_conf(sigmai).^2));
        
    end
end

h2 = figure;
set(gcf, 'Position', [500 500 400 300])
for sigmai = 1:length(sigma_conf)
    hold on
    plot(base, conf_R(sigmai,:), 'LineWidth', 2, 'Color', colors(sigmai,:))
end
for sigmai = 1:length(sigma_conf)
    hold on
    plot(base, conf_L(sigmai,:), 'LineWidth', 2, 'LineStyle', '--', 'Color', colors(sigmai,:))
end
set(gca, 'YLim',[0 1], 'XTick',[-1 0 1], 'YTick', [0 0.5 1], 'FontSize', 14, 'Position', [0.1 0.2 0.75 0.75]);
xlabel('X_{conf}', 'FontSize', 18)
ylabel('Confidence', 'FontSize', 18)
line([-1 1],[0.5 0.5], 'LineStyle', '--', 'Color', 'k')

% Legend
l1 = legend({num2str(sigma_conf(1)), num2str(sigma_conf(2))}, 'Location', 'East')
pos = get(l1, 'Position');
pos(1) = 0.65; pos(2) = 0.6;
set(l1, 'Position', pos);
text(0.5, 0.6, '\sigma', 'FontSize', 16);
legend boxoff

% Text for additional legend
line([0.7 0.95],[0.4 0.4], 'LineStyle', '--', 'Color', [0 0 0])
line([0.7 0.95],[0.32 0.32], 'LineStyle', '-', 'Color', [0 0 0])
text(0.6, 0.4, 'L', 'FontSize', 14);
text(0.6, 0.32, 'R', 'FontSize', 14);
box off
axis square

if savePlots
    export_fig([figDir 'postDecisional_conf_x.png'], '-transparent', '-painters', h2)
end

%% Second-order model confidence
sigma_conf = [0.5 0.75];
sigma_act = 1;
rho = 0.6;
for sigmai = 1:length(sigma_conf)
    for i = 1:length(base)
        
        conf_R(sigmai, i) = computeMetaConf(base(i), 1, sigma_act, sigma_conf(sigmai), rho);
        conf_L(sigmai, i) = computeMetaConf(base(i), -1, sigma_act, sigma_conf(sigmai), rho);
        
    end
end

h3 = figure;
set(gcf, 'Position', [500 500 400 300])
for sigmai = 1:length(sigma_conf)
    hold on
    plot(base, conf_R(sigmai,:), 'LineWidth', 2, 'Color', colors(sigmai,:))
end
for sigmai = 1:length(sigma_conf)
    hold on
    plot(base, conf_L(sigmai,:), 'LineWidth', 2, 'LineStyle', '--', 'Color', colors(sigmai,:))
end
set(gca, 'YLim',[0 1], 'XTick',[-1 0 1], 'YTick', [0 0.5 1], 'FontSize', 14, 'Position', [0.1 0.2 0.75 0.75]);
xlabel('X_{conf}', 'FontSize', 18)
ylabel('Confidence', 'FontSize', 18)
line([-1 1],[0.5 0.5], 'LineStyle', '--', 'Color', 'k')

% Legend
l1 = legend({num2str(sigma_conf(1)), num2str(sigma_conf(2))}, 'Location', 'East')
pos = get(l1, 'Position');
pos(1) = 0.65; pos(2) = 0.7;
set(l1, 'Position', pos);
text(0.35, 0.72, '\sigma_{conf} ', 'FontSize', 16);
legend boxoff

% Text for additional legend
line([0.7 0.95],[0.4 0.4], 'LineStyle', '--', 'Color', [0 0 0])
line([0.7 0.95],[0.32 0.32], 'LineStyle', '-', 'Color', [0 0 0])
text(0.6, 0.4, 'L', 'FontSize', 14);
text(0.6, 0.32, 'R', 'FontSize', 14);
box off
axis square

if savePlots
    export_fig([figDir 'secondOrder_conf_x.png'], '-transparent', '-painters', h3)
end

