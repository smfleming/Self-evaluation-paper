% Compare analytic and JAGS confidences across range of parameters

clear all
close all

addpath(genpath('~/Dropbox/Utils/matjags'));

% Set d
d = 1;
action = [0 1];
% Evaluate comparision for each parameter
sigma = 1;
sigma_a = 1.5;
sigma_p = 1;
xp_space = linspace(-2,2,20);

for acti = 1:2
    for xpi = 1:length(xp_space)
        
        % JAGS
        jags_conf(acti, xpi) = sampleMetaConf_singleTrial(xp_space(xpi), action(acti), sigma, sigma_a, sigma_p);
        
        % Analytic
        ana_conf(acti, xpi) = computeMetaConf_old(xp_space(xpi), action(acti), sigma, sigma_a, sigma_p);
        
    end
end

h = figure;
set(gcf, 'Position', [200 200 800 300]);
subplot(1,3,1)
plot(xp_space, jags_conf', 'LineWidth', 2);
legend({'a = L', 'a = R'}, 'Location', 'SouthEast');
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
xlabel('JAGS', 'FontSize', 14);
ylabel('Analytic', 'FontSize', 14);
axis square
set(gca, 'FontSize', 12);