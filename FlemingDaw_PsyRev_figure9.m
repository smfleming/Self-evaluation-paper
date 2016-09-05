% Simulated data for mismatch between generative and observation model
% parameters
%
% Steve Fleming 2014

clear all

addpath('~/Dropbox/Utils/graphics/export_fig/');
figDir = '~/Dropbox/Research/Metacognition/BN_model/selfSelf/figures/';
savePlots = 0;

N = 100000; % trials per simulation

%% Effect of observation sigma_act
gen_sigma_act = 1.5;
gen_sigma_conf = 1;
gen_rho = 0.6;
obs_sigma_act = [1 1.33 1.66 2];
obs_sigma_conf = 1;
obs_rho = 0.6;

for c = 1:length(obs_sigma_act)
    
    % Simulate x's, decisions and confidence for single route model
    conf = nan(1, N);
    d = nan(1, N);
    a = nan(1, N);
    
    for i = 1:N
        
        if rand > 0.5
            d(i) = 1;
        else
            d(i) = -1;
        end
        
        bigSigma = computeCov(gen_sigma_act, gen_sigma_conf, gen_rho);
        r = mvnrnd([d(i) d(i)], bigSigma, 1);
        xa = r(1);
        xp = r(2);
        
        if xa > 0
            a(i) = 1;
        else
            a(i) = -1;
        end
        
        % Compute self-confidence
        conf(i) = computeMetaConf(xp, a(i), obs_sigma_act(c), obs_sigma_conf, obs_rho);
        
    end
    
    confBins = linspace(0,1,11);
    for s = 1:length(confBins)-1
        mean_acc(c, s) = nanmean(d(conf > confBins(s) & conf <= confBins(s+1)) == a(conf > confBins(s) & conf <= confBins(s+1)));
        mean_conf(c, s) = nanmean(conf(conf > confBins(s) & conf <= confBins(s+1)));
    end
    disp(['Simulation complete for parameter ' num2str(c) ' out of ' num2str(length(obs_sigma_act))])
end

h1 = figure;
set(gcf, 'Position', [200 200 400 350])
map = colormap('copper');
points = linspace(min(obs_sigma_act), max(obs_sigma_act), length(map));
for i = 1:length(obs_sigma_act)
    ind = dsearchn(points', obs_sigma_act(i));
    plot(mean_conf(i,:), mean_acc(i,:), 'o-', 'LineWidth', 2, 'Color', map(ind,:));
    hold on
end
c = colorbar;
caxis([min(obs_sigma_act), max(obs_sigma_act)])
xlabel('Confidence', 'FontSize', 20);
ylabel('Performance','FontSize', 20);
ylabel(c,'\sigma_{act} (subject)','FontSize',16)
line([0 1], [0 1], 'Color', 'k', 'LineStyle', '--')
set(gca, 'XLim', [0 1], 'YLim', [0 1], 'FontSize', 18);
box off

if savePlots
    export_fig([figDir 'mismatchSim_sigmaAct.png'], '-transparent', '-painters', h1)
end

%% Effect of observation sigma_conf
gen_sigma_act = 1.5;
gen_sigma_conf = 1;
gen_rho = 0.6;
obs_sigma_act = 1.5;
obs_sigma_conf = [0.8 1 1.2 1.4];
obs_rho = 0.6;

for c = 1:length(obs_sigma_conf)
    
    % Simulate x's, decisions and confidence for single route model
    conf = nan(1, N);
    d = nan(1, N);
    a = nan(1, N);
    
    for i = 1:N
        
        if rand > 0.5
            d(i) = 1;
        else
            d(i) = -1;
        end
        
        bigSigma = computeCov(gen_sigma_act, gen_sigma_conf, gen_rho);
        r = mvnrnd([d(i) d(i)], bigSigma, 1);
        xa = r(1);
        xp = r(2);
        
        if xa > 0
            a(i) = 1;
        else
            a(i) = -1;
        end
        
        % Compute self-confidence
        conf(i) = computeMetaConf(xp, a(i), obs_sigma_act, obs_sigma_conf(c), obs_rho);
        
    end
    
    confBins = linspace(0,1,11);
    for s = 1:length(confBins)-1
        mean_acc(c, s) = nanmean(d(conf > confBins(s) & conf <= confBins(s+1)) == a(conf > confBins(s) & conf <= confBins(s+1)));
        mean_conf(c, s) = nanmean(conf(conf > confBins(s) & conf <= confBins(s+1)));
    end
    
     disp(['Simulation complete for parameter ' num2str(c) ' out of ' num2str(length(obs_sigma_conf))])
    
end

h2 = figure;
set(gcf, 'Position', [200 200 400 350])
map = colormap('copper');
points = linspace(min(obs_sigma_conf), max(obs_sigma_conf), length(map));
for i = 1:length(obs_sigma_conf)
    ind = dsearchn(points', obs_sigma_conf(i));
    plot(mean_conf(i,:), mean_acc(i,:), 'o-', 'LineWidth', 2, 'Color', map(ind,:));
    hold on
end
c = colorbar;
caxis([min(obs_sigma_conf), max(obs_sigma_conf)])
xlabel('Confidence', 'FontSize', 20);
ylabel('Performance','FontSize', 20);
ylabel(c,'\sigma_{conf} (subject)','FontSize',16)
line([0 1], [0 1], 'Color', 'k', 'LineStyle', '--')
set(gca, 'XLim', [0 1], 'YLim', [0 1], 'FontSize', 18);
box off

if savePlots
    export_fig([figDir 'mismatchSim_sigmaConf.png'], '-transparent', '-painters', h2)
end

%% Effect of observation sigma_p
gen_sigma_act = 1.5;
gen_sigma_conf = 1;
gen_rho = 0.6;
obs_sigma_act = 1.5;
obs_sigma_conf = 1;
obs_rho = [0.2 0.4 0.6 0.8];

for c = 1:length(obs_rho)
    
    % Simulate x's, decisions and confidence for single route model
    conf = nan(1, N);
    d = nan(1, N);
    a = nan(1, N);
    
    for i = 1:N
        
        if rand > 0.5
            d(i) = 1;
        else
            d(i) = -1;
        end
        
        bigSigma = computeCov(gen_sigma_act, gen_sigma_conf, gen_rho);
        r = mvnrnd([d(i) d(i)], bigSigma, 1);
        xa = r(1);
        xp = r(2);
        
        if xa > 0
            a(i) = 1;
        else
            a(i) = -1;
        end
        
        % Compute self-confidence
        conf(i) = computeMetaConf(xp, a(i), obs_sigma_act, obs_sigma_conf, obs_rho(c));
        
    end
    
    confBins = linspace(0,1,11);
    for s = 1:length(confBins)-1
        mean_acc(c, s) = nanmean(d(conf > confBins(s) & conf <= confBins(s+1)) == a(conf > confBins(s) & conf <= confBins(s+1)));
        mean_conf(c, s) = nanmean(conf(conf > confBins(s) & conf <= confBins(s+1)));
    end
    
     disp(['Simulation complete for parameter ' num2str(c) ' out of ' num2str(length(obs_rho))])
end

h3 = figure;
set(gcf, 'Position', [200 200 400 350])
map = colormap('copper');
points = linspace(min(obs_rho), max(obs_rho), length(map));
for i = 1:length(obs_rho)
    ind = dsearchn(points', obs_rho(i));
    plot(mean_conf(i,:), mean_acc(i,:), 'o-', 'LineWidth', 2, 'Color', map(ind,:));
    hold on
end
c = colorbar;
caxis([min(obs_rho), max(obs_rho)])
xlabel('Confidence', 'FontSize', 20);
ylabel('Performance','FontSize', 20);
ylabel(c,'\rho (subject)','FontSize',16)
line([0 1], [0 1], 'Color', 'k', 'LineStyle', '--')
set(gca, 'XLim', [0 1], 'YLim', [0 1], 'FontSize', 18);
box off

if savePlots
    export_fig([figDir 'mismatchSim_rho.png'], '-transparent', '-painters', h3)
end
