% Simulate changes in metacognitive accuracy under second-order model
%
% Requires toolbox for fitting meta-d' from http://www.columbia.edu/~bsm2105/type2sdt/
%
% Steve Fleming 2015

clear all
close all

addpath('~/Dropbox/Utils/graphics/export_fig/');
addpath('~/Dropbox/Utils/meta_d/');
figDir = '~/Dropbox/Research/Metacognition/BN_model/Self-evaluation-paper/figures/';

savePlots = 0;

%% Sigma_conf simulation
sigma_act = 1;
sigma_conf = [0.5 1 1.5];
rho = 0.5;
thresh = [0.2 0.4 0.6 0.8];
nratings = 5;
ntrials = 10000;
conf_criteria = linspace(0,1,20);

for sigma_pi = 1:length(sigma_conf)
    
    [nR_S1 nR_S2 data] = metaModelFit_generate_data(ntrials, nratings, sigma_act, sigma_conf(sigma_pi), rho, thresh);
    
    conf_criteria = linspace(0,1,20);
    
    k=length(conf_criteria);
    for j = 1:length(conf_criteria)
        H(sigma_pi, k) = sum(data.conf(data.correct == 1) >= conf_criteria(j))./sum(data.correct);
        FA(sigma_pi, k) = sum(data.conf(data.correct == 0) >= conf_criteria(j))./sum(~data.correct);
        k=k-1;
    end
    
    dprime(sigma_pi) = 2./sigma_act;
end

h1 = figure;
set(gcf, 'Position', [200 200 700 250]);

subplot(1,2,1);
pos = get(gca, 'Position');
pos(2) = 0.2;
pos(4) = 0.7;
set(gca, 'Position', pos);
map = colormap('copper');
points = linspace(min(sigma_conf), max(sigma_conf), length(map));
for i = 1:size(H,1)
    ind = dsearchn(points', sigma_conf(i));
    plot(FA(i,:), H(i,:), 'LineWidth', 2, 'Color', map(ind,:));
    hold on
end
set(gca, 'FontSize', 14)
c = colorbar;
caxis([min(sigma_conf), max(sigma_conf)])
xlabel('P(confidence|incorrect)')
ylabel('P(confidence|correct)')
ylabel(c,'\sigma_{conf}','FontSize',14)
text(0.5, 0.3, ['\rho = ' num2str(rho)],'FontSize',14)
text(0.5, 0.2, ['\sigma_{act} = ' num2str(sigma_act)],'FontSize',14)

%% Rho simulation
% Keep sum of variances constant
sigma_act = 1;
sigma_conf = 1;
rho = [0.1 0.5 1];
thresh = [0.2 0.4 0.6 0.8];
nratings = 5;
ntrials = 10000;
conf_criteria = linspace(0,1,20);

for sigma_pi = 1:length(rho)
    
    [nR_S1 nR_S2 data] = metaModelFit_generate_data(ntrials, nratings, sigma_act, sigma_conf, rho(sigma_pi), thresh);
    
    conf_criteria = linspace(0,1,20);
    
    k=length(conf_criteria);
    for j = 1:length(conf_criteria)
        H(sigma_pi, k) = sum(data.conf(data.correct == 1) >= conf_criteria(j))./sum(data.correct);
        FA(sigma_pi, k) = sum(data.conf(data.correct == 0) >= conf_criteria(j))./sum(~data.correct);
        k=k-1;
    end
    
    dprime(sigma_pi) = 2./sigma_act;
end

subplot(1,2,2);
pos = get(gca, 'Position');
pos(2) = 0.2;
pos(4) = 0.7;
set(gca, 'Position', pos);
map = colormap('copper');
points = linspace(min(rho), max(rho), length(map));
for i = 1:size(H,1)
    ind = dsearchn(points', rho(i));
    plot(FA(i,:), H(i,:), 'LineWidth', 2, 'Color', map(ind,:));
    hold on
end
set(gca, 'FontSize', 14)
c = colorbar;
caxis([min(rho), max(rho)])
xlabel('P(confidence|incorrect)')
ylabel('P(confidence|correct)')
ylabel(c,'\rho','FontSize',14)
text(0.5, 0.3, ['\sigma_{conf} = ' num2str(sigma_conf)],'FontSize',14)
text(0.5, 0.2, ['\sigma_{act} = ' num2str(sigma_act)],'FontSize',14)

if savePlots
    export_fig([figDir 'type2ROC.png'], '-transparent', '-painters', h1)
end

%% Meta-d vs d simulation
rho = 0.5;
thresh = [0.2 0.4 0.6 0.8];
nratings = 5;
ntrials = 1000;
conf_criteria = linspace(0,1,20);
nsim = 100;

for n = 1:nsim
    
    sigma_conf(n) = rand + 1.5;
    sigma_act(n) = rand + 1.5;
    sigma_ratio(n) = sigma_conf(n)./sigma_act(n);
    [nR_S1 nR_S2 data] = metaModelFit_generate_data(ntrials, nratings, sigma_act(n), sigma_conf(n), rho, thresh);
    
    % Obtain meta-d' using MLE fit
    padFactor = 1/(2*nratings);
    nR_S1 = nR_S1 + padFactor;
    nR_S2 = nR_S2 + padFactor;
    fit = fit_meta_d_MLE(nR_S1, nR_S2);
    metad(n) = fit.meta_da;
    dprime(n) = fit.da;
    
    % Compute proportion of detected errors
    err = (data.d  > 0) ~= data.a;
    p_error_detect(n) = sum(err & data.conf < 0.5)./sum(err);
    
    disp([num2str(n) ' simulations complete'])
end

h2 = figure;
set(gcf, 'Position', [200 200 800 300]);
subplot(1,2,1)
scatter(dprime, metad, [], sigma_ratio)
set(gca, 'XLim', [0 2], 'YLim', [0 2], 'FontSize', 16);
line([0 2], [0 2], 'LineStyle', '--', 'Color', 'k')
c = colorbar
xlabel('d''')
ylabel('meta-d''')
ylabel(c, '\sigma_{conf} / \sigma_{act}')
axis square
box off

%% Mratio against proportion of recognised errors
subplot(1,2,2)
mratio = metad./dprime;
plot(log(mratio), p_error_detect, 'kx ', 'MarkerSize', 8);
line([0 0], [0 0.5], 'Color', 'k', 'LineStyle', '--')
set(gca, 'FontSize', 16);
xlabel('Log meta-d''/d''');
ylabel('Proportion of detected errors')
axis square
box off

if savePlots
    export_fig([figDir 'dprime_metadprime.png'], '-transparent', '-painters', h2)
end

