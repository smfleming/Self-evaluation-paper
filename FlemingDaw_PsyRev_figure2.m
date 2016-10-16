% Make example plots of relationship between Xact and Xconf as function of
% parameters of Sigma
%
% Fig 2A, B, C

addpath('~/Dropbox/Utils/graphics/export_fig/')
figDir = '~/Dropbox/Research/Metacognition/BN_model/Self-evaluation-paper/figures/';

% how many points
nsamp = 500;

% Fig 1D (low rho)
mu = [1 1];
sigmaAct = 1;
sigmaConf = 1;
rho = 0.2;
bigSigma = computeCov(sigmaAct, sigmaConf, rho);
r = mvnrnd(mu, bigSigma, nsamp);

h1 = figure;
set(gcf, 'Position', [400 400 400 300]);
scatterhist(r(:,1), r(:,2), 'Kernel', 'on');
box off
xlabel('X_{act}');
ylabel('X_{conf}');
text(-3, 6, ['\sigma_{act} = ' num2str(sigmaAct)], 'FontSize', 14)
text(-3, 5, ['\sigma_{conf} = ' num2str(sigmaConf)], 'FontSize', 14)
text(-3, 4, ['\rho = ' num2str(rho)], 'FontSize', 14)
set(gca, 'FontSize', 16, 'XLim', [-4 6], 'YLim', [-4 6]);
axis square
export_fig([figDir 'XactXconf_scatter_1.png'], '-transparent', '-painters', h1)

% Fig 1E (high rho, difference in sigmas)
mu = [1 1];
sigmaAct = 1;
sigmaConf = 1.5;
rho = 0.6;
bigSigma = computeCov(sigmaAct, sigmaConf, rho);
r = mvnrnd(mu, bigSigma, nsamp);

figure;
set(gcf, 'Position', [400 400 400 300]);
h2 = scatterhist(r(:,1), r(:,2), 'Kernel', 'on');
box off
xlabel('X_{act}');
ylabel('X_{conf}');
text(-5, 6.6, ['\sigma_{act} = ' num2str(sigmaAct)], 'FontSize', 14)
text(-5, 5.3, ['\sigma_{conf} = ' num2str(sigmaConf)], 'FontSize', 14)
text(-5, 4.0, ['\rho = ' num2str(rho)], 'FontSize', 14)
set(gca, 'FontSize', 16, 'XLim', [-6 8], 'YLim', [-6 8]);
axis square
export_fig([figDir 'XactXconf_scatter_2.png'], '-transparent', '-painters', h2)

% Fig 1F (first-order case)
mu = [1 1];
sigmaAct = 1;
sigmaConf = 1;
rho = 1;
bigSigma = computeCov(sigmaAct, sigmaConf, rho);
r = mvnrnd(mu, bigSigma, nsamp);

figure;
set(gcf, 'Position', [400 400 400 300]);
h3 = scatterhist(r(:,1), r(:,2), 'Kernel', 'on');
box off
xlabel('X_{act}');
ylabel('X_{conf}');
text(-3, 6, ['\sigma_{act} = ' num2str(sigmaAct)], 'FontSize', 14)
text(-3, 5, ['\sigma_{conf} = ' num2str(sigmaConf)], 'FontSize', 14)
text(-3, 4, ['\rho = ' num2str(rho)], 'FontSize', 14)
set(gca, 'FontSize', 16, 'XLim', [-4 6], 'YLim', [-4 6]);
axis square
export_fig([figDir 'XactXconf_scatter_3.png'], '-transparent', '-painters', h3)