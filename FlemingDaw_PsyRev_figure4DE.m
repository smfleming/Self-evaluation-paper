% Simulations of how error awareness relates to second-order model
% parameters
%
% Steve Fleming 2015

clear all
close all

print_fig = 1;
addpath('~/Dropbox/Utils/graphics/export_fig/')
figDir = '~/Dropbox/Research/Metacognition/BN_model/selfSelf/figures/';

% Generate samples of Xact and Xconf with d = 1 (i.e. negative values of
% Xact are errors)
ntrials = 500;
mu = [1 1];
sigmaAct = 1;
sigmaConf = 1;
rho = 0.6;
bigSigma = computeCov(sigmaAct, sigmaConf, rho);
r = mvnrnd(mu, bigSigma, ntrials);
Xact = r(:,1);
Xconf = r(:,2);

% Step through trials, compute confidence conditional on Xconf and a
err = zeros(1,ntrials);   % intialise vector to store whether error occurs
err_aware = zeros(1,ntrials);   % intialise vector to store whether error is recognised or not
for i = 1:ntrials
    % compute action
    if Xact(i) > 0
        a = 1;
    else
        a = -1;
    end
    
    conf(i) = computeMetaConf(Xconf(i), a, sigmaAct, sigmaConf, rho);
    if a == -1
        err(i) = 1;
        if conf(i) < 0.5
            err_aware(i) = 1;
        end
    end
end

h1 = figure;
set(gcf, 'Position', [500 500 400 300])
plot(Xact(err==0), Xconf(err==0), 'o ');
hold on
plot(Xact(err==1), Xconf(err==1), 'o ');
plot(Xact(err_aware==1), Xconf(err_aware==1), 'o ');
set(gca, 'FontSize', 18);
xlabel('X_{act}');
ylabel('X_{conf}');
legend('Correct', 'Undetected error', 'Detected error', 'Location', 'SouthEast', 'FontSize', 14);
legend boxoff
box off
if print_fig
    export_fig([figDir 'error_awareness_scatter.png'], '-transparent', '-painters', h1)
end

% Loop over different parameter settings, compute p(detected errors), plot
% as heatmap
ntrials = 10000;
rho_vec = linspace(0.1,0.9,10);
sigmaConf_vec = linspace(0.5,1.5,10);
sigmaAct = 1;

for sigmai = 1:length(sigmaConf_vec)
    for rhoi = 1:length(rho_vec)
        bigSigma = computeCov(sigmaAct, sigmaConf_vec(sigmai), rho_vec(rhoi));
        r = mvnrnd(mu, bigSigma, ntrials);
        Xact = r(:,1);
        Xconf = r(:,2);
        err = zeros(1,ntrials);   % intialise vector to store whether error is recognised or not
        err_aware = zeros(1,ntrials);   % intialise vector to store whether error is recognised or not
        for i = 1:ntrials
            % compute action
            if Xact(i) > 0
                a = 1;
            else
                a = -1;
            end
            
            conf(i) = computeMetaConf(Xconf(i), a, sigmaAct, sigmaConf_vec(sigmai), rho_vec(rhoi));
            if a == -1
                err(i) = 1;
                if conf(i) < 0.5
                    err_aware(i) = 1;
                end
            end
        end
        p_detected_error(sigmai,rhoi) = sum(err_aware)./sum(err);
        fprintf('Simulation completed for sigmaConf = %d and rho = %d... \n', sigmaConf_vec(sigmai), rho_vec(rhoi));
    end
end

h2 = figure;
set(gcf, 'Position', [500 500 400 300])
imagesc(rho_vec, sigmaConf_vec, p_detected_error);
set(gca, 'FontSize', 18);
xlabel('\rho');
ylabel('\sigma_{conf}');
c = colorbar;
ylabel(c, 'P(detected error)')
box off
if print_fig
    export_fig([figDir 'error_awareness_heatmap.png'], '-transparent', '-painters', h2)
end
