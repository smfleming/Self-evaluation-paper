% Demonstrate fitting of second-order and same-order model to simulated
% data

Ntrials = 1000;
Nsamp = 100000;
thresholds = linspace(0, 1, 5);   % pre-specify thresholds for discrete ratings, you would get this from actual subject data
nratings = length(thresholds)-1;

%% Fit the data (second-order inference)
% Generate some simulated data
sigma = 1.5;
out = genConf(Ntrials, Nsamp, thresholds, sigma, 1.5, 0.7);

model = 'sigma+rho';
pArray = [log(2) log(0.7./(1-0.7))]; % initialise parameters
[params fit] = fitConf_fmin(out.allCounts+1, Nsamp, sigma, model, pArray);
% [params fit] = fitConf_gridsearch(out.allCounts+1, Nsamp, sigma, model);

figure;
plot(1:nratings, fit.simProbs(1:nratings), 'LineStyle', '--', 'LineWidth', 2);
hold on
plot(nratings+1:nratings*2, fit.simProbs(nratings+1:nratings*2), 'LineWidth', 2);
plot(out.allCounts./sum(out.allCounts), 'o ', 'MarkerSize', 8);
