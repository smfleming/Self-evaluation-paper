function out = genConf(Ntrials, Nsamp, thresholds, sigmaAct, sigmaConf, rho)
% Generate distributions of confidence ratings for 2-choice discrimination
% task, constant difficutly
%
% Generate confidence distributions for correct and error responses
%
% Steve Fleming 2016

a = [0 1];

% First sample from generative model of actions and internal states
Sigma = computeCov(sigmaAct, sigmaConf, rho);
r = mvnrnd([1 1], Sigma, Nsamp);
xact = r(:,1);
xconf = r(:,2);
xpos = xconf(xact > 0)';
xneg = xconf(xact <= 0)';

% Then compute distributions of confidence ratings conditional on
% model response
conf_a0 = computeMetaConf(xneg, a(1), sigmaAct, sigmaConf, rho);
conf_a1 = computeMetaConf(xpos, a(2),  sigmaAct,  sigmaConf, rho);

% Bin according to threshold and normalise
for b = 1:length(thresholds)-1
    bin_a0(b) = sum(conf_a0 > thresholds(b) & conf_a0 <= thresholds(b+1));
    bin_a1(b) = sum(conf_a1 > thresholds(b) & conf_a1 <= thresholds(b+1));
end

% Create single multinomial vector of length Nbins x choice x accuracy
simCounts = [bin_a0 bin_a1];
simProbs = simCounts./sum(simCounts);
out.allCounts = round(simProbs.*Ntrials);
