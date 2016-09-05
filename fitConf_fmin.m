function [params, out] = fitConf_fmin(data, Nsamp, sigmaAct, model, pArray)
% function [params, out] = fitConf_fmin(data, Nsamp, sigmaAct, model, pArray)
%
% Fit second-order model to data
% Requires data to be structured into 1 x N*2 array of rating counts, where
% N is the number of confidence ratings available to the subject and the
% first N bins contains the incorrect responses, and second N bins correct
% responses
%
% E.g. if there were 3 possible confidence levels available, a typical data
% array would look like [10 5 4 12 24 48] corresponding to:
%
% rating = 1, incorrect = 10 trials
% rating = 2, incorrect = 5 trials
% rating = 3, incorrect = 4 trials
% rating = 1, correct = 12 trials
% rating = 2, correct = 24 trials
% rating  = 3, correct = 48 trials
%
% sigmaAct is passed in as a fixed parameter and derived e.g. from the
% subject's type 1 d' (see parameter recovery script for an example of this)
% 
% Returns transformed parameter values for sigmaConf = exp(params(1)) and
% rho = 1./(1+exp(-params(2))
% 
% Steve Fleming 2016

% Calculate probabilities for ratings
nratings = length(data)./2;
totCounts = data(1:nratings) + data(nratings+1:end);
countsCDF = cumsum(totCounts)./sum(totCounts);

options = optimset('Display','Iter');
[params, negLogLik] = fminsearch(@fitfunc,pArray,options);
rng('default');

    function negLogLik = fitfunc(pArray)
        
        % Reset random number generator
        rng('default');
        
        % Distribute parameters
        switch model
            case 'sigma'
                sigmaConf = pArray(1);
                rho = 0.7;
            case 'sigma+rho'
                sigmaConf = exp(pArray(1));
                rho = 1./(1+exp(-pArray(2)));
        end
        
        a = [0 1];
        
        Sigma = computeCov(sigmaAct, sigmaConf, rho);
        r = mvnrnd([1 1], Sigma, Nsamp);
        xact = r(:,1);
        xconf = r(:,2);
        xpos = xconf(xact > 0)';
        xneg = xconf(xact <= 0)';
        
        % Compute distributions of confidence ratings conditional on
        % model samples for correct (a=1) and incorrect (a=0) responses
        conf_a0 = computeMetaConf(xneg, a(1), sigmaAct, sigmaConf, rho);
        conf_a1 = computeMetaConf(xpos, a(2),  sigmaAct,  sigmaConf, rho);
        
        % Specify thresholds from full distribution of confidence
        % irrespective of accuracy
        conf_full = [conf_a0 conf_a1];
        thresholds = quantile(conf_full, countsCDF);
        thresholds = [0 thresholds];
        
        % Bin according to threshold
        for b = 1:length(thresholds)-1
            bin_a0(b) = sum(conf_a0 > thresholds(b) & conf_a0 <= thresholds(b+1));
            bin_a1(b) = sum(conf_a1 > thresholds(b) & conf_a1 <= thresholds(b+1));
        end
        bin_a0(bin_a0 == 0) = 1;
        bin_a1(bin_a1 == 0) = 1;
        
        simCounts = [bin_a0 bin_a1];
        simProbs = simCounts./sum(simCounts); % create vector of normalised probabilities of confidence ratings for output
        
        % Get multinomial log-likelihood for this parameter setting against observed
        % data
        loglik = mnlogpdf(data, simProbs);
        negLogLik = -loglik;
        
        out.simProbs = simProbs;
        out.negLogLik = negLogLik;
    end
end