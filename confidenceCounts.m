function [errCounts, corCounts, sigma] = confidenceCounts(acc, conf, nratings)
% Take vectors of discrete confidence and accuracy and output error and
% correct confidence counts plus sigma from neutral-criterion SDT model
%
% SF 2016

for b = 1:nratings
    corCounts(b) = sum(acc(conf == b) == 1);
    errCounts(b) = sum(acc(conf == b) == 0);
end

% calculate sigma from simple SDT collapsing over choice
nI = sum(errCounts);
nC = sum(corCounts);
nTot = nI + nC;
d1 = norminv(nC./nTot) - norminv(nI./nTot);
sigma = 2/d1;
