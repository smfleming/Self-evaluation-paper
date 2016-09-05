% Simulate variance for sequential evidence accumulation, check calculation
% of 2sigma2 is correct through monte carlo simulation
%
% Steve Fleming 2016

ntrials = 10000;
sigma1 = 1.4;
sigma2 = 1.4;
varSeq = sigma1^2 + sigma2^2;

for i = 1:ntrials
    
    x1(i) = normrnd(1, sigma1);
    x2(i) = normrnd(1, sigma2);
    
    xtot(i) = x1(i) + x2(i);
    
end

figure;
subplot(1,3,1)
histogram(x1);
subplot(1,3,2)
histogram(x2);
subplot(1,3,3);
histogram(xtot);

mean(xtot)

% These should be the same
var(xtot)
varSeq