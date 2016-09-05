% Check equivalence between covariance and serial computation of confidence
%
% SF 2015

sigma_a = 0.5;
sigma_p = 0.75;
sigma_0 = 0.8;
a = 1;

xp = linspace(-2, 2, 1000);

conf1 = computeMetaConf_old(xp, a, sigma_0, sigma_a, sigma_p);
conf2 = computeMetaConf_alt(xp, a, sigma_0, sigma_a, sigma_p);

figure;
scatter(conf1, conf2);