function conf = computeFirstOrderConf(xconf, a, sigma_conf)
% function conf = computeFirstOrderConf(xconf, a, sigma_act, sigma_conf, rho)
% Compute confidence direct from xconf without inferring xact
%
% SF 2014

% Initialise
dhat = [-1 1];

% Conf = p(a=d)
if a == 1
    conf = normpdf(xconf,dhat(2),sigma_conf)./(normpdf(xconf,dhat(1),sigma_conf) + normpdf(xconf,dhat(2),sigma_conf));
else
    conf = normpdf(xconf,dhat(1),sigma_conf)./(normpdf(xconf,dhat(1),sigma_conf) + normpdf(xconf,dhat(2),sigma_conf));
end