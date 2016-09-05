function [mu_prod var_prod] = gaussian_product(mu_a, mu_b, var_a, var_b)
% Compute moments of products of two Gaussians
%
% SF 2014

mu_prod = (mu_a/var_a + mu_b/var_b)./(1/var_a + 1/var_b);
var_prod = 1/(1/var_a + 1/var_b);