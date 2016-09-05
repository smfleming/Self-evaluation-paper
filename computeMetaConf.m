function conf = computeMetaConf(xp, a, sigma_act, sigma_conf, rho)
% function conf = computeMetaConf(xp, a, sigma_act, sigma_conf, rho)
% Compute confidence from xp, a and covariance parameters
% If xp is a vector, gives the corresponding distribution of confidence
%
% SF 2014

% Initialise
dhat = [-1 1];
mu_x_xp_dhat = zeros(2,length(xp));
var_x_xp_dhat = zeros(1,length(xp));
rho_vec = repmat(rho,1,length(xp));
sigA_vec = repmat(sigma_act,1,length(xp));
sigP_vec = repmat(sigma_conf,1,length(xp));

Tol = 10e-4;

for dhati = 1:2
    
    dhat_vec = repmat(dhat(dhati),1,length(xp));
    
    mu_x_xp_dhat(dhati,:) = dhat_vec + (sigA_vec./sigP_vec).*rho_vec.*(xp-dhat_vec);
    var_x_xp_dhat = (1-rho_vec.^2).*sigA_vec.^2;
    
    if a == 1
        p_a_dhat_xp(dhati,:) = 1-normcdf(zeros(1,length(xp)), mu_x_xp_dhat(dhati,:), sqrt(var_x_xp_dhat));
    else
        p_a_dhat_xp(dhati,:) = normcdf(zeros(1,length(xp)), mu_x_xp_dhat(dhati,:), sqrt(var_x_xp_dhat));
    end
    
    lik_d(dhati,:) = normpdf(xp, dhat_vec, sigP_vec);
    
end

% Avoid underflow of probabilities
p_a_dhat_xp(p_a_dhat_xp < Tol) = Tol;
lik_d(lik_d < Tol) = Tol;

lik_d = lik_d./repmat(sum(lik_d,1),2,1);
p_dhat_xp_a = p_a_dhat_xp.*lik_d;
p_dhat_xp_a = p_dhat_xp_a./repmat(sum(p_dhat_xp_a,1),2,1);

% Conf = p(a=d)
if a == 1
    conf = p_dhat_xp_a(2,:);
else
    conf = p_dhat_xp_a(1,:);
end