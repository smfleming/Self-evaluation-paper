% Test product of three gaussians

base = linspace(-4,4,500);

mu_a = -0.7;
mu_b = 1.5;
mu_c = 2;
var_a = 0.3;
var_b = 0.7;
var_c = 0.2;

% Check
pdf_a = normpdf(base, mu_a, sqrt(var_a));
pdf_b = normpdf(base, mu_b, sqrt(var_b));
pdf_c = normpdf(base, mu_c, sqrt(var_c));
pdf_prod = pdf_a.*pdf_b.*pdf_c;

% Compute moments
[mu_prod_1 var_prod_1] = gaussian_product(mu_a, mu_b, var_a, var_b);
[mu_prod var_prod] = gaussian_product(mu_prod_1, mu_c, var_prod_1, var_c);
pdf_prod_exact = normpdf(base, mu_prod + 0.1, sqrt(var_prod));

figure;
plot(base, pdf_a./sum(pdf_a));
hold on
plot(base, pdf_b./sum(pdf_b), 'r');
plot(base, pdf_c./sum(pdf_c), 'g');
plot(base, pdf_prod./sum(pdf_prod), 'k');
plot(base, pdf_prod_exact./sum(pdf_prod_exact), 'k--')