% murange.m
% Shows range of mu_l for given ranges of mu_s and d
%
% Created by Jayanth Chennamangalam

%% ter5
%mu_d = 5.5;     % kpc
%sigma_d = 0.9;
% 47tuc
%mu_d = 4.69;     % kpc
%sigma_d = 0.17;
% m28
mu_d = 5.5;     % kpc
sigma_d = 0.3;

d = linspace((mu_d - 3 * sigma_d), (mu_d + 3 * sigma_d));

mu_s_min = -6.0;
mu_s_max = 2.0;
% Boyles' et al. priors
%mu_l_min = -1.19;
%mu_l_max = -1.04;
% my wide priors
%mu_l_min = -3.19;
%mu_l_max = +2.04;
% bagchi et al. range
mu_l_min = -2.0;
mu_l_max = +0.5;

mu_s = linspace(mu_s_min, mu_s_max);
mu_l = zeros(100, 100);
for i = 1: 100
    for j = 1: 100
        mu_l(i, j) = mu_s(i) + 2 * log10(d(j));
    end
end

figure(1);
imagesc(d, mu_s, mu_l);
set(gca,'YDir','normal');
hold on;
[c, h] = contour(d, mu_s, mu_l, [mu_l_min mu_l_min; mu_l_max mu_l_max], 'k');
set(h,'ShowText','on');
xlabel('d (kpc)');
ylabel('\mu_s (mJy kpc^2)');
hold off;

