# get_mu_range.py
# Shows range of mu_l for given ranges of mu_s and d
#
# Created by Jayanth Chennamangalam

import numpy as np
import matplotlib.pyplot as plt


# ter5
# mu_d = 5.5      # kpc
# sigma_d = 0.9
# 47tuc
# mu_d = 4.69      # kpc
# sigma_d = 0.17
# m28
mu_d = 5.5      # kpc
sigma_d = 0.3

d = np.linspace((mu_d - 3 * sigma_d), (mu_d + 3 * sigma_d), 100)

mu_s_min = -6.0
mu_s_max = 2.0
# Boyles' et al. priors
# mu_l_min = -1.19
# mu_l_max = -1.04
# my wide priors
# mu_l_min = -3.19
# mu_l_max = +2.04
# bagchi et al. range
mu_l_min = -2.0
mu_l_max = +0.5

mu_s = np.linspace(mu_s_min, mu_s_max, 100)
mu_l = np.zeros((100, 100))
for i in range(100):
    for j in range(100):
        mu_l[i, j] = mu_s[i] + 2 * np.log10(d[j])


plt.pcolormesh(d, mu_s, mu_l)
plt.contour(d, mu_s, mu_l, colors=['black'])
plt.xlabel('d (kpc)')
plt.ylabel('mu_s (mJy kpc^2)')
plt.show()