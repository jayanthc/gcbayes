# get_mu_range.py
# Shows range of mu_l for given ranges of mu_s and d
#
# Created by Jayanth Chennamangalam

import sys
import argparse
import numpy as np
import matplotlib.pyplot as plt


parser = argparse.ArgumentParser(
    description='Show range of mu_l for given ranges of mu_s and d'
)
parser.add_argument('--cluster', '-c', type=str, help='Cluster name')
parser.add_argument('--mu_l_min', '-n', type=float, default=-2.0,
                    help='Minimum of mu_l prior (default: %(default)s)')
parser.add_argument('--mu_l_max', '-x', type=float, default=0.5,
                    help='Maximum of mu_l prior (defaut: %(default)s)')
args = parser.parse_args()

if args.cluster == 'ter5':
    # ter5
    mu_d = 5.5      # kpc
    sigma_d = 0.9
elif args.cluster == '47tuc':
    # 47tuc
    mu_d = 4.69     # kpc
    sigma_d = 0.17
elif args.cluster == 'm28':
    # m28
    mu_d = 5.5      # kpc
    sigma_d = 0.3
elif args.cluster == 'omegacen':
    # omega cen
    mu_d = 5.19     # kpc
    sigma_d = 0.075
else:
    sys.stderr.write('ERROR: Invalid cluster name!\n')
    sys.exit(1)

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
# mu_l_min = -2.0
# mu_l_max = +0.5

# max. priors for omega cen
# mu_l_min = -1.75
# mu_l_max = -1.65

mu_s = np.linspace(mu_s_min, mu_s_max, 100)
mu_l = np.zeros((100, 100))
for i in range(100):
    for j in range(100):
        mu_l[i, j] = mu_s[i] + 2 * np.log10(d[j])

plt.pcolormesh(d, mu_s, mu_l)
plt.contour(d, mu_s, mu_l,
            levels=np.linspace(args.mu_l_min, args.mu_l_max, 2),
            colors=['black'])
plt.xlabel('d (kpc)')
plt.ylabel('mu_s (mJy kpc^2)')
plt.show()
