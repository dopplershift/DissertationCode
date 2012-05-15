import matplotlib.pyplot as plt
import numpy as np
import scattering
import dsd
import scipy.constants as consts
import quantities as pq

d = np.linspace(0.01, 8, 100) * pq.mm
#d0_lut = (np.linspace(0.01, .5, 40).reshape(1, -1, 1).astype(np.float32)
#        * consts.centi)
#nr_lut = np.logspace(-2, 6, 110).reshape(1, 1,-1).astype(np.float32)

trials = [(5 / pq.mm, 1e4 / pq.m**3, 'k'), (10 / pq.mm, 1e4 / pq.m**3, 'r'),
          (15 / pq.mm, 1e4 / pq.m**3, 'b'), (20 / pq.mm, 1e4 / pq.m**3, 'g'),
          (5 / pq.mm, 1e5 / pq.m**3, 'k--'), (10 / pq.mm, 1e5 / pq.m**3, 'r--'),
          (15 / pq.mm, 1e5 / pq.m**3, 'b--'), (20 / pq.mm, 1e5 / pq.m**3, 'g--')]
for lam,Nr,plot in trials:
    #dist = dsd.gamma(d, d0_lut, nr_lut, nu=-0.8)
    nu = dsd.constrained_gamma_shape(lam)
    dist = dsd.modified_gamma(d, lam, Nr, nu)
    plt.semilogy(d, dist, plot,
            label=r'$\Lambda:%.1f\, N_r:%.0e\, \nu:%.1f$' % (lam,Nr,nu))
    print dsd.lwc(d, dist).simplified, dsd.constrained_gamma_shape(lam)

plt.legend(loc='lower right')
plt.title('DSD Comparison')
plt.xlabel('Diameter (%s)' % d.dimensionality.latex)
plt.ylabel('Number (%s)' % dist.dimensionality.latex)
plt.ylim(1e1, None)
plt.grid()
plt.show()
