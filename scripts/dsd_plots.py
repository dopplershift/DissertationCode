import matplotlib.pyplot as plt
import numpy as np
import scattering
import dsd
import scipy.constants as consts
import quantities as pq

d = np.linspace(0.01, 0.8, 100) * pq.centimeter
#d0_lut = (np.linspace(0.01, .5, 40).reshape(1, -1, 1).astype(np.float32)
#        * consts.centi)
#nr_lut = np.logspace(-2, 6, 110).reshape(1, 1,-1).astype(np.float32)

trials = [(0.2, 1e4, -0.8, 'k'), (0.4, 1e4, -0.8, 'r'), (0.1, 1e4, -0.8, 'b'),
          (0.2, 1e5, -0.8, 'r--'), (0.2, 1e3, -0.8, 'b--'),
          (0.2, 1e4, 0.1, 'r:'), (0.2, 1e4, 0.8, 'b:')]
for d0,Nr,nu,plot in trials:
    #dist = dsd.gamma(d, d0_lut, nr_lut, nu=-0.8)
    dist = dsd.gamma(d, d0 * pq.centimeter, Nr / pq.meter**3, nu)
    plt.semilogy(d, dist, plot,
            label=r'$d_0:%.1f\, N_r:%.0e\, \nu:%.1f$' % (d0,Nr,nu))

plt.legend(loc='lower right')
plt.title('DSD Comparison')
plt.xlabel('Diameter (%s)' % d.dimensionality.latex)
plt.ylabel('Number (%s)' % dist.dimensionality.latex)
plt.ylim(1e-2, None)
plt.grid()
plt.show()
