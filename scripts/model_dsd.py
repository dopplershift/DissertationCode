import matplotlib.pyplot as plt
from disser import io

data = io.ModelData('/home/rmay/radar_sim_git/data/commas_wz_3600.nc')

mask = (data.qr > 1e-6) & (data.nr > 1e-3)

qr = data.qr[mask]
nr = data.nr[mask]

plt.hexbin(qr, nr, bins='log')
plt.xlabel(r'$q_r$ (%s)' % qr.dimensionality.latex)
plt.ylabel(r'$N_r$ (%s)' % nr.dimensionality.latex)
plt.colorbar()

plt.savefig('model-dsd-hexbin.png')
plt.show()
