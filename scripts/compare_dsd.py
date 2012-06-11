import numpy as np
import quantities as pq
from disser import io
import dsd

num = 100
ind = 1453433

data = io.ModelData('/home/rmay/radar_sim_git/data/commas_wz_3600.nc')
mask = (data.qr > 1e-6) & (data.nr > 1e-3)
qr = data.qr[mask]
nr = data.nr[mask]
d = np.linspace(0.01 * pq.mm, 20 * pq.mm, 300).reshape(-1,1)
dist = dsd.constrained_gamma_from_moments(nr[ind:ind+num], qr[ind:ind+num], d)

test_qr = dsd.lwc(d, dist).simplified
test_nr = np.trapz(dist, axis=0, x=d).simplified
print qr[ind:ind+num], test_qr, np.allclose(qr[ind:ind+num], test_qr)
print nr[ind:ind+num], test_nr, np.allclose(nr[ind:ind+num], test_nr)
