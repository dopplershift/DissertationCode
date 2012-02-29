import numpy as np
import matplotlib.pyplot as plt
from disser.io import NetCDFRadarData
from disser.plots import PPIPlot, get_cmap
from disser import datatypes

data = NetCDFRadarData('Sband_3600_20111010_210552.nc')

rings = np.arange(10, 60, 10)

datatypes.TypePlotInfo.set_defaults(cmap=get_cmap('Carbone42'))
datatypes.TypePlotInfo[datatypes.Reflectivity] = dict(
    norm=plt.Normalize(-30, 80), cmap=get_cmap('Carbone42'))

fig = plt.figure()
ax = fig.add_subplot(1, 1, 1)
ax.set_aspect('equal', 'datalim')
varinfo = data.fields.grab(datatypes.Reflectivity, pol='V', source='ts')
ppi = PPIPlot(data.fields, var=varinfo, ax=ax, rings=rings,
    labels=[])
plt.show()
