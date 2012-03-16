import matplotlib.pyplot as plt
from disser.io import NetCDFRadarData
from disser.plots import get_cmap
from disser import datatypes
import disser.plots.defaults as defaults

data = NetCDFRadarData('Sband_3600_20111010_210552.nc')

datatypes.TypePlotInfo[datatypes.Reflectivity].update(cmap=get_cmap('ScharfRefScaled'),
    norm=plt.Normalize(-20, 80))

moments = [data.fields.grab(moment, pol='H', source=src)
           for src in ('ts', 'average') for moment in (datatypes.Reflectivity,
           datatypes.DopplerVelocity, datatypes.SpectrumWidth)]

fig = plt.figure()

grid = defaults.multipanel_cbar_column(fig, (2,3), moments, data)
ax = grid[0]

ax.xaxis.set_major_locator(plt.MultipleLocator(10))
ax.set_ylim(0, 50)
ax.set_xlim(-20, 20)

plt.show()
