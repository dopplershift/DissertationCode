import matplotlib.pyplot as plt
from disser.io import NetCDFRadarData
from disser import datatypes
import disser.plots.defaults as defaults

data = NetCDFRadarData('Sband_3600_20111010_210552.nc')

moments = [data.fields.grab(moment, pol='H', source='ts')
           for moment in (datatypes.Power, datatypes.Reflectivity,
               datatypes.DopplerVelocity, datatypes.SpectrumWidth,
               datatypes.ZDR, datatypes.RhoHV, datatypes.KDP, datatypes.PhiDP)]

fig = plt.figure()

# Need some way to control the formatting of a datatype into a string
# Maybe using a "global" setting and a context maanger
grid = defaults.multipanel_cbar_each(fig, (2,4), moments, data)
ax = grid[0]

ax.xaxis.set_major_locator(plt.MultipleLocator(10))
ax.set_ylim(0, 50)
ax.set_xlim(-20, 20)

plt.show()
