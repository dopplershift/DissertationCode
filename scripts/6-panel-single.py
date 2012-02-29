import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import ImageGrid
from disser.io import NetCDFRadarData
from disser.plots import PPIPlot, LabelGenerator, get_cmap
from disser import datatypes
import disser.plots.defaults as defaults

data = NetCDFRadarData('Sband_3600_20111010_210552.nc')

datatypes.TypePlotInfo[datatypes.Reflectivity].update(cmap=get_cmap('ScharfRefScaled'),
    norm=plt.Normalize(-20, 80))

moments = [data.fields.grab(moment, pol='H', source=src)
           for src in ('ts', 'average') for moment in (datatypes.Reflectivity,
           datatypes.DopplerVelocity, datatypes.SpectrumWidth)]

fig = plt.figure()
grid = ImageGrid(fig, (1, 1, 1), nrows_ncols = (2, 3), direction='row',
    share_all=True, axes_pad=0.6, aspect=True, cbar_mode='each',
    cbar_location='right', cbar_pad=0.15, cbar_size='10%')

for m, ax, cax, panel_label in zip(moments, grid, grid.cbar_axes, LabelGenerator('a')):
    ppi = PPIPlot(data.fields, var=m, ax=ax, rings=defaults.rings)
    panel_label.patch.set_boxstyle("round, pad=0., rounding_size=0.2")
    ax.add_artist(panel_label)
    ax.set_title(m)
    cbar = cax.colorbar(ppi.mesh)
    cbar.set_label_text(data.fields[m].dimensionality)
    cbar.cbar_axis.labelpad = -4

ax.set_ylim(0, 50)
ax.set_xlim(-20, 20)

plt.show()
