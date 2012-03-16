import numpy as np
from functools import partial
import quantities as pq
pq.markup.format_units_latex = partial(pq.markup.format_units_latex,
    font='mathsf')

rings = np.arange(10, 60, 10) * pq.km

from .ctables import get_cmap
from .basic import LabelGenerator
from .ppi import PPIPlot
from disser import datatypes
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import ImageGrid

# Default colormaps and norms get set up on import
datatypes.TypePlotInfo.set_defaults(cmap=get_cmap('Carbone42'))
datatypes.TypePlotInfo[datatypes.Reflectivity].update(
    norm=plt.Normalize(-20, 70))
datatypes.TypePlotInfo[datatypes.DopplerVelocity].update(
    norm=plt.Normalize(-30, 30))
datatypes.TypePlotInfo[datatypes.SpectrumWidth].update(
    norm=plt.Normalize(0, 15))

# Helper for labelling the colorbar
def setup_cbar(cax, colorartist, units, pad=4):
    cbar = cax.colorbar(colorartist)
    if units and units != 'dimensionless':
        cbar.set_label_text(units)
    cbar.cbar_axis.labelpad = pad

# Helpers for multi-column plots
def multipanel_cbar_column(fig, layout, moments, data):

    grid = ImageGrid(fig, (1, 1, 1), nrows_ncols=layout, direction='row',
        share_all=True, axes_pad=0.45, aspect=True, cbar_mode='edge',
        cbar_location='bottom', cbar_pad=0.55, cbar_size='10%')

    for m, ax, cax, panel_label in zip(moments, grid, grid.cbar_axes,
            LabelGenerator('a')):
        ppi = PPIPlot(data.fields, var=m, ax=ax, rings=rings)
        panel_label.patch.set_boxstyle("round, pad=0., rounding_size=0.2")
        ax.add_artist(panel_label)
        ax.set_title(m)
        setup_cbar(cax, ppi.mesh, data.fields[m].dimensionality)

    return grid

def multipanel_cbar_row(fig, layout, moments, data):
    # TODO: The params here are not tuned
    grid = ImageGrid(fig, (1, 1, 1), nrows_ncols=layout, direction='row',
        share_all=True, axes_pad=0.45, aspect=True, cbar_mode='edge',
        cbar_location='right', cbar_pad=0.55, cbar_size='10%')

    for m, ax, cax, panel_label in zip(moments, grid, grid.cbar_axes,
            LabelGenerator('a')):
        ppi = PPIPlot(data.fields, var=m, ax=ax, rings=rings)
        panel_label.patch.set_boxstyle("round, pad=0., rounding_size=0.2")
        ax.add_artist(panel_label)
        ax.set_title(m)
        setup_cbar(cax, ppi.mesh, data.fields[m].dimensionality)

    return grid

def multipanel_cbar_each(fig, layout, moments, data):
    grid = ImageGrid(fig, (1, 1, 1), nrows_ncols=layout, direction='row',
        share_all=True, axes_pad=0.70, aspect=True, cbar_mode='each',
        cbar_location='right', cbar_pad=0.15, cbar_size='10%')

    for m, ax, cax, panel_label in zip(moments, grid, grid.cbar_axes,
            LabelGenerator('a')):
        ppi = PPIPlot(data.fields, var=m, ax=ax, rings=rings)
        panel_label.patch.set_boxstyle("round, pad=0., rounding_size=0.2")
        ax.add_artist(panel_label)
        ax.set_title(m)
        setup_cbar(cax, ppi.mesh, data.fields[m].dimensionality)

    return grid
