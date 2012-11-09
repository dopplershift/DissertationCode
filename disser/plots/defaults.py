import numpy as np
from matplotlib import rcParams
from functools import partial
from itertools import cycle
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
datatypes.TypePlotInfo[datatypes.Power].update(norm=plt.Normalize(-115, -25))
datatypes.TypePlotInfo[datatypes.ZDR].update(norm=plt.Normalize(-5, 5))
datatypes.TypePlotInfo[datatypes.PhiDP].update(norm=plt.Normalize(0, 75))
datatypes.TypePlotInfo[datatypes.KDP].update(norm=plt.Normalize(-5, 25))
datatypes.TypePlotInfo[datatypes.RhoHV].update(norm=plt.Normalize(0.98, 1.0))
datatypes.TypePlotInfo[datatypes.DiffAtten].update(norm=plt.Normalize(0, 5))
datatypes.TypePlotInfo[datatypes.SpecAttenuation].update(
        norm=plt.Normalize(0, 5))
datatypes.TypePlotInfo[datatypes.SpecDiffAtten].update(
        norm=plt.Normalize(0, 5))

# Set up some rcParams for figures
rcParams['savefig.dpi'] = 150
rcParams['font.size'] = 8
rcParams['figure.dpi'] = 107 # For laptop
rcParams['figure.figsize'] = (8, 4)
rcParams['figure.subplot.left'] = 0.05
rcParams['figure.subplot.right'] = 0.95
rcParams['figure.subplot.top'] = 0.98
rcParams['figure.subplot.bottom'] = 0.02

# Helper for labelling the colorbar
def setup_cbar(cax, colorartist, units, pad=4):
    cbar = cax.colorbar(colorartist)
    if units and units != 'dimensionless':
        cbar.set_label_text(units)
    cbar.cbar_axis.labelpad = pad

# Helper to turn a data object into a sequence if necessary
def make_data_iterator(data):
    if hasattr(data, 'fields'):
        return cycle([data])
    else:
        return data

# Helpers for multi-column plots
def multipanel_cbar_column(fig, layout, moments, data, rect=(1, 1, 1)):

    grid = ImageGrid(fig, rect, nrows_ncols=layout, direction='row',
        share_all=True, axes_pad=0.45, aspect=True, cbar_mode='edge',
        cbar_location='bottom', cbar_pad=0.55, cbar_size='10%')

    data = make_data_iterator(data)
    for d, m, ax, cax, panel_label in zip(data, moments, grid,
            grid.cbar_axes, LabelGenerator('a')):
        ppi = PPIPlot(d.fields, var=m, ax=ax, rings=rings)
        panel_label.patch.set_boxstyle("round, pad=0., rounding_size=0.2")
        ax.add_artist(panel_label)
        ax.set_title(m)
        setup_cbar(cax, ppi.mesh, d.fields[m].dimensionality)

    return grid

def multipanel_cbar_row(fig, layout, moments, data, rect=(1, 1, 1)):
    # TODO: The params here are not tuned
    grid = ImageGrid(fig, rect, nrows_ncols=layout, direction='column',
        share_all=True, axes_pad=0.15, aspect=True, cbar_mode='edge',
        cbar_location='right', cbar_pad=0.15, cbar_size='10%')

    data = make_data_iterator(data)
    for d, m, ax, cax, panel_label in zip(data, moments, grid,
            grid.cbar_axes, LabelGenerator('a')):
        ppi = PPIPlot(d.fields, var=m, ax=ax, rings=rings)
        panel_label.patch.set_boxstyle("round, pad=0., rounding_size=0.2")
        ax.add_artist(panel_label)
        setup_cbar(cax, ppi.mesh, '%s (%s)' % (m,
            d.fields[m].dimensionality))

    return grid

def multipanel_cbar_each(fig, layout, moments, data, rect=(1, 1, 1)):
    grid = ImageGrid(fig, rect, nrows_ncols=layout, direction='row',
        share_all=True, axes_pad=0.50, aspect=True, cbar_mode='each',
        cbar_location='right', cbar_pad=0.15, cbar_size='10%')

    data = make_data_iterator(data)
    for d, m, ax, cax, panel_label in zip(data, moments, grid,
            grid.cbar_axes, LabelGenerator('a')):
        ppi = PPIPlot(d.fields, var=m, ax=ax, rings=rings)
        panel_label.patch.set_boxstyle("round, pad=0., rounding_size=0.2")
        ax.add_artist(panel_label)
        ax.set_title(m)
        setup_cbar(cax, ppi.mesh, d.fields[m].dimensionality)

    return grid
