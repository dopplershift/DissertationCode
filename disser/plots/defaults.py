import numpy as np
from matplotlib import rcParams
from functools import partial, wraps
from itertools import cycle
from contextlib import contextmanager
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

# Allow for selecting a norm based on wavelength, which is set
# using a context manager elsewhere.
class wavelengthNorm(object):
    attribute_wrapper = True
    def __init__(self, X, C, S):
        self._xnorm = X
        self._cnorm = C
        self._snorm = S
        self._attrMap = datatypes.plotInfoAttrs

    def __call__(self):
        wavelength = self._attrMap['wavelength']
        if wavelength is None or wavelength < 0.04:
            return self._xnorm
        elif wavelength >= 0.08:
            return self._snorm
        else:
            return self._cnorm

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

_phiNorms = wavelengthNorm(X=plt.Normalize(0, 250),
        C=plt.Normalize(0, 200), S=plt.Normalize(0, 100))
datatypes.TypePlotInfo[datatypes.PhiDP].update(norm=_phiNorms)

_kdpNorms = wavelengthNorm(X=plt.Normalize(-5, 25), C=plt.Normalize(-5, 15),
        S=plt.Normalize(-5, 5))
datatypes.TypePlotInfo[datatypes.KDP].update(norm=_kdpNorms)
datatypes.TypePlotInfo[datatypes.RhoHV].update(norm=plt.Normalize(0.98, 1.0))

_attenNorms = wavelengthNorm(X=plt.Normalize(0, 100), C=plt.Normalize(0, 20),
        S=plt.Normalize(0, 5))
datatypes.TypePlotInfo[datatypes.Attenuation].update(norm=_attenNorms)

_diffANorms = wavelengthNorm(X=plt.Normalize(0, 20), C=plt.Normalize(0, 5),
        S=plt.Normalize(0, 2))
datatypes.TypePlotInfo[datatypes.DiffAtten].update(norm=_diffANorms)

_specAttenNorms = wavelengthNorm(X=plt.Normalize(0, 8), C=plt.Normalize(0, 2),
        S=plt.Normalize(0, 0.25))
datatypes.TypePlotInfo[datatypes.SpecAttenuation].update(norm=_specAttenNorms)

_specDANorms = wavelengthNorm(X=plt.Normalize(0, 1.0),
        C=plt.Normalize(0, 0.75), S=plt.Normalize(0, 0.1))
datatypes.TypePlotInfo[datatypes.SpecDiffAtten].update(norm=_specDANorms)

bidi_cmap = get_cmap('Carbone42')
datatypes.TypePlotInfo[datatypes.AttenDelta].update(
        norm=plt.Normalize(-10, 10), cmap=bidi_cmap)
datatypes.TypePlotInfo[datatypes.SpecAttenDelta].update(
        norm=plt.Normalize(-1, 1), cmap=bidi_cmap)
datatypes.TypePlotInfo[datatypes.DiffAttenDelta].update(
        norm=plt.Normalize(-2, 2), cmap=bidi_cmap)
datatypes.TypePlotInfo[datatypes.SpecDiffAttenDelta].update(
        norm=plt.Normalize(-0.5, 0.5), cmap=bidi_cmap)
datatypes.TypePlotInfo[datatypes.PhiDelta].update(norm=plt.Normalize(-50, 50),
        cmap=bidi_cmap)

# Set up some rcParams for figures
rcParams['savefig.dpi'] = 150
rcParams['font.size'] = 8
rcParams['figure.dpi'] = 107 # For laptop
rcParams['figure.figsize'] = (8, 4)
rcParams['figure.subplot.left'] = 0.05
rcParams['figure.subplot.right'] = 0.95
rcParams['figure.subplot.top'] = 0.98
rcParams['figure.subplot.bottom'] = 0.02

class AxisDefaults(object):
    def __init__(self):
        self.setup = lambda ax: None

    def __call__(self, func):
        @wraps(func)
        def wrapper(*args, **kwargs):
            grid = func(*args, **kwargs)
            self.setup(grid)
            return grid
        return wrapper

axisDefaults = AxisDefaults()

def default_ppi_axis(grid):
    ax = grid[0]
    ax.xaxis.set_major_locator(plt.MultipleLocator(10))
    ax.set_ylim(0, 50)
    ax.set_xlim(-20, 20)

axisDefaults.setup = default_ppi_axis

def defaultColorbarLabels(dt, units):
    return units if units != 'dimensionless' else ''

def sourceLabels(dt, units):
    return '%s (%s)' % (dt, units)

def algLabels(dt, units):
    abbr, src_str = dt.string_parts()
    return '%s from %s (%s)' % (abbr, src_str, units)

class ColorbarLabeller(object):
    def __init__(self):
        self.label = defaultColorbarLabels

    @contextmanager
    def __call__(self, new_label):
        try:
            # This way we ensure we aren't being nested
            if self.label == defaultColorbarLabels:
                saved = self.label
                self.label = new_label
            else:
                saved = None
            yield
        finally:
            if saved:
                self.label = saved
colorbarLabeller = ColorbarLabeller()

# Helper for labelling the colorbar
def setup_cbar(cax, colorartist, datatype, units, pad=4):
    cbar = cax.colorbar(colorartist)
    label = colorbarLabeller.label(datatype, units)
    if label:
        cbar.set_label_text(label)
    cbar.cbar_axis.labelpad = pad

# Helper to turn a data object into a sequence if necessary
def make_data_iterator(data):
    if hasattr(data, 'fields'):
        return cycle([data])
    else:
        return data

# Helpers for multi-column plots
@axisDefaults
def multipanel_cbar_column(fig, layout, moments, data, rect=(1, 1, 1)):

    grid = ImageGrid(fig, rect, nrows_ncols=layout, direction='row',
        share_all=True, axes_pad=0.45, aspect=True, cbar_mode='edge',
        cbar_location='bottom', cbar_pad=0.55, cbar_size='10%')

    use_labels = layout != (1,1)
    data = make_data_iterator(data)
    for d, m, ax, cax, panel_label in zip(data, moments, grid,
            grid.cbar_axes, LabelGenerator('a')):
        ppi = PPIPlot(d.fields, var=m, ax=ax, rings=rings)
        if use_labels:
            panel_label.patch.set_boxstyle("round, pad=0., rounding_size=0.2")
            ax.add_artist(panel_label)
        ax.set_title(m)
        setup_cbar(cax, ppi.mesh, m, d.fields[m].dimensionality)

    return grid

@axisDefaults
def multipanel_cbar_row(fig, layout, moments, data, rect=(1, 1, 1)):
    # TODO: The params here are not tuned
    grid = ImageGrid(fig, rect, nrows_ncols=layout, direction='column',
        share_all=True, axes_pad=0.15, aspect=True, cbar_mode='edge',
        cbar_location='right', cbar_pad=0.15, cbar_size='10%')

    use_labels = layout != (1,1)
    data = make_data_iterator(data)

    for d, m, ax, cax, panel_label in zip(data, moments, grid,
            grid.cbar_axes, LabelGenerator('a')):
        ppi = PPIPlot(d.fields, var=m, ax=ax, rings=rings)
        if use_labels:
            panel_label.patch.set_boxstyle("round, pad=0., rounding_size=0.2")
            ax.add_artist(panel_label)
        setup_cbar(cax, ppi.mesh, m, d.fields[m].dimensionality)

    return grid

@axisDefaults
def multipanel_cbar_each(fig, layout, moments, data, rect=(1, 1, 1)):
    grid = ImageGrid(fig, rect, nrows_ncols=layout, direction='row',
        share_all=True, axes_pad=0.65, aspect=True, cbar_mode='each',
        cbar_location='right', cbar_pad=0.15, cbar_size='10%')

    use_labels = layout != (1,1)
    data = make_data_iterator(data)
    for d, m, ax, cax, panel_label in zip(data, moments, grid,
            grid.cbar_axes, LabelGenerator('a')):
        ppi = PPIPlot(d.fields, var=m, ax=ax, rings=rings)
        if use_labels:
            panel_label.patch.set_boxstyle("round, pad=0., rounding_size=0.2")
            ax.add_artist(panel_label)
        ax.set_title(m)
        setup_cbar(cax, ppi.mesh, m, d.fields[m].dimensionality)

    return grid
