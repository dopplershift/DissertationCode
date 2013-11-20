# -*- coding: utf-8 -*-
import matplotlib.pyplot as plt
from matplotlib import rcParams
import numpy as np

from disser.io import DataCache
from disser.experiments import script_args, load_model_experiments
from disser import datatypes
import disser.plots.defaults as defaults

def spatial_sorter(k):
    # Map certain experiments to sort to front or end
    order = ['Control', 'Sidelobes', 'Beamwidth', 'Radial Width',
             'Range Resolution', 'Worst Case']
    return k[0], order.index(k[1])

def make_spatial_key(data):
    exp_key = 'Control'
    if data.runinfo.RadarAntennaCutoff == 'SideLobe1':
        exp_key = 'Sidelobes'
    if data.beamwidth > 1.0:
        exp_key = 'Beamwidth'
    if data.radialWidth > 1.9:
        exp_key = 'Radial Width'
    if data.gate_length > 125:
        exp_key = 'Range Resolution'
    if np.isnan(data.runinfo.FixedTemp):
        exp_key = 'Worst Case'
    return data.waveBand, exp_key

def load_spatial_experiments(data_dir, glob='*'):
    data_cache = DataCache(data_dir, make_spatial_key, ('band', 'exp'),
        pattern=glob)
    data_cache.key_sorter = spatial_sorter
    wavelengths,exps = data_cache.sub_keys()
    return data_cache

plt.rcParams['figure.figsize'] = (6, 6)
rcParams['figure.subplot.top'] = 0.95
rcParams['figure.subplot.bottom'] = 0.07

args = script_args().parse_args()

data_cache = load_spatial_experiments('spatial_runs', glob='Xband*')
data = data_cache[('X','Control')]

for band,exp in data_cache:
    data = data_cache[band,exp]
    fig_rect = [0.08, 0.08, 0.84, 0.84]
    title_text = '%s-band %s' % (data.waveBand, exp.title())
    basename = '%s_%s.png' % (band, exp)
    with datatypes.PlotInfoContext(wavelength=data.wavelength):
        moments = [data.fields.grab(moment, pol='H', source='ts')
                       for moment in (datatypes.Reflectivity,
                           datatypes.DopplerVelocity, datatypes.SpectrumWidth)]
        moments.append(data.fields.grab(datatypes.Attenuation, pol='H',
                source='calc'))
        grid = defaults.multipanel_cbar_each(plt.figure(), (2, 2), moments,
                data, rect=fig_rect)
        title = grid[0].figure.suptitle(title_text, fontsize=11)
    if args.save:
        grid[0].figure.savefig('single_' + basename)

    with datatypes.PlotInfoContext(wavelength=data.wavelength):
        moments = [data.fields.grab(moment, pol='H', source='ts')
                       for moment in (datatypes.ZDR, datatypes.RhoHV,
                           datatypes.PhiDP)]
        moments.append(data.fields.grab(datatypes.DiffAtten, pol='diff',
                source='calc'))
        grid = defaults.multipanel_cbar_each(plt.figure(), (2, 2), moments,
                data, rect=fig_rect)
        title = grid[0].figure.suptitle(title_text, fontsize=11)
    if args.save:
        grid[0].figure.savefig('dual_' + basename)

if not args.save:
    plt.show()
