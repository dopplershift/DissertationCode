'''
This is a module of things useful for making quick analysis script, but
are specific to my datasets.
'''

import numpy as np
import quantities as pq
import scipy.stats as ss

from .io import DataCache
from .atten import calc_specific_atten, attenAlgs
from . import datatypes

def regress_stats(truth, data):
    residuals = data - truth
    bias = residuals.mean()
    mse = (residuals * residuals).mean()
    corr = ss.pearsonr(truth, data)[0]
    return bias, mse, corr*corr

def script_args(desc=''):
    from argparse import ArgumentParser
    parser = ArgumentParser(description=desc)
    parser.add_argument('-s', '--save', action='store_true', help='Save plots',
            default=False)
    return parser

def exp_sorter(k):
    # Map certain experiments to sort to front or end
    return k[0], {'Control':'@', 'Combined':'}'}.get(k[1], k[1])

def make_spatial_key(data):
    exp_key = 'Control'
    if data.runinfo.AxisRatioCalc != 'Brandes':
        exp_key = 'Combined'
    elif data.runinfo.RadarAntennaCutoff != 'MainLobe':
        exp_key = 'Sidelobe'
    elif data.gate_length > 125:
        exp_key = 'Range Resolution'
    elif data.radialWidth > 1.5:
        exp_key = 'Radial Width'
    elif data.beamwidth > 1.5:
        exp_key = 'Beamwidth'
    return data.waveBand, exp_key

def make_model_key(data):
    diff_count = 0
    exp_key = 'Control'
    if np.abs(data.wavelength.rescale(pq.cm) -
            np.round(data.wavelength.rescale(pq.cm), 0)) < 0.1:
        exp_key = 'Wavelength'
        diff_count += 1
    if np.isnan(data.runinfo.FixedTemp):
        exp_key = 'Temperature'
        diff_count += 1
    if data.runinfo.CantingWidth > 10.0:
        exp_key = 'Canting'
        diff_count += 1
    if data.runinfo.AxisRatioCalc != 'Brandes':
        exp_key = 'Shape'
        diff_count += 1
    if diff_count > 1:
        exp_key = 'Combined'
    return data.waveBand, exp_key

def load_spatial_experiments(data_dir, glob='*'):
    data_cache = DataCache(data_dir, make_spatial_key, ('band', 'exp'),
        pattern=glob)
    data_cache.key_sorter = exp_sorter
    wavelengths,exps = data_cache.sub_keys()
    return data_cache

def load_model_experiments(data_dir, glob='*'):
    data_cache = DataCache(data_dir, make_model_key, ('band', 'exp'),
        pattern=glob)
    data_cache.key_sorter = exp_sorter
    wavelengths,exps = data_cache.sub_keys()
    return data_cache

def process_all_atten(data_cache, default_source='average'):
    # This way we can operate on an whole cache or a single data item
    if hasattr(data_cache, 'fields'):
        data = [data_cache]
    else:
        data = data_cache.values()

    for d in data:
        d.fields.default_keys['source'] = default_source
        attenAlgs.runAll(d, var='H')
        attenAlgs.runAll(d, var='V')
        attenAlgs.runAll(d, var='diff')
        calc_specific_atten(d, datatypes.Attenuation, pol='H')
        calc_specific_atten(d, datatypes.Attenuation, pol='V')
        calc_specific_atten(d, datatypes.DiffAtten, pol=None)
