# Generate PPI plots of Z and Attenuation for (calc, zphi, sc)
import matplotlib.pyplot as plt
from matplotlib import rcParams

from disser.experiments import (script_args, load_model_experiments,
        process_all_atten)
from disser.atten import attenAlgs
from disser.datatypes import Reflectivity, Attenuation, PlotInfoContext
from disser.plots.defaults import multipanel_cbar_each

rcParams['figure.subplot.top'] = 0.95
rcParams['figure.subplot.bottom'] = 0.07

args = script_args().parse_args()

data_cache = load_model_experiments('ref_runs', glob='Cband*')
data = data_cache[('C','Control')]
data.fields.default_keys['source'] = 'ts'
attenAlgs.run(['ZPHI', 'SC'], data, var='H')

data.fields.default_keys['pol'] = 'H'
with PlotInfoContext(wavelength=data.wavelength):
    for dt,src in [(Reflectivity, 'ts'), (Attenuation, 'calc'),
            (Attenuation, 'ZPHI'), (Attenuation, 'SC')]:
        fig = plt.figure(figsize=(8, 6))
        mom = data.fields.grab(dt, source=src)
        multipanel_cbar_each(fig, (1, 1), [mom], data)
        if args.save:
            fig.savefig('{.name}_{:s}.png'.format(dt, src))

if not args.save:
    plt.show()
