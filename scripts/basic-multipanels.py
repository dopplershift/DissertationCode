import matplotlib.pyplot as plt
from disser.io import NetCDFRadarData
from disser import datatypes
import disser.plots.defaults as defaults

def make_plots(fname, save_prefix=None):
    data = NetCDFRadarData(fname)

    moments_fig1 = [data.fields.grab(moment, pol='H', source='ts')
               for moment in (datatypes.Power, datatypes.Reflectivity,
                   datatypes.DopplerVelocity, datatypes.SpectrumWidth,
                   datatypes.ZDR, datatypes.RhoHV, datatypes.KDP,
                   datatypes.PhiDP)]

    moments_fig2 = [data.fields.grab(moment, pol='H', source='average')
               for moment in (datatypes.Power, datatypes.Reflectivity,
                   datatypes.DopplerVelocity, datatypes.SpectrumWidth,
                   datatypes.ZDR, datatypes.BackscatterPhase, datatypes.KDP,
                   datatypes.PhiDP)]

    figs = []
    for moms in [moments_fig1, moments_fig2]:
        figs.append(plt.figure(figsize=(8, 4), dpi=107))

        # Need some way to control the formatting of a datatype into a string
        # Maybe using a "global" setting and a context maanger
        grid = defaults.multipanel_cbar_each(figs[-1], (2,4), moms, data)
        ax = grid[0]

        ax.xaxis.set_major_locator(plt.MultipleLocator(10))
        ax.set_ylim(0, 50)
        ax.set_xlim(-20, 20)

    if save_prefix:
        for num,f in enumerate(figs):
            fig.savefig(save_prefix + '-%d.png' % num)
    else:
        plt.show()

if __name__ == '__main__':
    import sys
    import os.path
    if len(sys.argv) < 2:
        path = 'Sband_3600_20111010_210552.nc'
    else:
        path = sys.argv[1]

    if os.path.isdir(path):
        import glob
        for fname in glob.glob(os.path.join(path,'*.nc')):
            print 'Processing %s...' % os.path.split(fname)[-1],
            sys.stdout.flush()
            img_prefix = os.path.splitext(fname)[0]
            make_plots(fname, img_prefix)
            print 'done.'
    else:
        make_plots(path)
