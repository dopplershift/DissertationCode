from collections import namedtuple, defaultdict
from .plugintools import PluginRegistry
# For the plot information, we use a default dictionary that returns None
# for any unknown key, which results in using the matplotlib default.
# We can also add hooks for some attributes, wherein if the value stored
# for some plot element is a callable, it is called
class TypePlotInfoSource(defaultdict):
    def __init__(self, *args, **kwargs):
        defaultdict.__init__(self, *args, **kwargs)
        self.default_factory = lambda: None

    def __getitem__(self, key):
        item = defaultdict.__getitem__(self, key)
        if callable(item) and hasattr(item, 'attribute_wrapper'):
            return item()
        else:
            return item

plotInfoAttrs = defaultdict(lambda: None)

class PlotInfoContext(object):
    def __init__(self, **kwargs):
        self._mods = kwargs

    def __enter__(self):
        self._old_settings = plotInfoAttrs.copy()
        plotInfoAttrs.update(self._mods)
        return plotInfoAttrs

    def __exit__(self, type, value, tb):
        plotInfoAttrs.clear()
        plotInfoAttrs.update(self._old_settings)

# Used to contain default information for plotting data types.
# What we want is to able to look up information for any data type,
# and if it's not present, return some default plot information.
class TypePlotInfoLookup(defaultdict):
    def __init__(self, **kwargs):
        # Makes us return a new copy of the defaults for an unknown key
        defaultdict.__init__(self, self.make_default, **kwargs)
        self.defaults = TypePlotInfoSource()

    # Used to update the defaults
    def set_defaults(self, **kwargs):
        self.defaults.update(kwargs)

    # Creates and returns a copy of the default information
    def make_default(self):
        return self.defaults.copy()

TypePlotInfo = TypePlotInfoLookup()

MomentRegistry = PluginRegistry()

# Might want to eliminate auto-registry magic
class DataType(namedtuple('DataType', ['name', 'abbr'])):
    def __init__(self, *args, **kwargs):
        super(DataType, self).__init__(*args, **kwargs)
        MomentRegistry.register(self)
        self.sources = list()

    def __del__(self):
        if not MomentRegistry is None:
            MomentRegistry.unregister(self)

    def add_source(self, algorithm):
        self.sources.append(algorithm)


def getBestAlg(datatype, availTypes, unavailTypes=None):
    if unavailTypes == None:
        unavailTypes = list()
    else:
        unavailTypes = list(unavailTypes)

    for alg in datatype.sources:
        for src in alg.sources:
            if src not in availTypes:
                if getBestAlg(src, availTypes, unavailTypes + [datatype]) is None:
                    return None
        else:
            return alg
    else:
        return None


Power = DataType(name='Power', abbr='$P_r$')
SNR = DataType(name='Signal to Noise Ratio', abbr='SNR')
Reflectivity = DataType(name='Reflectivity', abbr='Z')
DopplerVelocity = DataType(name='Doppler Velocity', abbr=r'$V_r$')
SpectrumWidth = DataType(name='Spectrum Width', abbr=r'$\sigma_v$')

ZDR = DataType(name='Differential Reflectivity', abbr=r'$Z_{DR}$')
KDP = DataType(name='Specific Differential Phase', abbr=r'$K_{DP}$')
PhiDP = DataType(name='Differential Propagation Phase', abbr=r'$\Phi$')
RhoHV = DataType(name='Co-Polar Cross-Correlation Coefficient',
    abbr=r'$\rho_{HV}$')
BackscatterPhase = DataType(name='Differential Phase on Backscatter',
    abbr=r'$\delta$')

Attenuation = DataType(name='Attenuation', abbr='A')
DiffAtten = DataType(name='Differential Attenuation', abbr=r'$A_D$')

SpecAttenuation = DataType(name='Specific Attenuation', abbr=r'$\alpha$')
SpecDiffAtten = DataType(name='Specific Differential Attenuation',
        abbr=r'$\alpha_D$')
