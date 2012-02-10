from collections import namedtuple
from .plugintools import PluginRegistry

# Used to contain default information for plotting data types.
TypePlotInfo = dict()

MomentRegistry = PluginRegistry()

# Units could possibly go if we use Quantities
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


Reflectivity = DataType(name='Reflectivity', abbr='Z')
DopplerVelocity = DataType(name='Doppler Velocity', abbr=r'$V_r$')
SpectrumWidth = DataType(name='Spectrum Width', abbr=r'$\sigma_v$')

ZDR = DataType(name='Differential Reflectivity', abbr=r'$Z_{DR}$')
KDP = DataType(name='Specific Differential Phase', abbr=r'$K_{DP}$')
PhiDP = DataType(name='Differential Propagation Phase', abbr=r'$\Phi$')

Attenuation = DataType(name='Attenuation', abbr='A')
DiffAtten = DataType(name='Differential Attenuation', abbr=r'$A_D$')
