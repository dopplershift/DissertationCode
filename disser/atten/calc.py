import numpy as np
from .. import datatypes
def calc_specific_atten(data, dt=datatypes.Attenuation, pol='H'):
    destMap = {datatypes.Attenuation:datatypes.SpecAttenuation,
            datatypes.DiffAtten:datatypes.SpecDiffAtten}
    fields = data.fields.grabAll(dt, filt=lambda f: f.pol==pol)
    for f in fields:
        d = data.fields[f]
        spec = np.gradient(d, 1, data.gate_length)[1].rescale('dB/km')
        data.addField(spec, destMap[dt], pol=f.pol, source=f.source)
