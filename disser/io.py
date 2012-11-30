from collections import namedtuple
import numpy as np
import netCDF4
import quantities as pq
from quantities import sin, cos
from . import units, datatypes, calc
from .sigproc import auto_moments, auto_dual_pol, shift_phi


# This could be used to keep some information attached to an array. Maybe.
class MetadataArray(object):

    def __init__(self, data):
        self.data = data

    def __setattr__(self, attr, val):
        if hasattr(self.data, attr):
            setattr(self.data, attr, val)
        else:
            self.__dict__[attr] = val

    def __getattr__(self, attr):
        if attr in self.__dict__:
            return self.__dict__[attr]
        else:
            return getattr(self.data, attr)


class AttributeDict(dict):

    __getattr__ = dict.__getitem__
    __setattr__ = dict.__setitem__


# This is a class to be able to query for the best data match to a given
# set of metadata. This includes the data type, as well potentially other
# matches like the height of the ob (i.e. mesonet data with air temp at 1m and
# 10m) selecting radar data from time series vs. averaged and giving the
# polarization
class FieldStore(dict):
    def __init__(self, **kwargs):
        super(FieldStore, self).__init__(**kwargs)
        self.default_keys = {}

    def grabAll(self, datatype, filt=None, **keys):
        if filt is None:
            filt = lambda k: True
        keys.update(self.default_keys)
        potential = [k for k in self.keys() if k[0] is datatype]
        return filter(filt, sorted(potential, key=self.sorter(**keys)))

    def grab(self, datatype, filt=None, **keys):
        return self.grabAll(datatype, filt, **keys)[-1]

    def grabData(self, datatype, filt=None, **keys):
        return self[self.grab(datatype, **keys)]

    def sorter(self, **keys):
        return lambda key: sum(key[key._fields.index(field)] == v
            for ind, (field, v) in enumerate(keys.items())
            if field in key._fields)


class DataSet(AttributeDict):

    def __init__(self):
        self.fields = FieldStore()
        self._calc_cache = dict()
        self.metadata = []
        self.coordinates = []

    def addField(self, data, *args, **kwargs):
        key = self.fieldKey(*args, **kwargs)
        self.fields[key] = data

    def fieldKey(self, *args, **kwargs):
        if args:
            return args[0]
        elif kwargs:
            return kwargs[kwargs.keys()[0]]
        else:
            return None

# Is there some way that we could almost have a database here of the fields
# contained within our dataset. Each type matches a global type as well as
# having extra metadata that we could match if possible?

#    def getField(datatype, prefer=None):
#        if datatype in self.fields:
#            return self.fields[datatype]
#        elif (datatype,prefer) in self._calc_cache:
#            return self._calc_cache[datatype, prefer]
#        else:
#            alg = getBestAlg(datatype, self.fields.keys())
#            return alg(self)

class NetCDFData(DataSet):
    def __init__(self, fname):
        DataSet.__init__(self)
        nc = netCDF4.Dataset(fname, 'r')
        self.nc = nc

    def readVar(self, varname):
        var = self.nc.variables[varname]
        vals = var[:]

        # Check if we were given a masked array. For now, we'll fill masks
        # with NaNs
        if(hasattr(vals, 'fill_value')):
            vals.fill_value = np.nan
            vals = vals.filled()

        try:
            return pq.Quantity(vals, var.units.replace(' ', '_'))
        except LookupError:
            print 'Need to add support for: %s' % var.units
            return vals


class ModelData(NetCDFData):
    def __init__(self, fname):
        super(ModelData, self).__init__(fname)
        self.pressure = self.readVar('p')
        theta = self.readVar('theta')
        self.temp = calc.theta_to_temperature(theta, self.pressure)
        self.qr = calc.mixing_ratio_to_density(self.readVar('qr'),
                calc.air_density(self.temp, self.pressure))
        self.nr = self.readVar('nr')


# Make a base for moment info from a named tuple
class MomentInfo_(namedtuple('MomentInfo', ['type', 'pol', 'source'])):
    def __str__(self):
        abbr = self.type.abbr
        if self.pol:
            abbr = abbr.replace('$', '')
            if '_' in abbr:
                p1,p2 = abbr.split('_')
                abbr = (p1 + '_' + '{' + p2.replace('}', '').replace('{', '')
                    + self.pol + '}')
            else:
                abbr = abbr + '_{%s}' % self.pol
            name = '$' + abbr + '$'
        else:
            name = abbr

        if self.source:
            return name + ' (%s)' % self.source.capitalize()
        else:
            return name

MIDefault = MomentInfo_(*([None] * len(MomentInfo_._fields)))

# Extend from this to allow None as a default for some parameters. This is
# easier than subclassing to implement default args.
def MomentInfo(datatype, **kwargs):
    return MIDefault._replace(type=datatype, **kwargs)


class NetCDFRadarData(NetCDFData):
    def __init__(self, fname, freeIQ=True):
        NetCDFData.__init__(self, fname)
        self._band = None

        self.radar = self.nc.RadarName

        runinfo = AttributeDict()
        self.runinfo = runinfo
        for attr in ['CantingWidth', 'ClearAir', 'ConfigFilename',
            'DataFiles', 'DataType', 'DropModel', 'FixedTemp',
            'GitDate', 'GitHash', 'GitTimeStamp', 'RadarName',
            'RandomSeed', 'RunStarted', 'ScatteringModel', 'SweepType',
            'VersionNumber', 'AxisRatioCalc']:
            runinfo[attr] = getattr(self.nc, attr, None)

        self.wavelength = self.readVar('Wavelength')[0]
        self.az = self.readVar('Azimuth')
        self.gate_length = self.readVar('GateLength')
        self.pulse_length = self.readVar('PulseDuration')[0] * pq.c / 2.
        self.rng = (np.arange(len(self.nc.dimensions['gates']))
            * self.gate_length + self.readVar('RangeToFirstGate'))

        self.nyquist = self.readVar('NyquistVelocity')
        self.noise_pwr_db = self.readVar('MinimumDetectableSignal')
        self.noise_pwr = units.to_linear(self.noise_pwr_db)

        self.xlocs = (self.rng * sin(self.az[:, np.newaxis])).rescale('km')
        self.ylocs = (self.rng * cos(self.az[:, np.newaxis])).rescale('km')

        self.vel = self.readVar('Velocity')
        self.fields[MomentInfo(datatypes.DopplerVelocity, source='average',
            pol='H')] = self.vel
        self.spw = self.readVar('SpectrumWidth')
        self.fields[MomentInfo(datatypes.SpectrumWidth, source='average',
            pol='H')] = self.spw

        self.process_channel('H')
        self.process_channel('V')
        self.zdr_ts, self.rhohv, self.phidp_ts = auto_dual_pol(self.iq_H,
            self.iq_V, self.noise_pwr, self.noise_pwr)
        self.rhohv = self.rhohv * pq.dimensionless
        self.fields[MomentInfo(datatypes.RhoHV, source='ts')] = self.rhohv

        self.zdr_ts = units.dB * self.zdr_ts
        self.fields[MomentInfo(datatypes.ZDR, source='ts')] = self.zdr_ts

        self.phidp_ts = pq.radians * self.phidp_ts
        self.fields[MomentInfo(datatypes.PhiDP, source='ts')] = self.phidp_ts
        self.phidp_ts.units = pq.degrees
        shift_phi(self.phidp_ts)

        self.zdr = units.make_dB(self.ref_H - self.ref_V)
        self.fields[MomentInfo(datatypes.ZDR, source='average')] = self.zdr

        self.kdp = self.readVar('KDP')
        self.kdp.units = pq.degrees / pq.kilometer
        self.fields[MomentInfo(datatypes.KDP, source='average')] = self.kdp

        self.phidp = self.readVar('PhiDP')
        self.phidp.units = pq.degrees
        shift_phi(self.phidp)
#        self.phidp = remainder(self.phi_dp, 360.)
#        self.phidp2 = 2*self.kdp.cumsum(axis=1) * (self.pulse_length / kilo)
        self.fields[MomentInfo(datatypes.PhiDP, source='average')] = self.phidp

        thgrad, rgrad = np.gradient(self.phidp_ts, 1,
                self.gate_length.rescale(pq.kilometer))
        self.kdp_ts = rgrad / 2.
        self.fields[MomentInfo(datatypes.KDP, source='ts')] = self.kdp_ts

        self.diff_atten = (units.make_dB(self.unatten_pwr_H -
            self.unatten_pwr_V) - self.zdr)
        self.fields[MomentInfo(datatypes.DiffAtten,
            source='average')] = self.diff_atten

        self.diff_atten_ts = (units.make_dB(self.unatten_pwr_H -
            self.unatten_pwr_V) - self.zdr_ts)
        self.fields[MomentInfo(datatypes.DiffAtten,
            source='ts')] = self.diff_atten_ts

        mean_diff_atten = self.nc.variables['MeanDiffAtten'][:]
        # This way, when we take the mean, any set of pulses with a 0
        # ends up nan. This way empty attenuation pulses don't end up
        # producing a biased value, they just end up masked
        mean_diff_atten[mean_diff_atten == 0.0] = np.nan
        self.mean_diff_atten = -units.to_dB(mean_diff_atten.mean(axis=-1))
        self.fields[MomentInfo(datatypes.DiffAtten,
            source='calc')] = self.mean_diff_atten

        self.spec_mean_diff_atten = np.gradient(self.mean_diff_atten, 1,
                self.gate_length)[1].rescale('dB/km')
        self.fields[MomentInfo(datatypes.SpecDiffAtten,
            source='calc')] = self.spec_mean_diff_atten

        self.spec_diff_atten_ts = np.gradient(self.diff_atten_ts, 1,
                self.gate_length)[1].rescale('dB/km')
        self.fields[MomentInfo(datatypes.SpecDiffAtten,
            source='ts')] = self.spec_diff_atten_ts

        self.spec_diff_atten = np.gradient(self.diff_atten, 1,
                self.gate_length)[1].rescale('dB/km')
        self.fields[MomentInfo(datatypes.SpecDiffAtten,
            source='average')] = self.spec_diff_atten

        self.delta = self.readVar('Delta')
        self.delta.units = pq.degrees
        self.fields[MomentInfo(datatypes.BackscatterPhase,
            source='average')] = self.delta

        # TODO: Need to read in the diagnostic variables.
        self.fields['x'] = self.xlocs
        self.fields['y'] = self.ylocs

        # If we're allowed, delete the IQ data to reduce memory usage.
        if freeIQ:
            del self['iq_H']
            del self['iq_V']

    def process_channel(self, pol):
        # Read moments directly from file
        self['ref_' + pol] = self.readVar('Reflectivity_%s' % pol)
        self.fields[MomentInfo(datatypes.Reflectivity, pol=pol,
            source='average')] = self['ref_' + pol]

        self['pwr_' + pol] = units.dBW_to_dBm(self.readVar('Power_%s' % pol))
        self.fields[MomentInfo(datatypes.Power, pol=pol,
            source='average')] = self['pwr_' + pol]

        self['unatten_pwr_' + pol] = units.dBW_to_dBm(
            self.readVar('UnattenPower_%s' % pol))

        # Read time series data and calculate moments
        self['iq_' + pol] = (self.nc.variables['I_%s' % pol][:]
            + 1.0j * self.nc.variables['Q_%s' % pol][:])
        pwr_ts, vel_ts, spw_ts = auto_moments(self['iq_' + pol], self.nyquist)
        pwr_ts = pq.watt * pwr_ts
        pwr_ts_dbm = units.to_dBm(pwr_ts)
        self['pwr_ts_' + pol] = pwr_ts
        self['pwr_ts_%s_dbm' % pol] = pwr_ts_dbm
        self.fields[MomentInfo(datatypes.Power, pol=pol,
            source='ts')] = self['pwr_ts_%s_dbm' % pol]

        self['vel_ts_' + pol] = vel_ts
        self.fields[MomentInfo(datatypes.DopplerVelocity, pol=pol,
            source='ts')] = self['vel_ts_' + pol]

        self['spw_ts_' + pol] = spw_ts
        self.fields[MomentInfo(datatypes.SpectrumWidth, pol=pol,
            source='ts')] = self['spw_ts_' + pol]

        # Convert time series power estimate to reflectivity factor
        self['rad_const_%s_dB' % pol] = units.make_dB(
            self.readVar('RefCalibration_%s' % pol))
        self['rad_const_' + pol] = units.to_linear(
            self['rad_const_%s_dB' % pol]) * pq.mm**6 / (pq.m**5 * pq.watt)

        self['ref_ts_' + pol] = units.to_dBz(
            pwr_ts * self.rng**2 * self['rad_const_' + pol])
        self.fields[MomentInfo(datatypes.Reflectivity, pol=pol,
            source='ts')] = self['ref_ts_' + pol]

        # Calculate attenuation from fields
        self['atten_' + pol] = units.make_dB(
            self['unatten_pwr_' + pol] - self['pwr_' + pol])
        self.fields[MomentInfo(datatypes.Attenuation, pol=pol,
            source='average')] = self['atten_' + pol]

        self['atten_ts_' + pol] = units.make_dB(
            self['unatten_pwr_' + pol] - self['pwr_ts_%s_dbm' % pol])
        self.fields[MomentInfo(datatypes.Attenuation, pol=pol,
            source='ts')] = self['atten_ts_' + pol]

        mean_atten = self.nc.variables['MeanAtten'][:]
        # This way, when we take the mean, any set of pulses with a 0
        # ends up nan. This way empty attenuation pulses don't end up
        # producing a biased value, they just end up masked
        mean_atten[mean_atten == 0.0] = np.nan
        self['mean_atten_' + pol] = -units.to_dB(mean_atten.mean(axis=-1))
        self.fields[MomentInfo(datatypes.Attenuation, pol=pol,
            source='calc')] = self['mean_atten_' + pol]

        self['spec_mean_atten_' + pol] = np.gradient(self['mean_atten_' + pol],
                1, self.gate_length)[1].rescale('dB/km')
        self.fields[MomentInfo(datatypes.SpecAttenuation, pol=pol,
            source='calc')] = self['spec_mean_atten_' + pol]

        self['spec_atten_ts_' + pol] = np.gradient(self['atten_ts_' + pol], 1,
                self.gate_length)[1].rescale('dB/km')
        self.fields[MomentInfo(datatypes.SpecAttenuation, pol=pol,
            source='ts')] = self['spec_atten_ts_' + pol]

        self['spec_atten_' + pol] = np.gradient(self['atten_' + pol], 1,
                self.gate_length)[1].rescale('dB/km')
        self.fields[MomentInfo(datatypes.SpecAttenuation, pol=pol,
            source='average')] = self['spec_atten_' + pol]

    @property
    def waveBand(self):
        if self._band is None:
            if self.wavelength > 8 * pq.cm:
                band = 'S'
            elif self.wavelength > 4 * pq.cm:
                band = 'C'
            else:
                band = 'X'
            self._band = band
        return self._band

    def fieldKey(self, *args, **kwargs):
        return MomentInfo(*args, **kwargs)


class DataCache(dict):
    '''Class to simplify mass loading of data files into a cache.'''
    def __init__(self, dirpath, keygen, klass=NetCDFRadarData, pattern='*'):
        self._path = dirpath
        self._keygen = keygen
        self._dataClass = klass
        self._pattern = pattern
        self.load()
        self.key_sorter = lambda k: k

    def load(self):
        import glob
        import os.path
        self.clear()
        for datafile in glob.glob(os.path.join(self._path, self._pattern)):
            data = self._dataClass(datafile)
            key = self._keygen(data)
            self[key] = data

    def sub_keys(self):
        'Returns a list of the sets of unique sub keys in the cache.'
        # Initialize a list of items for each position
        keylists = []
        for i in range(len(self.keys()[0])):
            keylists.append([])

        # Sort all the keys and loop over, adding each tuple position's value
        # only once.
        for key in sorted(self.keys(), key=self.key_sorter):
            for ind, item in enumerate(key):
                if item not in keylists[ind]:
                    keylists[ind].append(item)

        return keylists
