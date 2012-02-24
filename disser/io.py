from collections import namedtuple
import numpy as np
import netCDF4
import quantities as pq
from quantities import sin, cos
from . import units, datatypes
from .sigproc import auto_moments, auto_dual_pol


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
    def grab(self, datatype, **keys):
        potential = [k for k in self.keys() if k[0] is datatype]
        return sorted(potential, key=self.sorter(**keys))[-1]
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

# Make a base for moment info from a named tuple
MIBase = namedtuple('MomentInfo', ['type', 'pol', 'source'])
MIDefault = MIBase(*([None] * len(MIBase._fields)))

# Extend from this to allow None as a default for some parameters. This is
# easier than subclassing to implement default args.
def MomentInfo(datatype, **kwargs):
    return MIDefault._replace(type=datatype, **kwargs)


class NetCDFRadarData(NetCDFData):
    def __init__(self, fname):
        NetCDFData.__init__(self, fname)

        self.radar = self.nc.RadarName

        runinfo = AttributeDict()
        self.runinfo = runinfo
        for attr in ['CantingWidth', 'ClearAir', 'ConfigFilename',
            'DataFiles', 'DataType', 'DropModel', 'FixedTemp',
            'GitDate', 'GitHash', 'GitTimeStamp', 'RadarName',
            'RandomSeed', 'RunStarted', 'ScatteringModel', 'SweepType',
            'VersionNumber']:
            runinfo[attr] = getattr(self.nc, attr)

        self.wavelength = self.readVar('Wavelength')[0]
        self.az = self.readVar('Azimuth')
        self.gate_length = self.readVar('GateLength')
        self.pulse_length = self.readVar('PulseDuration')[0] * units.c / 2.
        self.rng = (np.arange(len(self.nc.dimensions['gates']))
            * self.gate_length + self.readVar('RangeToFirstGate'))

        self.nyquist = self.readVar('NyquistVelocity')
        self.noise_pwr_db = self.readVar('MinimumDetectableSignal')
        self.noise_pwr = units.to_linear(self.noise_pwr_db)

        self.xlocs = self.rng * sin(self.az[:, np.newaxis])
        self.ylocs = self.rng * cos(self.az[:, np.newaxis])

        self.vel = self.readVar('Velocity')
        self.spw = self.readVar('SpectrumWidth')

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

        self.zdr = units.make_dB(self.ref_H - self.ref_V)
        self.fields[MomentInfo(datatypes.ZDR, source='average')] = self.zdr

        self.kdp = self.readVar('KDP')
        self.fields[MomentInfo(datatypes.KDP, source='average')] = self.kdp

        self.phidp = self.readVar('PhiDP')
#        self.phidp = remainder(self.phi_dp, 360.)
#        self.phidp2 = 2*self.kdp.cumsum(axis=1) * (self.pulse_length / kilo)
        self.fields[MomentInfo(datatypes.KDP, source='average')] = self.phidp

        thgrad, rgrad = np.gradient(self.phidp_ts, 1, self.gate_length)
        self.kdp_ts = rgrad / 2.
        self.fields[MomentInfo(datatypes.KDP, source='ts')] = self.kdp_ts

        self.diff_atten = (units.make_dB(self.unatten_pwr_H -
            self.unatten_pwr_V) - self.zdr)
        self.fields[MomentInfo(datatypes.DiffAtten,
            source='average')] = self.diff_atten

        # TODO: Need to read in the diagnostic variables.
        self.fields['x'] = self.xlocs
        self.fields['y'] = self.ylocs

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
