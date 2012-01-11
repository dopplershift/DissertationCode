import numpy as np
from numpy.ma import masked_array
from scipy.constants import degree, kilo, centi
import netCDF4
import quantities as pq
from quantities import sin, cos
import .units
from .sigproc import auto_moments


class AttributeDict(dict): 
    __getattr__ = dict.__getitem__
    __setattr__ = dict.__setitem__

class DataSet(object):

    def __init__(self):
        self.fields = []
        self._calc_cache = dict()
        self.metadata = []
        self.coordinates = []

#    def getField(datatype, prefer=None):
#        if datatype in self.fields:
#            return self.fields[datatype]
#        elif (datatype,prefer) in self._calc_cache:
#            return self._calc_cache[datatype, prefer]
#        else:
#            alg = getBestAlg(datatype, self.fields.keys())
#            return alg(self)

class NetCDFData(DataSet):
    def __init__(self, fname
        DataSet.__init(self)
        nc = netCDF4.Dataset(fname, 'r')
        self.nc = nc

    def readVar(varname):
        var = self.nc.variables[varname]
        vals = var[:]
        # Check if we were given a masked array. For now, we'll fill masks
        # with NaNs
        if(hasattr(vals, 'fill_value')):
            vals.fill_value = np.nan
            vals = vals.filled()

        try:
            return pq.Quantity(vals, var.units)
        except LookupError:
            return vals

class NetCDFRadarData():

    def __init__(self, fname):
        NetCDFData.__init__(self, fname)

        self.radar = nc.RadarName

        runinfo = AttributeDict()
        self.runinfo = runinfo
        for attr in ['CantingWidth', 'ClearAir', 'ConfigFilename',
            'DataFiles', 'DataType', 'DropModel', 'FixedTemp',
            'GitDate', 'GitHash', 'GitTimeStamp', 'RadarName',
            'RandomSeed', 'RunStarted', 'ScatteringModel', 'SweepType',
            'VersionNumber']:
            runinfo[attr] = getattr(nc, attr)

        self.wavelength = readVar('Wavelength')[0]
        self.az = readVar('Azimuth')
        self.gate_length = readVar('GateLength')
        self.pulse_length = readVar('PulseDuration')[0] * c / 2.
        self.rng = (np.arange(len(nc.dimensions['gates']))
            * self.gate_length + readVar('RangeToFirstGate'))

        self.nyquist = readVar('NyquistVelocity')
        self.noise_pwr_db = readVar('MinimumDetectableSignal')
        self.noise_pwr = to_linear(self.noise_pwr_db)

        self.xlocs = self.rng * sin(self.az[:, np.newaxis])
        self.ylocs = self.rng * cos(self.az[:, np.newaxis])

        self.vel = readVar('Velocity')
        self.spw = readVar('SpectrumWidth')

        self.process_channel('H')
        self.process_channel('V')
        self.zdr_ts, self.rhohv, self.phidp_ts = auto_dual_pol(self.H.iq,
            self.V.iq, self.noise_pwr, self.noise_pwr)
        self.zdr_ts *= units.dB
        self.phidp_ts *= pq.radians

        self.zdr = make_dB(self.H.ref - self.V.ref)
        self.kdp = readVar('KDP')
        self.phidp = readVar('PhiDP')
#        self.phidp = remainder(self.phi_dp, 360.)
#        self.phidp2 = 2 * self.kdp.cumsum(axis=1) * (self.pulse_length / kilo)
        thgrad, rgrad = np.gradient(self.phidp_ts, 1, self.gate_length)
        self.kdp_ts = rgrad / 2.
        self.diff_atten = (make_dB(self.H.unatten_pwr - self.V.unatten_pwr)
            - self.zdr)

        def process_channel(self, pol):
            # Create container to hold all the moments
            data = Bag()
            setattr(self, pol, data)

            # Read moments directly from file
            data.ref = readVar('Reflectivity_%s' % pol)
            data.pwr = dBW_to_dBm(readVar('Power_%s' % pol))
            data.unatten_pwr = dBW_to_dBm(readVar('UnattenPower_%s' % pol))

            # Read time series data and calculate moments
            data.iq = (self.nc.variables['I_%s' % pol][:]
                + 1.0j * self.nc.variables['Q_%s' % pol][:])
            data.pwr_ts, data.vel_ts, data.spw_ts = auto_moments(data.iq,
                self.nyquist)
            data.pwr_ts_dbm = to_dBm(data.pwr_ts)

            # Convert time series power estimate to reflectivity factor
            data.rad_const_dB = readVar('RefCalibration_%s' % pol)
            data.rad_const = to_linear(data.rad_const_dB)
            data.ref_ts = to_dBz(data.pwr_ts * self.rng**2 * data.rad_const)

            # Calculate attenuation from fields
            data.atten = make_dB(self.unatten_pwr - self.pwr)

