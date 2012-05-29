from numpy import pi
import quantities as pq
R = pq.constants.R

#Earth
Re = earth_avg_radius = 6.37e6 * pq.m
g = earth_avg_gravity = 9.81 * pq.m / pq.s**2
omega = earth_avg_angular_vel = 2 * pi * pq.rad / pq.day
d = earth_sfc_avg_dist_sun = 1.50e11 * pq.m
S = earth_solar_irradiance = 1.38e3 * pq.W / pq.m**2

spec_heat_units = pq.unit_registry['J / (K * kg)']

#Water
Mw = water_molecular_weight = 18.016 * pq.g / pq.mol
Rv = water_gas_constant = (R / Mw).rescale(spec_heat_units)
# Nominal density of liquid water at 0C
rho_l = density_water = 1e3 * pq.kg / pq.m**3
Cp_v = water_vapor_specific_heat_press = 1952. * spec_heat_units
Cv_v = water_vapor_specific_heat_vol = 1463. * spec_heat_units
Cp_l = water_specific_heat = 4218. * spec_heat_units # at 0C
Lv = water_latent_heat_vaporization = 2.5e6 * pq.J / pq.kg #0C
Lf = water_latent_heat_fustion = 3.34e5 * pq.J / pq.kg #0C
Cp_i = ice_specific_heat = 2106 * spec_heat_units # at 0C
rho_i = density_ice = 917 * pq.kg / pq.m**3 # at 0C

#Dry air
Md = dry_air_molecular_weight = 28.97 * pq.g / pq.mol # at the sfc
Rd = dry_air_gas_constant = (R / Md).rescale(spec_heat_units)
Cp_d = dry_air_spec_heat_press = 1004. * spec_heat_units
Cv_d = dry_air_spec_heat_vol = 717. * spec_heat_units
dry_air_density_stp = 1.275 * pq.kg / pq.m**3 # at 0C 1000mb

#General meteorology constants
kappa = poisson_exponent = (Rd / Cp_d).rescale('dimensionless')
gamma_d = dry_adiabatic_lapse_rate = (g / Cp_d).rescale('K/km')
epsilon = molecular_weight_ratio = Mw / Md

def theta_to_temperature(pot_temp, pressure):
    '''Calculates the temperature in K from the given values of of theta (K)
       and pressure (Pa).'''
    return pot_temp * (pressure / (100000. * pq.Pa))**kappa

def spec_humidity_to_pressure(spec_humid, pressure):
    '''Calculates the vapor pressure given the values of specific humidity and
       pressure.'''
    return spec_humid * pressure / (epsilon + (1 - epsilon) * spec_humid)

def air_density(temperature, pressure):
    '''Calculates air density given values of temperature and pressure.'''
    return (pressure / (dry_air_gas_constant * temperature)).simplified

def adjust_fall_speed(vt, air_density):
    '''Adjust terminal fall speeds for variations in air density.  The
    correction comes from Foote and Du Toit (1969).  The refernce density is
    for 20C and 1013mb pressure.'''
    reference_density = 1.204 * pq.kg / pq.m**3
    vt *= (reference_density / air_density).simplified**0.4

def mixing_ratio_to_density(mixing_ratio, air_density):
    '''Calculates density given values of mixing ratio and air density.'''
    return mixing_ratio * air_density
