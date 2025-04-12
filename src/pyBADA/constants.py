"""
Common constants module
"""

# Gravitational acceleration
g = 9.80665

# Real gas constant for air
R = 287.05287

# ISA temperature gradient with altitude below the tropopause
temp_h = 0.0065

# Standard atmospheric temperature at MSL
temp_0 = 288.15

# Standard atmospheric temperature at tropopause (11km)
temp_11 = 216.65

# Tropopause geopotatial pressure altitude [m]
h_11 = 11000.0

# Standard atmospheric pressure at tropopause (11km)
p_11 = 22632.04

# Standard atmospheric pressure at MSL
p_0 = 101325.0

# Standard atmospheric density at MSL
rho_0 = 1.225

# Speed of sound
a_0 = 340.294

# Adiabatic index of air
Agamma = 1.4

# Adiabatic index of air - ratio
Amu = (Agamma - 1) / Agamma

# Average Earth radius [km] - https://en.wikipedia.org/wiki/Earth_radius#Mean_radius
AVG_EARTH_RADIUS_KM = 6371

# semi-major axis (WGS-84) [m]
a = 6378137.0

# semi-minor axis (WGS-84) [m]
b = 6356752.314245

# flattening (WGS-84)
f = (a - b) / a
