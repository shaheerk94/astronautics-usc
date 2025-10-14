import numpy as np
import pint
import space_functions as sf
import orbital_mechanics as omf
import globals as g


ureg = pint.UnitRegistry()
## Homework 3

rs = 8.71 * ureg.light_year
s_mag = sf.apparent_magnitude(rs, 1.41)
print(s_mag)

## homework 4
# part 1
deg_day = 360/365.24
days_equinox_to_perihelion = 102.94 * deg_day

# Part 2
Bs = -39.61
lambda_s = 104.08
ra_c, dec_c = 219.90, -60.84

ra_s, dec_s = sf.ecliptic_to_celestial(lambda_s, Bs)
print(ra_s, dec_s)

lambda_c, Bc = sf.celestial_to_ecliptic(ra_c, dec_c)

angle = sf.planet_angle(lambda_c, Bc, lambda_s, Bs)
print(angle * 180/np.pi)

## homework 5
# 2b
bm, cd, a = 40, 2.15, 14 * g.ft2_to_m2
m = bm * cd * a
print(m)

# 3
m_a = 2.5 # g/cm^2, areal density
rho_al = 2.70 # g/cm^2
t_cm = 2.5/rho_al # cm
t_mil = t_cm / g.mil_to_cm
print(t_mil)

# 3b

#TLE
TLE1 = "1 00005U 58002B 06352.91214340 +.00000075 +00000-0 +98410-4 0 01338"
TLE2 = "2 00005 034.2426 280.7256 1851086 148.7108 223.8259 10.83943630670210"
elements = omf.parse_tle(TLE1, TLE2)
print(elements)