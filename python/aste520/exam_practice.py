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
TLE2 = "2 00005 034.2426 280.7256 1851086 148.7108 223.8259 10.83943630 670210"
elements = omf.parse_tle(TLE1, TLE2)
print(elements)

## Problem 6
# Part 1
tm = 0.5 * 365.24 * 24 * 3600 # s
V0 = 7226 # m/s
N0 = 1.1e9 * (100**3) # m^-3
F0 = N0 * V0 * tm
F0_cm = F0 / (100**2) # cm^-2
print(f"F0 = {F0_cm:.3e}")

# Part 2
ve = omf.orbital_velocity(g.mu_sun, g.AU, g.AU)
print(ve)

# Part 3
V2 = 2 * g.mu_earth/g.AU
V1 = ve


## problem 7
hp = 600.0 # n.m.
e = 0.820

rp = omf.km_from_naut_miles(hp) + g.r_earth

## Problem 8
TLE1 = "1 03955U 69046E 24105.44078525 -.00000639 00000-0 00000-0 0 9993"
TLE2 = "2 03955 34.3121 112.7933 8505845 114.4265 359.4801 0.21375648 5306"
elements = omf.parse_tle(TLE1, TLE2)
print(elements)

## Problem 9

