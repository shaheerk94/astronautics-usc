import numpy as np
import space_functions as sf
import orbital_mechanics as omf
import globals as g

## Problem 9

# Part A
h = omf.km_from_naut_miles(150.0)
r1 = h + g.r_earth
r_geo = 42164 # km
delta_ic = sf.dms_to_decimal(28, 30)
delta_ic_rad = delta_ic * np.pi / 180

print(f'Height: {h} km, init oribit radius: {r1}')
print(f'Inclination Change: {delta_ic} degrees')
print(f'Inclination Change: {delta_ic_rad} rads')

v1 = omf.orbital_velocity(g.mu_earth, r1, r1)
print(f'Initial Velocity of Circular orbit: {v1}')

rp2 = r1 # perigee radius = circular orbit radius
ra2 = r_geo # apogee radius = geo orbit radius
a2 = (rp2 + ra2)/2 # semi-major axis of transfer orbit
vp2 = omf.orbital_velocity(g.mu_earth, rp2, a2)
print(f'Perigee Velocity Transfer Orbit: {vp2}')

dv1 = vp2 - v1
print(f'DV first burn: {dv1}')

