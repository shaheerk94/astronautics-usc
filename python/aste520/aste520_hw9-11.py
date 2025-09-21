import numpy as np
import space_functions as sf
import orbital_mechanics as omf
import globals as g
from launch_sites import LAUNCH_SITES

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

## Circ orbit 1 --> Transfer Orbit
rp2 = r1 # perigee radius = circular orbit radius
ra2 = r_geo # apogee radius = geo orbit radius
a2 = (rp2 + ra2)/2 # semi-major axis of transfer orbit
vp2 = omf.orbital_velocity(g.mu_earth, rp2, a2)
print(f'Perigee Velocity Transfer Orbit: {vp2}')

dv1 = vp2 - v1
print(f'DV first burn: {dv1}')

## Inclination change @ apogee of transfer orbit and circularize
r3 = r_geo
va2 = omf.orbital_velocity(g.mu_earth, ra2, a2)
vr3 = omf.orbital_velocity(g.mu_earth, r3, r3)
dv2 = omf.dv_inclination_change(va2, vr3, delta_ic_rad)

print(f'Final circular velocity: {vr3} km/s; dv2 = {dv2} km/s')

## total DV
dv = dv1 + dv2
print(f'Total DV: {dv} km/s')

## Part C
sites = [
    LAUNCH_SITES['KOUROU'],
    LAUNCH_SITES['BAIKONUR'],
    LAUNCH_SITES['SRIHARIKOTA'],
    LAUNCH_SITES['TANEGASHIMA']
]

for site in sites:
    hohmann_transfer_dv = omf.hohmann_with_inclination(r1, r_geo, g.mu_earth, site['latitude_deg'])
    print(f"For site {site['name']} with inclination {site['latitude_deg']}:\n\t dv = {hohmann_transfer_dv['dv_total']} km/s")

