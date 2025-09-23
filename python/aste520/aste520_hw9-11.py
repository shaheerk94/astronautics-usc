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
    print(f"For site {site['name']} with inclination {site['latitude_deg']}:\n\t "
          f"dv = {hohmann_transfer_dv['dv_total']} km/s")


## Problem 10
h1 = 1000 # km
r1 = h1 + g.r_earth
P2 = 98.5 * 60 # orbital period, seconds

r2 = (((P2/(2*np.pi))**2) * g.mu_earth) **(1/3)
h2 = r2 - g.r_earth
print(f'Final orbit radius: {r2} km; altitude = {h2} km')

n = 2 * np.pi / P2
n_deg = n * 180/np.pi
print(f'N = {n_deg} degrees/s, N = {n} rad/s')

v = np.sqrt(g.mu_earth/r2)
print(f'Final orbit velocity: {v} km/s')

dvi = 2 * v * np.sin(8.13 * np.pi / 180 / 2)
print(f'DV for inclination change: {dvi} km/s')

# hohmann transfer
dv_hohmann = omf.hohmann_dv(g.mu_earth, r1, r2)
print(f'dv_hohmann: dv1 = {dv_hohmann["dv1"]}, dv2 = {dv_hohmann["dv2"]}, total dv = {dv_hohmann["dv_total"]} km/s')

## Homework 11
P0 = 150 * 60 # seconds
K = 1.6
dv = 0.5 # km/s

a1 = (g.mu_earth * ((P0/(2*np.pi))**2)) ** (1/3)
print(f'a1 = {a1} km')

e = np.roots([2.6, 2, -0.6])[1]

rp1 = a1 * (1-e)
ra1 = a1 * (1+e)
va1 = omf.orbital_velocity(g.mu_earth, ra1, a1)
v_total = np.sqrt(va1**2 + dv**2)
a2 = -g.mu_earth / ((v_total**2)/2 - (g.mu_earth/ra1)) / 2
P2 = 2 * np.pi * np.sqrt((a2**3) / g.mu_earth)
dp = P2 - P0

h = ra1 * va1
print(f'Angular momentum: {h} km^2/s')
e2 = np.sqrt(1 - (h**2) / (g.mu_earth * a2))
print(f'e2 = {e2}')
rp2 = a2 * (1 - e2)
hp2 = rp2 - g.r_earth
ra2 = a2 * (1 + e2)
ha2 = ra2 - g.r_earth
print(f'Perigee radius: {rp2} km; Perigee altitude: {hp2} km')
print(f'Apogee radius: {ra2} km; Apogee altitude: {ha2} km')
