# =====================================================
# Universal Constants
# =====================================================
import numpy as np

G_km = 6.67430e-20  # km^3 kg^-1 s^-2
AU = 1.495978707e8        # Astronomical unit, km
day_sec = 86400           # Seconds in a day
year_sec = 365.25 * day_sec
c_m = 2.998e8 # m
c_km = 2.998e5 # km

# Conversion helpers
mil_to_cm = 0.00254
ft2_to_m2 = 0.092903
ft_to_m = 0.3048
lbf_to_N = 4.44822

# =====================================================
# Planetary Constants
# =====================================================

# --- Sun ---
mu_sun = 1.32712440018e11   # km^3/s^2 (GM of Sun)
m_sun = 1.98847e30          # kg
r_sun = 696340.0            # km (mean radius)

# --- Earth ---
mu_earth = 3.986004418e5    # km^3/s^2
m_earth = 5.97219e24        # kg
r_earth = 6378.137          # km (equatorial radius)
g0 = 9.80665                # m/s^2, standard gravity
v_esc_earth = 11.186        # km/s, escape velocity from surface
v_orbit_earth = 7.905       # km/s, circular orbit at LEO (~200 km)
obliquity_earth = np.deg2rad(23.439)
n_earth = 1.9909836747685184e-7

# --- Moon ---
mu_moon = 4902.800066       # km^3/s^2
m_moon = 7.342e22           # kg
r_moon = 1737.4             # km
a_moon = 384400.0           # km (Earth–Moon distance, semi-major axis)

# --- Mercury ---
mu_mercury = 2.2032e4       # km^3/s^2
m_mercury = 3.3011e23       # kg
r_mercury = 2439.7          # km

# --- Venus ---
mu_venus = 3.24859e5        # km^3/s^2
m_venus = 4.8675e24         # kg
r_venus = 6051.8            # km

# --- Mars ---
mu_mars = 4.282837e4        # km^3/s^2
m_mars = 6.4171e23          # kg
r_mars = 3389.5             # km

# --- Jupiter ---
mu_jupiter = 1.26686534e8   # km^3/s^2
m_jupiter = 1.8982e27       # kg
r_jupiter = 69911.0         # km

# --- Saturn ---
mu_saturn = 3.7931187e7     # km^3/s^2
m_saturn = 5.6834e26        # kg
r_saturn = 58232.0          # km

# --- Uranus ---
mu_uranus = 5.793939e6      # km^3/s^2
m_uranus = 8.6810e25        # kg
r_uranus = 25362.0          # km

# --- Neptune ---
mu_neptune = 6.836529e6     # km^3/s^2
m_neptune = 1.02413e26      # kg
r_neptune = 24622.0         # km

# =====================================================
# Other Useful Values
# =====================================================

# Standard Earth orbital parameters
earth_orbital_speed = 29.78  # km/s around the Sun
earth_orbital_period = 365.25 * day_sec  # s

# Sphere of influence radii (approximate, km)
SOI_earth = 924000.0
SOI_moon = 66000.0

# Propulsion
q = 1.6022e-19 # Coulombs; electric charge
amu = 1.6605e-27 # kg; atomic mass unit

