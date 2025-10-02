import numpy as np
import space_functions as sf
import orbital_mechanics as omf
import globals as g

## Problem 12

def problem_12():
    ra = 20000 + g.r_earth # apogee distance
    rp = 800 + g.r_earth # perigee distance
    a = (rp + ra)/2
    e = (ra - rp)/(2*a)

    f_asc = 360 - 200
    f_desc = f_asc + 180

    E_asc = omf.true_to_eccentric(f_asc, e)
    E_desc = omf.true_to_eccentric(f_desc, e)

    E_asc_deg = np.rad2deg(E_asc)
    E_desc_deg = np.rad2deg(E_desc)

    M_asc = omf.eccentric_to_mean(E_asc, e)
    M_desc = omf.eccentric_to_mean(E_desc, e)

    t_asc = omf.time_from_mean_anomaly(M_asc, a, g.mu_earth)
    t_desc = omf.time_from_mean_anomaly(M_desc, a, g.mu_earth)

    t_total = omf.orbital_period(g.mu_earth, a)
    t_asc_to_desc = t_desc - t_asc
    t_desc_to_asc = (t_total - t_desc) + t_asc


def problem_13():
    # constants
    mu = 398600.4418  # km^3/s^2
    Re = 6378.137  # km (equatorial radius)
    J2 = 1.08262668e-3
    i = np.radians(97.63)  # inclination in radians

    # sun apparent motion (rad/s)
    deg_per_day = 360.0 / 365.2422  # deg/day
    sun_rate = np.radians(deg_per_day) / 86400.0  # rad/s

    # solve for semi-major axis a
    a = ((1.5 * J2 * np.sqrt(mu) * Re ** 2 * abs(np.cos(i))) / sun_rate) ** (2 / 7)
    h = a - Re

    print(a, h)

    # Part B
    m = 160 # kg
    A = 7.14 # m^2
    Cd = 2.24

    Bm = m/(Cd * A)

    ## Part C
    Rv = 6371.0 # km
    r = Rv + h
    dh = np.sqrt(r**2 - Rv**2)
    print(dh)

    ## Part D
    theta_0 = np.acos(Rv/r)
    print(theta_0)

    Sw = 2 * theta_0 * Rv
    print(Sw)

    ## Part E
    Tmax = 2 * theta_0 * np.sqrt(r**3 / g.mu_earth)
    print(Tmax/60)

    ## Part F
    n = np.sqrt(g.mu_earth/ (r**3))
    vg = n * Rv
    print(vg)

    ## part G
    e = np.deg2rad(24)
    rho = np.asin(Rv/r)
    n = np.asin(np.cos(e) * np.sin(rho))

    theta = (np.pi/2) - n - e
    Tmax2 = 2 * theta * np.sqrt(r**3 / g.mu_earth)
    print(Tmax2/60)


def problem_14():
    # Part a
    P = 86164.0905 # s

    # Part b
    a = (g.mu_earth * ((P/(2 * np.pi))**2)) **(1/3)
    r_geo = a
    h = a - g.r_earth
    print(f'a: {a} km')
    print(f'h: {h} km')

    # Part C
    theta = np.sin(g.r_earth/a)
    print('theta:', theta*2)
    print(f'theta: {theta * 180/np.pi * 2} deg')

    # Part D
    theta_0 = np.acos(g.r_earth/a)
    print('theta_0:', theta_0*180/np.pi, 'deg')

    # Part F
    theta = sf.dms_to_decimal(61, 13, 0, 'N')

    x = np.sqrt(g.r_earth**2 + r_geo**2 - (2 * g.r_earth * r_geo * np.cos(theta * np.pi / 180)))
    print('x:', x, 'km')

    y = ((r_geo**2) - (g.r_earth**2) - (x**2))/(-2 * g.r_earth * x)
    e = np.acos(y) - np.pi/2
    print('e:', e*180/np.pi, 'deg')

    # Part G
    c = 299774 # km/s
    T = 2*x/c
    print(f'T: {T} s')


def main():
    # problem_12()
    # problem_13()
    problem_14()


if __name__ == '__main__':
    main()



