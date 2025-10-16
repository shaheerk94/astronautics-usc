import numpy as np
import pint

import space_functions
import space_functions as sf
import orbital_mechanics as omf
import globals as g


ureg = pint.UnitRegistry()
## Homework 3

def archive():
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
    r1 = g.r_earth + omf.km_from_naut_miles(150.0)
    r2 = 42164 #km
    ic = space_functions.dms_to_decimal(28, 30)

    hohmann_transfer = omf.hohmann_with_inclination(r1, r2, g.mu_earth, ic)
    print(hohmann_transfer)

    ## Problem 10
    P2 = 98.5 * 60 # s
    r1 = 1000.0 + g.r_earth # km
    a2 = (((P2/(2*np.pi))**2)*g.mu_earth)**(1/3)
    h2 = a - g.r_earth

    print(h2)

    v = np.sqrt(g.mu_earth/a)

    inc = omf.sun_sync_inclination(a)
    print(inc * 180/np.pi)

    dv = omf.dv_inclination_change(v, v, inc-(np.pi/2))

    ## Problem 11
    P0 = 150 *60 # s
    a = omf.semimajor_axis_from_orbital_period(g.mu_earth, P0)
    print(a)

    # b
    e = np.roots([2.6, 2, -0.6])[1]
    ra = a * (1+e)
    va = omf.orbital_velocity(g.mu_earth, ra, a)
    print(ra, va)
    v2 = np.sqrt(va**2 + (0.5**2))
    a2 = (-1 * g.mu_earth / (((v2**2)/2) - (g.mu_earth/ra)))/2
    print(a2)
    p2 = omf.orbital_period(g.mu_earth, a2)
    print(p2)

    ## Problem 12
    rp = 800 + g.r_earth
    ra = 20000 + g.r_earth
    a = (rp + ra)/2
    e = 1 - (rp/a)
    f1 = (360 - 200)
    f2 = (160 + 180)
    e1 = omf.true_to_eccentric(f1, e)
    e2 = omf.true_to_eccentric(f2, e)
    m1 = omf.eccentric_to_mean(e1, e)
    m2 = omf.eccentric_to_mean(e2, e)
    t1 = omf.time_from_mean_anomaly(m1, a, g.mu_earth)
    t2 = omf.time_from_mean_anomaly(m2, a, g.mu_earth)
    t_total = omf.orbital_period(g.mu_earth, a)
    asc_to_desc = t2 -t1
    desc_to_asc = (t_total - t2) + t1
    print(f1, e1, m1, t1)
    print(f1, e1*180/np.pi, m1, t1)
    print(f2, e2, m2, t2)
    print(f2, e2*180/np.pi, m2, t2)
    print(t2 - t1)
    print(desc_to_asc)

def problem13():
    i = 97.63
    a = omf.sun_sync_semi_major_axis(incl_deg=i, e=0.00)
    h = a - g.r_earth
    print(a)
    print(h)
    m = 160 # kg
    A = 7.14 # m^2
    Cd = 2.24
    Bm = m/(Cd * A)
    print(Bm)

    rho = np.asin(g.r_earth/a)
    print(rho)
    print(rho * 180/np.pi)
    Dh = a * np.cos(rho)
    print(Dh)

    theta0 = np.acos(g.r_earth/a)
    Sw = 2 * theta0 * g.r_earth
    print(theta0, Sw)

    Tmax = 2 * theta0 * np.sqrt(a**3/g.mu_earth)
    P0 = omf.orbital_period(g.mu_earth, a)
    Tmax = 2 * theta0 / (2 * np.pi) * P0
    print(P0, Tmax)

    elev = 24 * np.pi/180
    n = np.asin(np.cos(elev) * np.sin(rho))
    print(n)
    print(n*180/np.pi)
    lamb = np.pi/2 - elev - n
    print(lamb)
    print(lamb * 180/np.pi)
    tmax2 = 2 * lamb * np.sqrt(a**3/g.mu_earth)
    print(tmax2)

def main():
    problem13()


if __name__ == "__main__":
    main()
