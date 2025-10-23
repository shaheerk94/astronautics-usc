#!/usr/bin/env python3
"""
Exam Solutions – 9 Problems, Total 1000 points

Each problem has its own function (problem1() … problem9()).
The main() routine runs all problems sequentially.
"""

from datetime import datetime, timedelta
import numpy as np
import orbital_mechanics as omf
import space_functions as sf
import attitude_control as atc
import globals as g


# ------------------------------------------------------------
# Problem 1 (300 pt)
# ------------------------------------------------------------
def problem1():
    """Placeholder for Problem 1 logic."""
    print("Running Problem 1 (300 pt)…")
    Ps = 24.47 * 24 * 60 * 60 # s
    a = omf.semimajor_axis_from_orbital_period(g.mu_sun, Ps)
    print('a:', a)

    # part bi
    r0 = g.AU
    r1 = a
    i = 7.5 * np.pi/180
    a = (r0 + r1)/2
    v1  = omf.orbital_velocity(g.mu_sun, r0, r0)
    v2 = omf.orbital_velocity(g.mu_sun, r0, a)

    dv1 = np.sqrt(v1**2 + v2**2 - 2* v1 * v2 * np.cos(i))
    print('v1:', v1, 'v2:', v2, 'dv1:', dv1)

    v3 = omf.orbital_velocity(g.mu_sun, r1, a)
    v4 = omf.orbital_velocity(g.mu_sun, r1, r1)
    dv2 = v3 - v4
    print('v3:', v3, 'v4:', v4, 'dv2:', dv2)
    dvtotal = dv2 + dv1
    print('dvtotal:', dvtotal)

    # part bii
    c3 = dv1**2
    print('c3:', c3)

    # part biii
    Ttransfer = omf.orbital_period(g.mu_sun, a)/2
    print('Ttransfer:', Ttransfer)
    print('Ttransfer:', sf.format_time(Ttransfer))

    # part biv
    e = (r0 - r1)/(r0 + r1)
    b = a * np.sqrt(1-e**2)
    print(f"e: {e:.6f}, b: {b:.3e}")
    print('b:', b/g.AU, 'AU')

    result = None
    return result


# ------------------------------------------------------------
# Problem 2 (40 pt)
# ------------------------------------------------------------
def problem2():
    print("Running Problem 2 (40 pt)…")
    i = np.pi/2
    P0 = 90.5 * 60 # s
    a = omf.semimajor_axis_from_orbital_period(g.mu_earth, P0)
    print('a:', a)
    hpole = a - 6357
    heq = a - 6378
    print('hpole:', hpole, 'heq:', heq)
    return None


# ------------------------------------------------------------
# Problem 3 (60 pt)
# ------------------------------------------------------------
def problem3():
    print("Running Problem 3 (60 pt)…")
    mv = 0.03
    Mv = 0.58
    R = 10**((Mv - mv -5)/-5)
    print('Mv:', mv, 'Mv:', Mv, 'R:', R, 'parsec', 'R:', R*3.26156, 'ly')
    return None


# ------------------------------------------------------------
# Problem 4 (40 pt)
# ------------------------------------------------------------
def problem4():
    print("Running Problem 4 (40 pt)…")
    # TODO: Add your solution here
    return None


# ------------------------------------------------------------
# Problem 5 (50 pt)
# ------------------------------------------------------------
def problem5():
    print("Running Problem 5 (50 pt)…")
    r = 400 + g.r_earth
    V = omf.orbital_velocity(g.mu_earth, r, r) * 1000 # m/s
    m = 410 * 1000 # kg
    A = 1400 # m^2
    Cd = 2.02
    rho = 8e-15 * (100**3) / 1000 # kg/m^3
    ad = (-1/2) * Cd * A * rho * V**2/m
    print('ad:', ad)
    dv = ad * 365.24 * 24 * 3600
    print('dv:', dv)
    return None


# ------------------------------------------------------------
# Problem 6 (50 pt)
# ------------------------------------------------------------
def problem6():
    print("Running Problem 6 (50 pt)…")
    Bm = 10
    A = 8.0 * g.ft2_to_m2
    Cd = 2.15
    m = Bm * A * Cd
    print('A:', A, 'B:', Bm, 'Cd:', Cd, 'm:', m)
    return None


# ------------------------------------------------------------
# Problem 7 (60 pt)
# ------------------------------------------------------------
def problem7():
    print("Running Problem 7 (60 pt)…")
    ds = 440 # m
    A = np.pi *((ds/2)**2)
    cs = 0.88
    ca = 0.12
    R = 0.2 #km
    Fs = 1362.0 # W/m^2
    Fr = Fs / (R**2)
    fx, fy = atc.solar_radiation_force(0, ca, cs, 0, Fr, A)
    print('f:', fx)

    return None


# ------------------------------------------------------------
# Problem 8 (140 pt)
# ------------------------------------------------------------
def problem8():
    print("Running Problem 8 (140 pt)…")
    Vg = 6.6 # km/s
    Rv = 6371 # km
    P0 = 2 * np.pi * Rv/Vg
    a = omf.semimajor_axis_from_orbital_period(g.mu_earth, P0)
    print('P0:', P0, 'a:', a)
    return None


# ------------------------------------------------------------
# Problem 9 (260 pt)
# ------------------------------------------------------------
def problem9():
    print("Running Problem 9 (260 pt)…")
    TLE1 = "1 38752U 12046A 25251.97208406 .00356591 82059-6 22905-2 0 9998"
    TLE2 = "2 38752 9.6138 247.0804 5914676 302.4155 13.5078 4.29177350 137897"
    elements = omf.parse_tle(TLE1, TLE2)
    a = elements['semi_major_axis_km']
    b = a * np.sqrt(1 - elements['eccentricity']**2)
    rp = a * (1 - elements['eccentricity'])
    hp = rp - g.r_earth
    ra = a * np.sqrt(1 + elements['eccentricity'])
    vp = omf.orbital_velocity(g.mu_earth, rp, a)
    va = omf.orbital_velocity(g.mu_earth, ra, a)
    M = elements['mean_anomaly_deg']
    E = omf.mean_to_true_anomaly(M, elements['eccentricity'])
    print('b:', b, 'rp:', rp, 'hp:', hp, 'ra:', ra, 'vp:', vp, 'va:', va, 'E:', E)
    return None


# ------------------------------------------------------------
# Main driver
# ------------------------------------------------------------
def main():
    print("========== Exam Execution Start ==========")
    start_time = datetime.now()

    # problem1()
    problem3()
    problem4()
    problem5()
    problem6()
    problem7()
    problem8()
    problem9()
    # Execute all problems in order
    # results = {
    #     "Problem 1": problem1(),
    #     "Problem 2": problem2(),
    #     "Problem 3": problem3(),
    #     "Problem 4": problem4(),
    #     "Problem 5": problem5(),
    #     "Problem 6": problem6(),
    #     "Problem 7": problem7(),
    #     "Problem 8": problem8(),
    #     "Problem 9": problem9(),
    # }

    end_time = datetime.now()
    elapsed = end_time - start_time

    # print("========== Exam Execution Complete ==========")
    print(f"Total runtime: {elapsed}")
    print("\nResults Summary:")
    # for k, v in results.items():
    #     print(f"  {k}: {v}")

if __name__ == "__main__":
    main()
