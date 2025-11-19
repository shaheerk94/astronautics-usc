import numpy as np
import globals as g
import orbital_mechanics as omf

def problem_21():
    Mp = 40 # kg
    t_day = 50 * 60  # s/day
    T_mission = 12 * 365 # days
    print('## Homework 17 ##')

    # Part A
    print('## Part A')
    total_firing_time = t_day * T_mission # total firing duration in minutes
    mass_flow_rate = Mp / total_firing_time

    print(f'Total firing time: {total_firing_time:E} seconds')
    print(f'Mass flow rate: {mass_flow_rate:E} kg/s')

    # Part B
    Pb = 4500 # Watts
    mxe = 131.0 * g.amu # kg

    Ua = mxe * Pb/(g.q * mass_flow_rate)
    print(f'Ua: {Ua:E} v')

    ## Part C
    ve = np.sqrt(2 * Pb / mass_flow_rate)
    print(f'Ve: {ve:E} m/s')
    isp = ve / g.g0
    print(f'Isp: {isp:E}')



def problem_22():
    print('## Homework 18 Part 22##')
    hp = 960 #km
    f = 2.8e9 # Hz
    P = 32 # W
    n = 0.54
    T = 12 * 3600 # 12 hours in seconds

    # Part a
    rp = hp + g.r_earth
    a = omf.semimajor_axis_from_orbital_period(g.mu_earth, T)
    print(f'a: {a} km')
    e = 1 - (rp/a)
    print(f'e: {e}')
    ra = a * (1+e)
    print(f'ra: {ra} km')

    rho = np.asin(g.r_earth/ra)
    print(f'rho: {rho} rad')
    print(f'rho: {rho*180/np.pi} deg')

    # Part B
    omega = 2 * np.pi * (1 - np.cos(rho))
    print(f'omega: {omega} rad')
    print(f'omega: {omega*180/np.pi} deg')

    # Part C
    Gd = 4 * np.pi/omega
    print(f'Gd: {Gd}')
    print(f'Gd: {10 * np.log10(Gd)} db')

    # Part D
    Gp = n * Gd
    print(f'Gp: {Gp}')
    print(f'Gp: {10 * np.log10(Gp)} db')


    # Part E
    wv = g.c_m/f # m
    print(f'wv: {wv} m')
    D = 1.22 * wv/0.30
    print(f'D: {D} m')

    # Part F
    dr = 4.0
    n = 0.52
    la_db = 6 # DB
    Gdr = (np.pi * dr/wv) **2
    print(f'Gd: {Gdr}')
    print(f'Gd: {10 * np.log10(Gdr)} db')
    Grp = n * Gdr
    print(f'Grp: {Grp}')
    print(f'Grp: {10 * np.log10(Grp)} db')

    # Part G
    la = 10 ** (-6/10)
    print(f'la: {la}')
    h_m = (ra - g.r_earth) * 1000
    Pr = (P * Gp * la)/(4* np.pi * ((h_m)**2))
    print(f'Pr: {Pr} W')


def main():
    # problem_21()
    problem_22()

if __name__ == "__main__":
    main()