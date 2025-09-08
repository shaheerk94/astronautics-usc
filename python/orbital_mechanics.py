import numpy as np

def convert_naut_miles(miles):
    return miles * 1.852

def convert_miles(miles):
    return miles * 1.60934

def circular_delta_a(a, u, delta_v):
    return 2 * np.sqrt(a**3 / u) * delta_v

def orbital_period(u, a):
    return 2 * np.pi / np.sqrt(u) * (a**1.5)

def orbital_velocity(u, r, a):
    return np.sqrt(u * (2/r - 1/a))

def perigee_velocity(e, u, p):
    return (1 + e) * np.sqrt(u/p)

def apogee_velocity(e, u, p):
    return (1 - e) * np.sqrt(u/p)

def first_hohmann_dv(u, r1, a):
    return np.sqrt(u * (2/r1 - 1/a))

def second_hohmann_dv(u, r2, a):
    return np.sqrt(u * (2/r2 - 1/a)) - np.sqrt(u/r2)

def total_hohman_dv(u, r1, r2, at):
    return (np.sqrt(u/r2) - np.sqrt(u * (2/r2 - 1/at)) +
            np.sqrt(u * (2/r1 - 1/at)) - np.sqrt(u/r1))

def time_of_flight_hohmann(at, u):
    return np.pi * np.sqrt(at**3 / u)

def dv_inclination_change(v, thetad):
    return 2 * v * np.sin(np.radians(thetad) / 2)

def dtheta_chaser(ac, at):
    return 2 * np.pi * (1 - (ac/at)**1.5)

def arc_length_from_phase_change(ac, at):
    return 2 * np.pi * at * (1 - (ac/at)**1.5)

def angular_speed(u, a):
    return np.sqrt(u / (a**3))

def phase_angle_rendezvous_begin(ac, at):
    return np.pi * (1 - ((ac + at) / (2*at))**1.5)

def wait_time(thetaf, thetai, wt, wc, k):
    return (thetaf - thetai) / (wt - wc) + (2 * np.pi * k) / (wt - wc)

def time_of_flight_rendezvous(ac, at, u):
    atrans = (ac + at) / 2
    return np.pi * np.sqrt(atrans**3 / u)

def ac_coplanar(u, thetai, wt):
    return (u * ((thetai / (2 * np.pi * wt))**2))**(2/3)

def dv_initial_coplanaar(u, at, ac):
    return np.sqrt(u * (2/at - 1/ac)) - np.sqrt(u/at)

def dv_total_coplanar(u, ac, at):
    return 2 * np.abs(np.sqrt(u * (2/at - 1/ac)) - np.sqrt(u/at))
