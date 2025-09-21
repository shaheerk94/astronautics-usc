import numpy as np
import orbital_mechanics as omf
from globals import mu_earth, r_earth

def test_convert_naut_miles():
    assert np.isclose(omf.km_from_naut_miles(1), 1.852)

def test_convert_miles():
    assert np.isclose(omf.km_from_miles(1), 1.60934)

def velocity_circular_orbit():
    r = r_earth + 100
    out = omf.velocity_circular_orbit(r, mu_earth)
    assert abs(out - 7.844108) < 1e-6

def test_orbital_period():
    a = 7000
    P = omf.orbital_period(mu_earth, a)
    # should be about 5828 s for LEO
    assert np.isclose(P, 5828, rtol=1e-3)

def test_orbital_velocity():
    a = 7000
    r = 7000
    v = omf.orbital_velocity(mu_earth, r, a)
    # circular velocity ~7.55 km/s
    assert np.isclose(v, 7.55, atol=0.05)

def test_perigee_velocity():
    e, p = 0.1, 7000
    vp = omf.perigee_velocity(e, mu_earth, p)
    assert vp > 0

def test_apogee_velocity():
    e, p = 0.1, 7000
    va = omf.apogee_velocity(e, mu_earth, p)
    assert va > 0
    assert va < omf.perigee_velocity(e, mu_earth, p)

def test_first_hohmann_dv():
    r1, r2 = 7000, 10000
    a = (r1 + r2)/2
    dv1 = omf.first_hohmann_dv(mu_earth, r1, a)
    assert dv1 > 0

def test_second_hohmann_dv():
    r1, r2 = 7000, 10000
    a = (r1 + r2)/2
    dv2 = omf.second_hohmann_dv(mu_earth, r2, a)
    assert dv2 < 0  # by construction in your formula

def test_total_hohman_dv():
    r1, r2 = 7000, 10000
    at = (r1 + r2)/2
    dv_tot = omf.total_hohman_dv(mu_earth, r1, r2, at)
    assert dv_tot > 0

def test_time_of_flight_hohmann():
    r1, r2 = 7000, 10000
    at = (r1 + r2)/2
    tof = omf.time_of_flight_hohmann(at, mu_earth)
    assert tof > 0

def test_dv_inclination_change():
    v, thetad = 7.5, 30
    dv = omf.dv_inclination_change(v, thetad)
    assert np.isclose(dv, 2*v*np.sin(np.radians(thetad)/2))

def test_dtheta_chaser():
    ac, at = 7000, 10000
    dtheta = omf.dtheta_chaser(ac, at)
    assert dtheta > 0

def test_arc_length_from_phase_change():
    ac, at = 7000, 10000
    ds = omf.arc_length_from_phase_change(ac, at)
    assert ds > 0

def test_angular_speed():
    w = omf.angular_speed(mu_earth, 7000)
    assert np.isclose(w, 0.001078, rtol=1e-3)

def test_phase_angle_rendezvous_begin():
    ac, at = 7000, 10000
    thetaf = omf.phase_angle_rendezvous_begin(ac, at)
    assert thetaf > 0

def test_wait_time():
    WT = omf.wait_time(1.0, 0.5, 0.001, 0.0009, 1)
    assert isinstance(WT, float) or isinstance(WT, np.ndarray)

def test_time_of_flight_rendezvous():
    ac, at = 7000, 10000
    tof = omf.time_of_flight_rendezvous(ac, at, mu_earth)
    assert tof > 0

def test_ac_coplanar():
    ac = omf.ac_coplanar(mu_earth, 0.5, 0.001)
    assert ac > 0

def test_dv_initial_coplanaar():
    ac, at = 7000, 10000
    dv = omf.dv_initial_coplanaar(mu_earth, at, ac)
    assert isinstance(dv, float) or isinstance(dv, np.ndarray)

def test_dv_total_coplanar():
    ac, at = 7000, 10000
    dv_tot = omf.dv_total_coplanar(mu_earth, ac, at)
    assert dv_tot >= 0
