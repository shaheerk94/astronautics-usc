import numpy as np
from datetime import datetime, timedelta
from scipy.optimize import fsolve

# =====================================================
# Unit Conversions
# =====================================================

def km_from_naut_miles(miles):
    """
    Convert nautical miles to kilometers.

    Args:
        miles (float): Distance in nautical miles.

    Returns:
        float: Distance in kilometers.
    """
    return miles * 1.852


def km_from_miles(miles):
    """
    Convert miles to kilometers.

    Args:
        miles (float): Distance in miles.

    Returns:
        float: Distance in kilometers.
    """
    return miles * 1.60934


# =====================================================
# Orbital Basics
# =====================================================

def orbital_period(u, a):
    """
    Compute orbital period.

    Args:
        u (float): Gravitational parameter μ (km^3/s^2).
        a (float): Semi-major axis (km).

    Returns:
        float: Orbital period (s).
    """
    return 2 * np.pi * np.sqrt((a ** 3)/u)


def semimajor_axis_from_orbital_period(u, T):
    """
    Compute semimajor axis from orbital period.
    Args:
        u (float): Gravitational parameter μ (km^3/s^2).
        a (float): Semi-major axis (km).

    Returns:
        float: Orbital period (s).
    """
    return (((T/(2*np.pi))**2)*u)**(1/3)


def orbital_velocity(u, r, a):
    """
    Orbital velocity from the vis-viva equation.

    Args:
        u (float): Gravitational parameter μ (km^3/s^2).
        r (float): Distance from central body (km).
        a (float): Semi-major axis (km).

    Returns:
        float: Velocity (km/s).
    """
    return np.sqrt(u * (2 / r - 1 / a))


def perigee_velocity(e, u, p):
    """
    Velocity at perigee of an elliptical orbit.

    Args:
        e (float): Eccentricity.
        u (float): Gravitational parameter μ (km^3/s^2).
        p (float): Semi-latus rectum (km).

    Returns:
        float: Perigee velocity (km/s).
    """
    return (1 + e) * np.sqrt(u / p)


def apogee_velocity(e, u, p):
    """
    Velocity at apogee of an elliptical orbit.

    Args:
        e (float): Eccentricity.
        u (float): Gravitational parameter μ (km^3/s^2).
        p (float): Semi-latus rectum (km).

    Returns:
        float: Apogee velocity (km/s).
    """
    return (1 - e) * np.sqrt(u / p)


def angular_speed(u, a):
    """
    Mean angular speed of an orbit.

    Args:
        u (float): Gravitational parameter μ (km^3/s^2).
        a (float): Semi-major axis (km).

    Returns:
        float: Angular speed (rad/s).
    """
    return np.sqrt(u / (a ** 3))


def mean_to_true_anomaly(M_deg, e, tol=1e-10, max_iter=50):
    """
    Convert mean anomaly to true anomaly for an elliptical orbit.

    Args:
        M_deg (float): Mean anomaly in degrees
        e (float): Eccentricity (0 <= e < 1)
        tol (float): Convergence tolerance (radians)
        max_iter (int): Maximum iterations for Newton's method

    Returns:
        float: True anomaly in degrees (0–360)
    """
    # Convert mean anomaly to radians, normalize to [-pi, pi]
    M = np.radians(M_deg) % (2 * np.pi)
    if M > np.pi:
        M -= 2 * np.pi

    # Initial guess
    E = M if e < 0.8 else np.pi  # rule of thumb

    # Newton-Raphson iteration
    for _ in range(max_iter):
        f = E - e * np.sin(E) - M
        fprime = 1 - e * np.cos(E)
        delta = f / fprime
        E -= delta
        if abs(delta) < tol:
            break

    # Convert eccentric anomaly to true anomaly
    nu = 2 * np.arctan2(np.sqrt(1 + e) * np.sin(E / 2), np.sqrt(1 - e) * np.cos(E / 2))
    nu_deg = np.degrees(nu) % 360.0
    return nu_deg


def true_to_eccentric(f_deg, e):
    """
    Convert true anomaly to eccentric anomaly robustly.

    Args:
        f_deg (float): True anomaly in degrees.
        e (float): Eccentricity (0 <= e < 1).

    Returns:
        float: Eccentric anomaly in radians (0 to 2π).
    """
    f = np.radians(f_deg)
    cosE = (e + np.cos(f)) / (1 + e * np.cos(f))
    sinE = (np.sqrt(1 - e ** 2) * np.sin(f)) / (1 + e * np.cos(f))
    E = np.arctan2(sinE, cosE)
    if E < 0:
        E += 2 * np.pi
    return E


def eccentric_to_mean(E, e):
    """
    Convert eccentric anomaly to mean anomaly.

    Args:
        E (float): Eccentric anomaly [radians].
        e (float): Eccentricity.

    Returns:
        float: Mean anomaly [radians], normalized to [0, 2π).
    """
    M = E - e * np.sin(E)
    # wrap into [0, 2π)
    M = np.mod(M, 2*np.pi)
    return M


def time_from_mean_anomaly(M, a, u):
    return np.sqrt((a**3)/u) * M

# =====================================================
# Hohmann Transfers
# =====================================================

def first_hohmann_dv(u, r1, a):
    """
    Delta-V for the first burn of a Hohmann transfer.

    Args:
        u (float): Gravitational parameter μ (km^3/s^2).
        r1 (float): Initial orbit radius (km).
        a (float): Semi-major axis of transfer orbit (km).

    Returns:
        float: Delta-V (km/s).
    """
    return np.sqrt(u * (2 / r1 - 1 / a))


def second_hohmann_dv(u, r2, a):
    """
    Delta-V for the second burn of a Hohmann transfer.

    Args:
        u (float): Gravitational parameter μ (km^3/s^2).
        r2 (float): Final orbit radius (km).
        a (float): Semi-major axis of transfer orbit (km).

    Returns:
        float: Delta-V (km/s).
    """
    return np.sqrt(u * (2 / r2 - 1 / a)) - np.sqrt(u / r2)


import numpy as np

def hohmann_dv(mu: float, r1: float, r2: float):
    """
    Two-burn Hohmann transfer ΔV from circular radius r1 to r2.

    Args:
        mu (float): Gravitational parameter μ (km^3/s^2).
        r1 (float): Initial circular-orbit radius (km).
        r2 (float): Final circular-orbit radius (km).

    Returns:
        dict: {
            'dv1': ΔV at perigee (km/s),
            'dv2': ΔV at apogee (km/s),
            'dv_total': dv1 + dv2 (km/s)
        }
    """
    a_t = 0.5 * (r1 + r2)

    v1 = np.sqrt(mu / r1)                       # initial circular speed
    v2 = np.sqrt(mu / r2)                       # final circular speed
    v_p = np.sqrt(mu * (2/r1 - 1/a_t))          # transfer speed at perigee
    v_a = np.sqrt(mu * (2/r2 - 1/a_t))          # transfer speed at apogee

    dv1 = abs(v_p - v1)
    dv2 = abs(v2 - v_a)
    return {'dv1': dv1, 'dv2': dv2, 'dv_total': dv1 + dv2}



def time_of_flight_hohmann(at, u):
    """
    Time of flight for a Hohmann transfer.

    Args:
        at (float): Semi-major axis of transfer orbit (km).
        u (float): Gravitational parameter μ (km^3/s^2).

    Returns:
        float: Transfer time (s).
    """
    return np.pi * np.sqrt(at ** 3 / u)


def hohmann_with_inclination(r1: float, r2: float, mu: float, delta_i_deg: float = 0.0):
    """
    Compute Δv for a two-burn Hohmann transfer with optional inclination change
    at apogee of the transfer ellipse.

    Args:
        r1 (float): Radius of initial circular orbit (km).
        r2 (float): Radius of final circular orbit (km).
        mu (float): Gravitational parameter μ (km^3/s^2).
        delta_i_deg (float): Inclination change at apogee (degrees). Default 0.

    Returns:
        dict: {
            'dv1': first burn at perigee (km/s),
            'dv2': second burn at apogee including plane change (km/s),
            'dv_total': total Δv (km/s),
            'va': velocity at apogee before burn (km/s),
            'vcirc2': circular velocity at r2 (km/s)
        }
    """
    delta_i = np.deg2rad(delta_i_deg)

    # Semi-major axis of transfer ellipse
    a_t = 0.5 * (r1 + r2)

    # Initial circular speed at r1
    v1_circ = np.sqrt(mu / r1)

    # Transfer ellipse perigee velocity
    v_p = np.sqrt(mu * (2 / r1 - 1 / a_t))

    # First burn Δv
    dv1 = v_p - v1_circ

    # Transfer ellipse apogee velocity
    v_a = np.sqrt(mu * (2 / r2 - 1 / a_t))

    # Final circular velocity
    v2_circ = np.sqrt(mu / r2)

    # Combined burn with plane change (law of cosines on velocity triangle)
    dv2 = np.sqrt(v_a**2 + v2_circ**2 - 2 * v_a * v2_circ * np.cos(delta_i))

    return {
        'dv1': dv1,
        'dv2': dv2,
        'dv_total': dv1 + dv2,
        'va': v_a,
        'vcirc2': v2_circ
    }

# =====================================================
# Plane Change Maneuvers
# =====================================================

def dv_inclination_change(v1, v2, di):
    """
    Delta-V required for an inclination change.

    Args:
        v1 (float): Initial orbital velocity at maneuver point (km/s).
        v2 (float): Final orbital velocity at maneuver point (km/s).
        di (float): Inclination change (rad).

    Returns:
        float: Delta-V (km/s).
    """
    return np.sqrt(v1**2 + v2**2 - (2 * v1 * v2 * np.cos(di)))


# =====================================================
# Rendezvous and Phasing
# =====================================================

def dtheta_chaser(ac, at):
    """
    Phase angle change per orbit between chaser and target.

    Args:
        ac (float): Chaser semi-major axis (km).
        at (float): Target semi-major axis (km).

    Returns:
        float: Phase angle change (rad).
    """
    return 2 * np.pi * (1 - (ac / at) ** 1.5)


def arc_length_from_phase_change(ac, at):
    """
    Arc length corresponding to phase change.

    Args:
        ac (float): Chaser semi-major axis (km).
        at (float): Target semi-major axis (km).

    Returns:
        float: Arc length (km).
    """
    return 2 * np.pi * at * (1 - (ac / at) ** 1.5)


def phase_angle_rendezvous_begin(ac, at):
    """
    Initial phase angle required to begin rendezvous.

    Args:
        ac (float): Chaser semi-major axis (km).
        at (float): Target semi-major axis (km).

    Returns:
        float: Required phase angle (rad).
    """
    return np.pi * (1 - ((ac + at) / (2 * at)) ** 1.5)


def wait_time(thetaf, thetai, wt, wc, k):
    """
    Waiting time until rendezvous.

    Args:
        thetaf (float): Final phase angle (rad).
        thetai (float): Initial phase angle (rad).
        wt (float): Target angular velocity (rad/s).
        wc (float): Chaser angular velocity (rad/s).
        k (int): Number of extra revolutions.

    Returns:
        float: Waiting time (s).
    """
    return (thetaf - thetai) / (wt - wc) + (2 * np.pi * k) / (wt - wc)


def time_of_flight_rendezvous(ac, at, u):
    """
    Time of flight for a rendezvous transfer.

    Args:
        ac (float): Chaser semi-major axis (km).
        at (float): Target semi-major axis (km).
        u (float): Gravitational parameter μ (km^3/s^2).

    Returns:
        float: Transfer time (s).
    """
    atrans = (ac + at) / 2
    return np.pi * np.sqrt(atrans ** 3 / u)


def ac_coplanar(u, thetai, wt):
    """
    Chaser semi-major axis for coplanar rendezvous.

    Args:
        u (float): Gravitational parameter μ (km^3/s^2).
        thetai (float): Initial phase angle (rad).
        wt (float): Target angular velocity (rad/s).

    Returns:
        float: Chaser semi-major axis (km).
    """
    return (u * ((thetai / (2 * np.pi * wt)) ** 2)) ** (2 / 3)


def dv_initial_coplanaar(u, at, ac):
    """
    Initial delta-V for coplanar rendezvous.

    Args:
        u (float): Gravitational parameter μ (km^3/s^2).
        at (float): Target semi-major axis (km).
        ac (float): Chaser semi-major axis (km).

    Returns:
        float: Delta-V (km/s).
    """
    return np.sqrt(u * (2 / at - 1 / ac)) - np.sqrt(u / at)


def dv_total_coplanar(u, ac, at):
    """
    Total delta-V for coplanar rendezvous.

    Args:
        u (float): Gravitational parameter μ (km^3/s^2).
        ac (float): Chaser semi-major axis (km).
        at (float): Target semi-major axis (km).

    Returns:
        float: Total delta-V (km/s).
    """
    return 2 * np.abs(np.sqrt(u * (2 / at - 1 / ac)) - np.sqrt(u / at))


def parse_tle(line1: str, line2: str):
    """
    Parse TLE lines that are separated by spaces rather than fixed-width columns.
    Returns orbital elements and derived semi-major axis.
    """

    # Split the lines by any whitespace
    parts1 = line1.strip().split()
    parts2 = line2.strip().split()

    # Example: line1 fields
    # 0: '1', 1: '25544U', 2: '98067A', 3: '21275.51167824', ...
    satnum = parts1[1]  # e.g., '25544U'

    # Epoch field in YYDDD.DDDDDDDD format
    epoch_field = parts1[3]
    year = int(epoch_field[0:2])
    year += 2000 if year < 57 else 1900
    day_of_year = float(epoch_field[2:])
    epoch = datetime(year, 1, 1) + timedelta(days=day_of_year - 1)

    # Example: line2 fields
    # 0:'2', 1:'25544', 2:'51.6442', 3:'247.4627', 4:'0006719',
    # 5:'130.5360', 6:'325.0288', 7:'15.48954326', 8:'99999'
    inclination = float(parts2[2])      # degrees
    raan = float(parts2[3])             # degrees
    eccentricity = float(f"0.{parts2[4]}")
    arg_perigee = float(parts2[5])      # degrees
    mean_anomaly = float(parts2[6])     # degrees
    mean_motion = float(parts2[7])      # rev/day
    rev_number = int(parts2[8])

    # Compute semi-major axis (from mean motion)
    mu_earth = 398600.4418  # km^3/s^2
    n_rad_s = mean_motion * 2 * np.pi / 86400.0
    semi_major_axis = (mu_earth / (n_rad_s ** 2)) ** (1/3)

    return {
        "satnum": satnum,
        "epoch": epoch,
        "inclination_deg": inclination,
        "raan_deg": raan,
        "eccentricity": eccentricity,
        "arg_perigee_deg": arg_perigee,
        "mean_anomaly_deg": mean_anomaly,
        "mean_motion_rev_per_day": mean_motion,
        "rev_number": rev_number,
        "semi_major_axis_km": semi_major_axis
    }



def nodal_regression_rate(a, e, i, J2=1.08263e-3, Re=6378.137, mu=398600.4418):
    """
    Computes the rate of nodal regression (RAAN precession) due to Earth's oblateness (J2 effect).

    Parameters
    ----------
    a : float
        Semi-major axis [km]
    e : float
        Eccentricity
    i : float
        Inclination [radians]
    J2 : float, optional
        Earth's J2 coefficient (default: 1.08263e-3)
    Re : float, optional
        Earth's mean equatorial radius [km] (default: 6378.137)
    mu : float, optional
        Earth's gravitational parameter [km^3/s^2] (default: 398600.4418)

    Returns
    -------
    float
        Rate of change of RAAN (radians per second)
    """
    n = np.sqrt(mu / a ** 3)
    dOmega_dt = -1.5 * n * (Re / a) ** 2 * J2 * np.cos(i) / (1 - e ** 2) ** 2
    return dOmega_dt


def sun_sync_inclination(a, e=0, retrograde=True, J2=1.08263e-3, Re=6378.137, mu=398600.4418):
    """
    Computes the inclination for a Sun-synchronous orbit given semi-major axis.

    Parameters
    ----------
    a : float
        Semi-major axis [km]
    e : float, optional
        Eccentricity (default: 0)
    retrograde : bool, optional
        If True, returns retrograde inclination (~98°). If False, prograde (~82°)
    J2, Re, mu : float, optional
        Earth constants

    Returns
    -------
    float
        Inclination [radians]
    """
    target_rate = 1.991e-7  # rad/s (≈ 360° per year)
    if not retrograde:
        target_rate *= -1  # opposite direction for prograde

    def f(i):
        return nodal_regression_rate(a, e, i, J2, Re, mu) - target_rate

    i_guess = np.deg2rad(98 if retrograde else 82)
    i_sol = fsolve(f, i_guess)[0]
    return i_sol


def sun_sync_semi_major_axis(incl_deg, e=0.0):
    """
    Compute the semi-major axis (km) of a Sun-synchronous orbit
    for a given inclination (degrees) and eccentricity.

    Args:
        incl_deg (float): Inclination in degrees
        e (float): Orbital eccentricity (default=0)

    Returns:
        float: Semi-major axis in kilometers
    """
    # Constants
    mu = 398600.4418  # km^3/s^2 (Earth)
    J2 = 1.08262668e-3
    Re = 6378.137  # km (Earth radius)
    T_year = 365.2422 * 24 * 3600  # s

    # Desired precession rate (rad/s)
    omega_dot_target = 2 * np.pi / T_year

    # Inclination in radians
    i = np.deg2rad(incl_deg)

    # Solve for mean motion n using J2 relation
    # omega_dot = -(3/2)*J2*(Re/a)^2 * n * cos(i)/(1 - e^2)^2
    # rearrange -> n = omega_dot * (2/3) * (a/Re)^2 * (1 - e^2)^2 / (J2 * cos(i))
    # but n also = sqrt(mu/a^3)
    # solve iteratively for a

    def f(a):
        n = np.sqrt(mu / a ** 3)
        return abs((3 / 2) * J2 * (Re ** 2) * n * np.cos(i) / (a ** 2 * (1 - e ** 2) ** 2)) - omega_dot_target

    # Simple numerical solve by sweeping a reasonable range (6500–8000 km)
    a_vals = np.linspace(6500, 8000, 10000)
    f_vals = np.abs([f(a) for a in a_vals])
    a_best = a_vals[np.argmin(f_vals)]
    return a_best
