import numpy as np


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
    return 2 * np.pi / np.sqrt(u) * (a ** 1.5)


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


def total_hohman_dv(u, r1, r2, at):
    """
    Total delta-V for a Hohmann transfer.

    Args:
        u (float): Gravitational parameter μ (km^3/s^2).
        r1 (float): Initial orbit radius (km).
        r2 (float): Final orbit radius (km).
        at (float): Semi-major axis of transfer orbit (km).

    Returns:
        float: Total delta-V (km/s).
    """
    return (np.sqrt(u / r2) - np.sqrt(u * (2 / r2 - 1 / at)) +
            np.sqrt(u * (2 / r1 - 1 / at)) - np.sqrt(u / r1))


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
