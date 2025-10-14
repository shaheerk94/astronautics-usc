import pint
import numpy as np
import globals as g

ureg = pint.UnitRegistry()
Distance = pint.Quantity


def apparent_magnitude(r: Distance, abs_magnitude: float) -> float:
    """
    Calculate the apparent magnitude of an object.

    Args:
        r (Quantity): Distance from the object to the observer. Must have units of length (e.g., parsecs).
        abs_magnitude (float): Absolute magnitude of the object.

    Returns:
        float: Apparent magnitude.
    """
    r_pc = r.to(ureg.parsec).magnitude  # convert to parsecs
    return abs_magnitude - 5 + 5 * np.log10(r_pc)


def absolute_magnitude(r: Distance, apparent_magnitude: float) -> float:
    """
    Calculate the absolute magnitude of an object given apparent magnitude and distance in km.
    :param r:  distance in km.
    :param apparent_magnitude: magnitude from observer.
    :return: absolute magnitude.
    """
    r_pc = r.to(ureg.parsec).magnitude
    return apparent_magnitude + 5 - 5 * np.log10(r_pc)


def synodic_period(P_sidereal_days, P_earth_days=365.25):
    """
    Compute the synodic (apparent) rotation period as seen from Earth.

    Parameters
    ----------
    P_sidereal_days : float
        Sidereal rotation period (days)
    P_earth_days : float, optional
        Earth's orbital period (days), default = 365.25

    Returns
    -------
    P_synodic_days : float
        Synodic rotation period (days)
    """
    return 1.0 / (1.0 / P_sidereal_days - 1.0 / P_earth_days)


def dms_to_decimal(degrees: int, minutes: int = 0, seconds: float = 0.0, hemisphere: str = None) -> float:
    """
    Convert degree-minute-second (DMS) notation to decimal degrees.

    Args:
        degrees (int): Degrees part of the coordinate.
        minutes (int, optional): Minutes part of the coordinate. Defaults to 0.
        seconds (float, optional): Seconds part of the coordinate. Defaults to 0.0.
        hemisphere (str, optional): Hemisphere indicator ('N', 'S', 'E', 'W').
                                    If given, applies correct sign. Defaults to None.

    Returns:
        float: Decimal degrees.
    """
    decimal = abs(degrees) + minutes / 60 + seconds / 3600

    if hemisphere:
        hemisphere = hemisphere.upper()
        if hemisphere in ['S', 'W']:
            decimal = -decimal

    return decimal


def planet_distance(lam1_deg, beta1_deg, r1, lam2_deg, beta2_deg, r2):
    """
    Compute the heliocentric distance between two planets given in
    ecliptic coordinates using the spherical law of cosines.

    Parameters
    ----------
    lam1_deg, lam2_deg : float
        Ecliptic longitudes (degrees)
    beta1_deg, beta2_deg : float
        Ecliptic latitudes (degrees)
    r1, r2 : float
        Heliocentric distances (AU)

    Returns
    -------
    S_AU : float
        Distance between the two planets (AU)
    S_km : float
        Distance between the two planets (km)
    """

    # Convert to radians
    lam1 = np.deg2rad(lam1_deg)
    beta1 = np.deg2rad(beta1_deg)
    lam2 = np.deg2rad(lam2_deg)
    beta2 = np.deg2rad(beta2_deg)

    # Angular separation
    cos_psi = planet_angles(lam1_deg, beta1_deg, lam2_deg, beta2_deg)

    # Straight-line distance
    S_AU = np.sqrt(r1**2 + r2**2 - 2 * r1 * r2 * cos_psi)

    # Convert AU → km
    AU_to_km = 149597870
    S_km = S_AU * AU_to_km

    return S_AU, S_km

# --------------------------------------------------------------------
# Ecliptic → Celestial (RA/Dec)
# --------------------------------------------------------------------
def ecliptic_to_celestial(lambda_deg, beta_deg):
    lam = np.deg2rad(lambda_deg)
    bet = np.deg2rad(beta_deg)
    # ecliptic unit vector
    x = np.cos(lam) * np.cos(bet)
    y = np.sin(lam) * np.cos(bet)
    z = np.sin(bet)
    # rotate about x by +ε to equatorial
    y_eq = y*np.cos(g.obliquity_earth) - z*np.sin(g.obliquity_earth)
    z_eq = y*np.sin(g.obliquity_earth) + z*np.cos(g.obliquity_earth)
    x_eq = x
    # recover RA, Dec
    alpha = np.arctan2(y_eq, x_eq) % (2*np.pi)
    delta = np.arcsin(z_eq)
    return np.rad2deg(alpha), np.rad2deg(delta)

# --------------------------------------------------------------------
# Celestial → Ecliptic
# --------------------------------------------------------------------
def celestial_to_ecliptic(alpha_deg, delta_deg):
    alp = np.deg2rad(alpha_deg)
    dec = np.deg2rad(delta_deg)
    # equatorial unit vector
    x = np.cos(alp) * np.cos(dec)
    y = np.sin(alp) * np.cos(dec)
    z = np.sin(dec)
    # rotate about x by -ε to ecliptic
    y_ec = y*np.cos(g.obliquity_earth) + z*np.sin(g.obliquity_earth)
    z_ec = -y*np.sin(g.obliquity_earth) + z*np.cos(g.obliquity_earth)
    x_ec = x
    # recover λ, β
    lam = np.arctan2(y_ec, x_ec) % (2*np.pi)
    bet = np.arcsin(z_ec)
    return np.rad2deg(lam), np.rad2deg(bet)
# --------------------------------------------------------------------
# Ecliptic spherical → XYZ Cartesian
# --------------------------------------------------------------------
def ecliptic_to_xyz(r, lambda_deg, beta_deg):
    """
    Convert ecliptic spherical coords to heliocentric Cartesian (AU).

    Parameters
    ----------
    r : float
        Radius (AU)
    lambda_deg : float
        Ecliptic longitude [deg]
    beta_deg : float
        Ecliptic latitude [deg]

    Returns
    -------
    X, Y, Z : float
        Cartesian coords (same length units as r)
    """
    lam = np.deg2rad(lambda_deg)
    beta = np.deg2rad(beta_deg)
    X = r * np.cos(lam) * np.cos(beta)
    Y = r * np.sin(lam) * np.cos(beta)
    Z = r * np.sin(beta)
    return X, Y, Z

# --------------------------------------------------------------------
# Celestial spherical → XYZ Cartesian
# --------------------------------------------------------------------
def celestial_to_xyz(r, alpha_deg, delta_deg):
    """
    Convert celestial (RA, Dec) to Cartesian coords.

    Parameters
    ----------
    r : float
        Radius (AU)
    alpha_deg : float
        Right ascension [deg]
    delta_deg : float
        Declination [deg]

    Returns
    -------
    X, Y, Z : float
        Cartesian coords
    """
    alpha = np.deg2rad(alpha_deg)
    delta = np.deg2rad(delta_deg)
    X = r * np.cos(alpha) * np.cos(delta)
    Y = r * np.sin(alpha) * np.cos(delta)
    Z = r * np.sin(delta)
    return X, Y, Z


# --------------------------------------------------------------------
# Angle between Stars, given ecliptic coords
# --------------------------------------------------------------------
def planet_angle(lam1_deg, beta1_deg, lam2_deg, beta2_deg):
    """
    Calculate angles between celestial objects given ecliptic coordinates.

    Parameters
    ----------
    lam1_deg : float
        ecliptic longitude of object 1 [deg]
    beta1_deg : float
        ecliptic latitude of object 2 [deg]
    lam2_deg : float
        ecliptic longitude of object 1 [deg]
    beta2_deg : float
        ecliptic latitude of object 2 [deg]

    Returns
    -------
    psi : float
        radians
    """
    # Convert to radians
    lam1 = np.deg2rad(lam1_deg)
    beta1 = np.deg2rad(beta1_deg)
    lam2 = np.deg2rad(lam2_deg)
    beta2 = np.deg2rad(beta2_deg)

    # Angular separation
    psi = np.acos(np.sin(beta1) * np.sin(beta2) +
               np.cos(beta1) * np.cos(beta2) * np.cos(lam1 - lam2))
    return psi