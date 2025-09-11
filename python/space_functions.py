import pint
import numpy as np

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
    cos_psi = (np.sin(beta1) * np.sin(beta2) +
               np.cos(beta1) * np.cos(beta2) * np.cos(lam1 - lam2))

    # Straight-line distance
    S_AU = np.sqrt(r1**2 + r2**2 - 2 * r1 * r2 * cos_psi)

    # Convert AU → km
    AU_to_km = 149597870
    S_km = S_AU * AU_to_km

    return S_AU, S_km

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


