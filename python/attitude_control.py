import numpy as np
from scipy.constants import c

import numpy as np


def solar_radiation_force(th, CA, CS, CD, Fs, A):
    """
    Radiation-pressure force on a flat plate for arbitrary optical fractions.

    Args:
        theta (float): Incidence angle between beam and surface normal.
                       0° = normal incidence; 90° = grazing. (radians)
        CA (float): Absorbed fraction.
        CS (float): Specularly reflected fraction.
        CD (float): Diffusely (Lambertian) reflected fraction.
        Fs (float): Solar flux [W/m^2].
        A  (float): Plate area [m^2].

    Returns:
        (fx, fy): Tuple of force components [N].
                  fx is along +X (surface normal into plate),
                  fy is along +Y (projection of beam in-plane).
    """
    # Optional: re-normalize tiny sum errors
    s = CA + CS + CD
    if not np.isclose(s, 1.0):
        CA, CS, CD = CA/s, CS/s, CD/s

    mu = np.cos(th)
    if mu <= 0:  # back side or grazing beyond 90° -> no illumination
        return 0.0, 0.0

    sin_th = np.sqrt(max(0.0, 1.0 - mu**2))
    P0 = Fs / 299792458.0  # N/m^2 (radiation pressure scale)
    K = P0 * A

    # Normal (+X) and tangential (+Y) components
    fx = K * (mu**2) * (CA + 2.0*CS + (2.0/3.0)*CD)
    fy = K * (mu * sin_th) * (CA + CD)

    return fx, fy

