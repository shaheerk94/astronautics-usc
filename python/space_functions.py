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
