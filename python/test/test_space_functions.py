import numpy as np
import pytest
from pint import UnitRegistry
from space_functions import apparent_magnitude, planet_distance, synodic_period

ureg = UnitRegistry()

def test_apparent_magnitude_at_10pc():
    # At 10 parsecs, apparent magnitude = absolute magnitude
    r = 10 * ureg.parsec
    abs_mag = 5.0
    m = apparent_magnitude(r, abs_mag)
    assert np.isclose(m, abs_mag, atol=1e-8)

def test_apparent_magnitude_closer():
    # At 1 parsec, apparent magnitude should be brighter (smaller number)
    r = 1 * ureg.parsec
    abs_mag = 5.0
    m = apparent_magnitude(r, abs_mag)
    assert np.isclose(m, 0.0, atol=1e-8)

def test_apparent_magnitude_farther():
    # At 100 parsecs, apparent magnitude = abs_mag + 5
    r = 100 * ureg.parsec
    abs_mag = 2.0
    m = apparent_magnitude(r, abs_mag)
    assert np.isclose(m, 7.0, atol=1e-8)

def test_planet_distance_same_point():
    # Identical coordinates should give distance = 0
    d_au, d_km = planet_distance(0, 0, 1, 0, 0, 1)
    assert np.isclose(d_au, 0.0, atol=1e-10)
    assert np.isclose(d_km, 0.0, atol=1e-2)

def test_planet_distance_opposite():
    # Two planets opposite each other at 1 AU → 2 AU apart
    d_au, d_km = planet_distance(0, 0, 1, 180, 0, 1)
    assert np.isclose(d_au, 2.0, atol=1e-8)
    assert np.isclose(d_km, 2.0 * 149597870, atol=1e-2)

def test_synodic_period_simple_case():
    # Example from solar rotation: 24.5-day sidereal → ~26.3-day synodic
    P_syn = synodic_period(24.5)
    assert np.isclose(P_syn, 26.26, rtol=1e-3)

def test_synodic_period_equal_orbits():
    # If sidereal = Earth's orbital period, synodic should be infinite (denominator ~0)
    with pytest.raises(ZeroDivisionError):
        synodic_period(365.25)
