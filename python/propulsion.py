import math
from typing import Optional

def _tau(M: float, gamma: float) -> float:
    """Helper: 1 + (γ−1)/2 * M²"""
    return 1.0 + 0.5 * (gamma - 1.0) * M**2


def T_static(M: float, gamma: float, T0: float) -> float:
    """Static temperature from stagnation temperature."""
    return T0 / _tau(M, gamma)


def P_static(M: float, gamma: float, P0: float) -> float:
    """Static pressure from stagnation pressure."""
    return P0 / (_tau(M, gamma) ** (gamma / (gamma - 1.0)))


def rho_static(
    M: float,
    gamma: float,
    *,
    rho0: Optional[float] = None,
    P0: Optional[float] = None,
    T0: Optional[float] = None,
    R: Optional[float] = None,
) -> float:
    """
    Static density. Two paths
      1) provide rho0
      2) provide P0, T0, R
    """
    if rho0 is not None:
        return rho0 / (_tau(M, gamma) ** (1.0 / (gamma - 1.0)))

    if P0 is not None and T0 is not None and R is not None:
        P = P_static(M, gamma, P0)
        T = T_static(M, gamma, T0)
        return P / (R * T)

    raise ValueError("Provide rho0 or P0, T0, and R.")


def a_from_T(gamma: float, R: float, T: float) -> float:
    """Speed of sound from static temperature."""
    return math.sqrt(gamma * R * T)


def a_static(M: float, gamma: float, T0: float, R: float) -> float:
    """Speed of sound using stagnation temperature."""
    return a_from_T(gamma, R, T_static(M, gamma, T0))
