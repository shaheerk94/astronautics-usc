def oxygen_fluence(n0, V0, T):
    """
    
    Args:
        n0: oxygen number density, cm^-3
        V0: spacecraft initial velocity, m/s
        T: time interval, s

    Returns:
        o_fluence: oxygen fluence
    """
    V_cm2 = V0 * 100 # cm^2/s
    return n0 * V_cm2 * T
