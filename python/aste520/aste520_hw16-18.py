import numpy
import numpy as np

import globals as g


def problem_16():
    print('## Homework 16 ##')
    d0 = 8.0 * g.ft_to_m # m
    m0 = 1600 # kg
    w0 = 24 * 2 * np.pi / 60 # rad/s
    F0 = 2.5 * g.lbf_to_N # N

    # Part A
    print('# Part A')
    Iz = m0 * d0**2 /8
    print("Iz", Iz)

    # Part B
    alpha = 3 * F0 * d0/ (2 * Iz)
    print('Angular Acceleration', alpha, 'rad/s')
    print('Angular Acceleration', alpha*180/np.pi, 'deg/s')

    # Part C
    t0 = w0/alpha
    print('Time', t0, 'seconds')

    # Part D
    t = 20 # s
    w1 = alpha * t
    print('w1', w1, 'rad/s')
    print('w1', w1*180/np.pi, 'deg/s')
    print('w1', w1 * 60/(2*np.pi), 'rpm')

    theta1 = 1/2 * alpha * t**2
    print('theta1', theta1, 'radians')
    print('theta1', theta1 * 180/np.pi, 'degrees')
    print('theta1', theta1 / (2*np.pi), 'revolutions')

    # Part E
    Hp = Iz * w0
    print('Angular Momentum', Hp, 'kg m^2/s')
    E = (1/2) * Iz * w0 **2
    print('Energy', E, 'Joules')

    # Part F
    T1 = 5e-4 # Nm
    omega = T1/Hp
    print('Precession', omega, 'rad/s')
    print('Precession', omega*180/np.pi, 'deg/s')
    print('Precession', omega*180/np.pi * 3600 * 24, 'deg/day')






def main():
    problem_16()

if __name__ == "__main__":
    main()