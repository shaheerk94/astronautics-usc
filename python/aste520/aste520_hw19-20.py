import globals as g

def problem_19():

    print("### Problem 19 ###")
    print("### Part 1 ###")
    ft = 16.5 * 1000 * g.lbf_to_N # Newtons
    isp = 444 # seconds

    ueq = g.g0 * isp
    print('ueq', ueq, 'm/s')
    mdot = ft / (isp * g.g0)
    print('mdot', mdot, 'kg/s')

    print("### Part 2 ###")
    isp = 194 # s
    ft = 2.5 * 3 * g.lbf_to_N # N
    t = 73.5 # s
    mdot = ft / (isp * g.g0)
    print('mdot', mdot, 'kg/s')
    m_total = mdot * t
    print('m_total', m_total, 'kg')

def problem_20():
    pass


def main():
    problem_19()
    problem_20()


if __name__ == '__main__':
    main()