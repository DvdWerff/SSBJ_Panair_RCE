import math
def ISACalculator(h):
    '''
    Determines the temperature, pressure, density and speed of sound at a given altitude.
    Based on the first three layers (for now) of the ISA atmosphere model
    :param h: altitude [m]
    :return: temperature T in K, pressure p in Pa, density rho in kg/m³, speed of sound a in m/s
    '''
    def CtoK(T):
        return T + 273.15

    def atmosphere(h, a, h0, T0, p0):
        hdiff = h - h0
        if a == 0:
            T = T0
            p = atmosphere_constantT(hdiff, T, p0)
        else:
            T, p = atmosphere_gradientT(hdiff, a, T0, p0)

        rho = p / (R * T)
        SoS = math.sqrt(gamma * R * T)  # speed of sound
        return T, p, rho, SoS

    def atmosphere_constantT(hdiff, T, p0):
        p = p0 * math.exp(-g0 / (R * T) * hdiff)
        return p

    def atmosphere_gradientT(hdiff, a, T0, p0):
        T = T0 + hdiff * a
        p = p0 * math.pow(T / T0, -(g0 / (R * a)))
        return T, p

    # Constants
    g0 = 9.80665
    R = 287.0
    gamma = 1.4

    # Base values
    h0, a0, T0, p0, rho0 = 0, -0.0065, 288.15, 101325.9, 1.225
    h1, a1, T1, p1, rho1 = 11000, 0, CtoK(-56.5), 22632, 0.3639
    h2, a2, T2, p2, rho2 = 20000, 0.0010, CtoK(-56.5), 5474.9, 0.088
    h3 = 32000

    # Calculations
    if 0 <= h < h1:
        T, p, rho, a = atmosphere(h, a0, h0, T0, p0)
        return T, p, rho, a

    if h1 < h < h2:
        T, p, rho, a = atmosphere(h, a1, h1, T1, p1)
        return T, p, rho, a

    if h2 < h < h3:
        T, p, rho, a = atmosphere(h, a2, h2, T2, p2)
        return T, p, rho, a

    else:
        raise ValueError('Requested altitude is currently not implemented.')
