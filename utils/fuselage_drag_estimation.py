
import numpy as np

def fuselage_drag(M, T, rho, a, l, d, n=10):
    """
    Estimate fuselage drag based on assumptions

    ASSUMPTIONS:
    -negligible interference effects --> Q=1
    -Von Karman Ogive body


    :param h: cruise height (m)
    :param M: cruise Mach
    :param l: fuselage length (m)
    :param d: fuselage width (m)
    :param F: Form factor
    :param n: number of fuselage discretisation segments
    :return:
    """
    def vonKarmanOgive(l, d, n):
        """
        Returns the radii of two sequential and mirror von karman ogive bodies.
        :param l: Length of the total body
        :param d: Maximum diameter
        :param n: Number of lengthwise segments
        :return:
        """
        x = np.linspace(0, l, n+1)

        nose = x[:int(np.ceil((n+1) / 2))]

        x_centered = nose - l/4
        x_centered_nondim = 4*x_centered / l

        r_x_over_r_0_sqrd = 1 / np.pi * (
                    x_centered_nondim * np.sqrt(1 - x_centered_nondim ** 2) + np.arccos(-x_centered_nondim))
        d_x = d * np.sqrt(r_x_over_r_0_sqrd)

        VKO = np.zeros(n+1)
        VKO[:int(np.ceil((n+1) / 2))] = d_x
        VKO[int(np.floor((n+1) / 2)):] = np.flip(d_x)
        return x, VKO

    feettom = 0.3048

    V = a*M
    #V = 1925.7*feettom
    q = 1/2*rho*V**2

    mu = 1.458e-6*T**1.5/(T+110.4)
    v =  mu/rho

    #von karman ogive
    x, d_x = vonKarmanOgive(l,d,n)

    p_x = np.pi*d_x

    S_x = p_x*(l/n)
    surface_area = sum(S_x)

    A_max = np.pi*d**2/4

    Re_x = V*x[1:]/v

    C_f = 0.455/(np.log10(Re_x)**2.58*(1+0.144*M**2)**0.65)

    friction_drag = sum(C_f*q*S_x[1:])

    CD_w = 4*A_max/(np.pi*l**2)

    CD_w = 4*(d/l)**2


    fineness_factor_inv = l/d
    form_factor = 1 + 60/fineness_factor_inv**3 + fineness_factor_inv/400

    wave_drag = q*A_max*CD_w

    CD_0S = (friction_drag + wave_drag)/q
    #To be decided by the main wing surface area
    return CD_0S

if __name__ == '__main__':
    from ISACalculator import ISACalculator
    feettom = 0.3048

    h = 55000 #ft
    M = 2.1 #m
    l = 126*feettom #m
    d = 9*feettom #m

    T, p, rho, a = ISACalculator(h * feettom)

    fuselage_drag(M,l,d, T, rho, a)
