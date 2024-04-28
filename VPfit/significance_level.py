import numpy as np

def cos_siglevel(W, wavelength, b, snpix=None, snx=None, disp=None, bin=None, xopt=None, fcx=None):
    if snpix is None and snx is None:
        print('SIGLEVEL: either SNPIX or SNX must be specified!!!')
        return -1

    if disp is None:
        disp = 9.97 if wavelength <= 1425 else 12.23

    if bin is None:
        bin = 1

    dx = 1e3 * b * wavelength / (2.998e5 * disp)

    xopt = 1.605 * dx + 5.1 * dx ** (-0.25)
    eta = 0.15 + xopt ** 0.37
    fcx = 0.743 - 0.185 * np.exp(-dx / 11.6)

    sn1 = snpix if bin == 1 else snpix / (0.15 + bin ** 0.37)

    siglevel = snx * (W / disp) * fcx / xopt if snx is not None else sn1 * (W / disp) * eta * fcx / xopt

    xopt /= bin

    return siglevel
