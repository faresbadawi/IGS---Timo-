# gauss_integration.py
import numpy as np

def gauss_punkt_integrationsdaten(a, b, p):
    num_punkte = p + 1
    gauss_punkte, gewichte = np.polynomial.legendre.leggauss(num_punkte)
    abgebildete_punkte = 0.5 * (b - a) * gauss_punkte + 0.5 * (a + b)
    return gauss_punkte, gewichte, abgebildete_punkte

def det_Jξ(a, b):
    return (b - a) / 2
