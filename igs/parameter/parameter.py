# Material and calculation parameters.py
import numpy as np


# Material- und Berechnungsparameter
def berechne_steuergroessen():
    # Berechnungsparameter
    p = 3  # Grad der Ansatzfunktionen
    n_elements = 2  # Anzahl der Elemente im Modell

    # Materialparameter
    kontrollpunkte_vektor = [0, 0.16667, 0.5, 0.83333, 1]  # X-Koordinaten der Kontrollpunkte des Balkens (in Metern)
    E = 1e7  # Elastizitätsmodul des Materials (in Kilonewton pro Quadratmeter, kN/m²)
    b = 0.1  # Breite des Balkens (in Metern)
    h = 0.01  # Höhe des Balkens (in Metern)
    nu = 0.2  # Poisson-Verhältnis des Materials (dimensionslos)
    alpha = 5 / 6  # Ein Materialparameter (dimensionslos)
    k_b = (E * b * h ** 3) / 12  # Biegesteifigkeit des Balkens (in kN·m²)
    k_s = alpha * (E * b * h) / (2 * (1 + nu))  # Schubsteifigkeit des Balkens (in kN/m)

    return p, n_elements, kontrollpunkte_vektor, k_b, k_s


# Auflagerbedingungen
def erstelle_auflagerbedingungen():
    auflagerbedingungen = [1, 0, 0, 0, 1]  # Beispiel: 1 für festes Lager, 0 für freies Lager
    return auflagerbedingungen
