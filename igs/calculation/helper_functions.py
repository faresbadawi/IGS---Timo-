# helper_functions.py
import numpy as np

def generiere_knotenvektor(p, n_elements):
    # Berechnung des Knotenabstands basierend auf der Anzahl der Elemente
    knotenabstand = 1.0 / n_elements

    # Generieren des Knotenvektors im Bereich [0, 1] mit 3 Dezimalstellen
    knotenvektor = [round(i * knotenabstand, 3) for i in range(n_elements + 1)]

    # Wiederholen des Anfangs (0) p + 1 Mal und des Endes (1) p + 1 Mal
    knotenvektor = [0.0] * (p) + knotenvektor + [1.0] * (p)

    # Formatieren des Knotenvektors mit 3 Dezimalstellen und eckigen Klammern
    formatierter_knotenvektor = "[" + ", ".join([format(e, '.3f') for e in knotenvektor]) + "]"
    print("Knotenvektor:", formatierter_knotenvektor)
    return knotenvektor

def integral_bereiche_festlegen(n_elements):
    schritt = 1 / n_elements
    bereiche = []

    for i in range(n_elements):
        a = i * schritt
        b = (i + 1) * schritt
        bereich_name = f"Element {i+1}"  # Benenne den Bereich als "Element 1", "Element 2", usw.
        bereiche.append((bereich_name, (a, b)))
    return bereiche