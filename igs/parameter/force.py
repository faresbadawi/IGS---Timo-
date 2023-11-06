import numpy as np

from calculation.helper_functions import generiere_knotenvektor, integral_bereiche_festlegen
from calculation.gauss_integration import gauss_punkt_integrationsdaten, det_Jξ
from calculation.basis_functions import basisfunktion
from calculation.matrix_operations import determinante_und_inverse_berechnen

# Konstanten
q_0 = 3.0  # Einheit in KN/m
L = 6.0  # Einheit in m
n_elements = 2
p = 3

# Kontrollpunkte
kontrollpunkte_vektor = [0, 0.1666667, 0.5, 0.83333333, 1]

# Setze den Knotenvektor und p
knotenvektor = generiere_knotenvektor(p, n_elements)

# Festlegung der Integralbereiche
integral_ranges = integral_bereiche_festlegen(n_elements)

# Initialisiere die Gesamtsummenmatrix
gesamt_result_matrix_bending_T = np.zeros((1, len(kontrollpunkte_vektor))).reshape(-1, 1)

# Verwende die integral_ranges in gauss_point_integration_data und dann in ableitung_kter_ordnung
for bereich_name, (a, b) in integral_ranges:
    gauss_points, weights, mapped_points = gauss_punkt_integrationsdaten(a, b, p)
    basisfunktion_bending = np.zeros((1, len(kontrollpunkte_vektor)))  # Initialisiere die Gesamtsumme-Matrix

    det_Jx_values = []  # Liste zur Speicherung der det_Jx-Werte für jedes t
    inverse_Jx_values = []  # Liste zur Speicherung der inverse_Jx-Werte für jedes t

    for i, (t, weight) in enumerate(zip(mapped_points, weights)):
        # Einsetzen der t-Werte und Umwandeln in eine Matrix (bending part)
        basisfunktion_array_bending = [basisfunktion(t, knotenvektor, i, p) for i in range(len(knotenvektor) - p - 1)]
        basisfunktion_matrix_bending = np.array(basisfunktion_array_bending).reshape(1, -1)
        basisfunktion_matrix_bending_T = np.array(basisfunktion_array_bending).reshape(-1, 1)

        # Berechne det_Jx und inverse_Jx für das aktuelle t
        det_ξ_value = det_Jξ(a, b)

        # Berechnung von det_Jx und inverse_Jx
        det_Jx, inverse_Jx = determinante_und_inverse_berechnen(basisfunktion_matrix_bending, kontrollpunkte_vektor)

        # Multipliziere die ableitung_kter_ordnung_matrix mit dem Gewicht, det_ξ und inverse_Jx
        result_matrix_bending_T = basisfunktion_matrix_bending_T * weight * det_ξ_value * det_Jx

        # Addiere die Basisfunktionsmatrix für diesen Gauss-Punkt zur Gesamtsummenmatrix
        gesamt_result_matrix_bending_T += result_matrix_bending_T

        # Speichern der Werte
        det_Jx_values.append(det_Jx)
        inverse_Jx_values.append(inverse_Jx)

        # Ausgabe für die Kontrolle
        print(f"\nBereich ({a}, {b}), Gauss-Punkt {i + 1} (t = {t}):")
        print("basisfunktion_Matrix_bending:")
        print(basisfunktion_matrix_bending_T)
        print("Gewicht (Weight) =", weight)
        print("Gauss-Punkt =", gauss_points[i])
        print("Abgebildeter Punkt =", t)
        print("det_Jx =", det_Jx)
        print("inverse_Jx =", inverse_Jx)

# Gesamtsummenmatrix anzeigen
print("\nGesamtsummenmatrix der Basisfunktionen:")
print(gesamt_result_matrix_bending_T)

#Gesamtsummenmatrix Speichern
def berechne_kraftvektor():
    kraftvektor = gesamt_result_matrix_bending_T.reshape(1, -1)

    return kraftvektor
