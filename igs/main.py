# Importieren der erforderlichen Module
import numpy as np
from scipy.linalg import solve
from calculation.helper_functions import generiere_knotenvektor, integral_bereiche_festlegen
from calculation.gauss_integration import gauss_punkt_integrationsdaten, det_Jξ
from calculation.basis_functions import ableitungsbasisfunktion, ableitung_kter_ordnung
from parameter.parameter import berechne_steuergroessen, erstelle_auflagerbedingungen
from parameter.force import berechne_kraftvektor

# Berechnung der Steuergrößen und Materialparameter
p, n_elements, kontrollpunkte_vektor, k_b, k_s = berechne_steuergroessen()

# Setze den Knotenvektor und p
knotenvektor = generiere_knotenvektor(p, n_elements)

# Festlegung der Integralbereiche
integral_ranges = integral_bereiche_festlegen(n_elements)

# Erstellen der Auflagerbedingungen
auflagerbedingungen = erstelle_auflagerbedingungen()

# Aufrufen des Kraftvektors

Z = berechne_kraftvektor()

# Erstelle Dictionaries zur Speicherung der Gesamtsummen der Ergebnismatrizen und B-Parameter-Matrizen
gesamtsummen_dict_bending = {}
gesamtsummen_dict_shear = {}
b_parameter_dict_bending = {}
b_parameter_dict_shear = {}

# Verwende die integral_ranges in gauss_punkt_integrationsdaten und dann in ableitung_kter_ordnung
for bereich_name, (a, b) in integral_ranges:
    gauss_points, weights, mapped_points = gauss_punkt_integrationsdaten(a, b, p)
    gesamtsumme_matrix_bending = np.zeros((1, len(kontrollpunkte_vektor)))  # Initialisiere die Gesamtsumme-Matrix
    gesamtsumme_matrix_shear = np.zeros((1, len(kontrollpunkte_vektor)))  # Initialisiere die Gesamtsumme-Matrix

    det_Jx_values = []  # Liste zur Speicherung der det_Jx-Werte für jedes t
    inverse_Jx_values = []  # Liste zur Speicherung der inverse_Jx-Werte für jedes t

    for i, (t, weight) in enumerate(zip(mapped_points, weights)):
        # Einsetzen der t-Werte und Umwandeln in eine Matrix (bending part)
        ableitung_kter_ordnung_array_bending = [ableitung_kter_ordnung(t, knotenvektor, i, p, 1) for i in range(len(knotenvektor) - p - 1)]
        ableitung_kter_ordnung_matrix_bending = np.array(ableitung_kter_ordnung_array_bending).reshape(1, -1)

        # Einsetzen der t-Werte und Umwandeln in eine Matrix (shear part)
        ableitung_kter_ordnung_array_shear = [ableitung_kter_ordnung(t, knotenvektor, i, p, 2) for i in range(len(knotenvektor) - p - 1)]
        ableitung_kter_ordnung_matrix_shear = np.array(ableitung_kter_ordnung_array_shear).reshape(1, -1)

        # Einsetzen der t-Werte und Umwandeln in einer Matrix, um det_Jx und inverse_Jx zu berechnen
        ableitungsbasisfunktion_array = [ableitungsbasisfunktion(t, knotenvektor, i, p) for i in range(len(knotenvektor) - p - 1)]
        ableitungsbasisfunktion_matrix = np.array(ableitungsbasisfunktion_array).reshape(1, -1)

        # Berechne det_Jx und inverse_Jx für das aktuelle t
        det_ξ_value = det_Jξ(a, b)
        det_Jx = np.dot(ableitungsbasisfunktion_matrix, kontrollpunkte_vektor)
        inverse_Jx = 1 / det_Jx

        # Multipliziere die ableitung_kter_ordnung_matrix mit dem Gewicht, det_ξ und inverse_Jx
        result_matrix_bending = ableitung_kter_ordnung_matrix_bending * weight * det_ξ_value * inverse_Jx * det_Jx**0.5 * k_b**0.5
        result_matrix_shear = ableitung_kter_ordnung_matrix_shear * weight * det_ξ_value * inverse_Jx * det_Jx**0.5 * (k_b**2 / k_s)**0.5

        # Ausgabe für die Kontrolle
        print(f"\n{bereich_name}, Gauss-Punkt {i+1} (t = {t}):")
        print("Ableitung_kter_ordnung_Matrix_bending:")
        print(ableitung_kter_ordnung_matrix_bending)
        print("Ableitung_kter_ordnung_Matrix_shear:")
        print(ableitung_kter_ordnung_matrix_shear)
        print("Gewicht (Weight) =", weight)
        print("Gauss-Punkt =", gauss_points[i])
        print("Abgebildeter Punkt =", t)
        print("det_ξ =", det_ξ_value)
        print("det_Jx =", det_Jx)
        print("inverse_Jx =", inverse_Jx)

        # Füge die Ergebnismatrix zur Gesamtsumme hinzu
        gesamtsumme_matrix_bending += result_matrix_bending
        gesamtsumme_matrix_shear += result_matrix_shear

        # Füge die det_Jx- und inverse_Jx-Werte zur entsprechenden Liste hinzu
        det_Jx_values.append(det_Jx)
        inverse_Jx_values.append(inverse_Jx)

    # Speichere die Gesamtsumme der Ergebnismatrix für das aktuelle Element
    gesamtsummen_dict_bending[bereich_name] = gesamtsumme_matrix_bending
    gesamtsummen_dict_shear[bereich_name] = gesamtsumme_matrix_shear

    # Berechne die B-Parameter-Matrix für das aktuelle Element
    def B_Parameter_Matrix(B_matrix, gesamtsumme_matrix):
        num_params = len(B_matrix[0])
        result = np.zeros((num_params, num_params))

        for i in range(num_params):
            for j in range(num_params):
                result[i][j] = B_matrix[0][i] * B_matrix[0][j]

        return result

    b_parameter_matrix_bending = B_Parameter_Matrix(gesamtsumme_matrix_bending, gesamtsumme_matrix_bending)
    b_parameter_dict_bending[bereich_name] = b_parameter_matrix_bending

    b_parameter_matrix_shear = B_Parameter_Matrix(gesamtsumme_matrix_shear, gesamtsumme_matrix_shear)
    b_parameter_dict_shear[bereich_name] = b_parameter_matrix_shear

    # Zeige die det_Jx- und inverse_Jx-Werte für jedes t in diesem Bereich
    print(f"\ndet_Jx-Werte für {bereich_name}:")
    print(det_Jx_values)
    print(f"\ninverse_Jx-Werte für {bereich_name}:")
    print(inverse_Jx_values)

# Zeige die Gesamtsummen der Ergebnismatrizen für jedes Element
for bereich_name, gesamtsumme_matrix in gesamtsummen_dict_bending.items():
    print(f"\nGesamtsumme der Ergebnismatrizen_bending für {bereich_name}:")
    print(gesamtsumme_matrix)
for bereich_name, gesamtsumme_matrix in gesamtsummen_dict_shear.items():
    print(f"\nGesamtsumme der Ergebnismatrizen_shear für {bereich_name}:")
    print(gesamtsumme_matrix)

# Initialisiere die Gesamt-B-Parameter-Matrix für bending und shear
gesamt_b_parameter_matrix_bending = np.zeros((len(kontrollpunkte_vektor), len(kontrollpunkte_vektor)))
gesamt_b_parameter_matrix_shear = np.zeros((len(kontrollpunkte_vektor), len(kontrollpunkte_vektor)))

# Iteriere durch die Elemente und addiere die B-Parameter-Matrizen für bending und shear für jedes Element zur Gesamtmatrix
for bereich_name, b_parameter_matrix_bending in b_parameter_dict_bending.items():
    gesamt_b_parameter_matrix_bending += b_parameter_matrix_bending

for bereich_name, b_parameter_matrix_shear in b_parameter_dict_shear.items():
    gesamt_b_parameter_matrix_shear += b_parameter_matrix_shear

# Zeige die Gesamt-B-Parameter-Matrix für bending
print("Gesamt-B-Parameter-Matrix-bending:")
print(gesamt_b_parameter_matrix_bending)

# Zeige die Gesamt-B-Parameter-Matrix für shear
print("Gesamt-B-Parameter-Matrix-shear:")
print(gesamt_b_parameter_matrix_shear)

# Initialisiere die Gesamt-B-Parameter-Matrix für shear und bending
gesamt_b_parameter_matrix_shear_bending = gesamt_b_parameter_matrix_shear + gesamt_b_parameter_matrix_bending

# Zeige die Gesamt-B-Parameter-Matrix für shear und bending zusammengefasst
print("Gesamt-B-Parameter-Matrix-shear und bending zusammengefasst:")
print(gesamt_b_parameter_matrix_shear_bending)

# Erzeugen Sie eine Kopie der gesamt_b_parameter_matrix_shear_bending
K_Matrix = np.copy(gesamt_b_parameter_matrix_shear_bending)

# Iterieren Sie durch die Auflagerbedingungen und setzen Sie die entsprechenden Zeilen und Spalten auf null
for i, auflager in enumerate(auflagerbedingungen):
    if auflager == 1:
        K_Matrix[i, :] = 0  # Setze die Zeile auf null
        K_Matrix[:, i] = 0  # Setze die Spalte auf null

# Entfernen der Nullzeilen
K_Matrix = K_Matrix[~np.all(K_Matrix == 0, axis=1)]

# Entfernen der Nullspalten
K_Matrix = K_Matrix[:, ~np.all(K_Matrix == 0, axis=0)]

# Zeigen Sie die K-Matrix
print("K-Matrix:")
print(K_Matrix)

# Erstellen Sie einen Vektor für die aufgebrachten Kräfte
f = np.zeros(K_Matrix.shape[0])# Initialisieren Sie f mit Nullen, basierend auf der aktualisierten K-Matrix

# Fügen Sie die aufgebrachten Kräfte hinzu (Beispiel: F1 an Knotenpunkt 2 und F2 an Knotenpunkt 5)
F1 = 1.0 # Beispiel: Wert der ersten aufgebrachten Kraft
F2 = 1.0  # Beispiel: Wert der zweiten aufgebrachten Kraft
F3 = 1.0  # Beispiel: Wert der zweiten aufgebrachten Kraft

f[0,] = F1  # Setzen Sie die aufgebrachte Kraft F1 am entsprechenden Knotenpunkt
f[1,] = F2  # Setzen Sie die aufgebrachte Kraft F2 am entsprechenden Knotenpunkt
f[2,] = F1  # Setzen Sie die aufgebrachte Kraft F3 am entsprechenden Knotenpunkt

# Zeigen Sie den aktualisierten Kraftvektor f
print("Aktualisierter Kraftvektor f:")
print(f)

# Lösen Sie das Gleichungssystem K_Matrix * W = f
W = solve(K_Matrix, f)

# Zeigen Sie die Lösung, d.h., die Werte für die Verschiebungen der Knotenpunkte
print("Lösung (Verschiebungen der Knotenpunkte) W:")
print(W)

print(k_b)
print(k_s)

# Laden der gespeicherten Matrix

print(Z)