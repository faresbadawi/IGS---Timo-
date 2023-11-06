# matrix_operations.py

import numpy as np
def determinante_und_inverse_berechnen(ableitungsbasisfunktion_matrix, kontrollpunkte_vektor):
    if len(ableitungsbasisfunktion_matrix[0]) != len(kontrollpunkte_vektor):
        raise ValueError("Die Anzahl der Spalten in der Basisfunktionsmatrix stimmt nicht mit der Anzahl der Kontrollpunkte überein.")
    
    det_Jx = np.dot(ableitungsbasisfunktion_matrix, kontrollpunkte_vektor)
    inverse_Jx = 1 / det_Jx
    return det_Jx, inverse_Jx