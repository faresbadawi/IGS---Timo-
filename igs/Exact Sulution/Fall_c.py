import numpy as np
import matplotlib.pyplot as plt

# Gegebene Parameter
E = 1e7  # Elastizitätsmodul (Young's Modulus)
b = 0.1  # Breite des Balkens
h = 0.01  # Dicke des Balkens
nu = 0.2  # Poisson's Ratio
alpha = 5/6  # Schubkorrektionsfaktor

# Berechnung der Biegesteifigkeit Kb
K_b = (E * b * h**3) / (12)

# Berechnung der Schubsteifigkeit Ks
K_s = alpha * (E * b * h) / 2*(1 + nu)

# Ausgabe der berechneten Werte
print("Biegesteifigkeit Kb:", K_b)
print("Schubsteifigkeit Ks:", K_s)

# Erstellen eines Arrays von x-Werten
x = np.linspace(0, 1, 100)

# Berechnung der Funktion w(x)
w = -(1 / (24 * K_b)) * (x**4 - 4 * x**3 + 6 * x**2) - (1 / (2 * K_s)) * (-x**2 + 2 * x)

# Finden des Minimums und seiner Position
min_w = np.min(w)
min_x = x[np.argmin(w)]

# Runden auf fünf Dezimalstellen
rounded_min_w = round(min_w, 5)

# Plot
fig, ax = plt.subplots()
ax.plot(x, w)
ax.set_xlabel('x')
ax.set_ylabel('W(x)')
ax.set_title('Eingespannter Träger')

# Markierung des Minimums als Punkt mit gerundetem Wert
ax.scatter(min_x, min_w, color='red', label=f'Min: {rounded_min_w:.5f}')
ax.legend()

# Gestrichelte Linie bei y = 0 hinzufügen
ax.axhline(0, color='gray', linestyle='--')

# Speichern des Bildes in hoher Qualität (DPI erhöht)
plt.savefig(r'C:\Users\Fares Badawi\Desktop\4. Semester\Studienarbeit\Bilder\Exact sultion\Exact_solution_c.png', dpi=500, bbox_inches='tight')
plt.show()

