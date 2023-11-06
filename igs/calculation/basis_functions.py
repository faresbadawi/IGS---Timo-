# basis_functions.py

def basisfunktion(t, knotenvektor, i, p):
    if knotenvektor[-1] == 1.0 and i == len(knotenvektor) - p - 1:
        return 1.0
    else:
        if p == 0:
            return 1 if knotenvektor[i] <= t < knotenvektor[i + 1] else 0
        else:
            nenner1 = knotenvektor[i + p] - knotenvektor[i]
            nenner2 = knotenvektor[i + p + 1] - knotenvektor[i + 1]
            a = 0 if nenner1 == 0 else (t - knotenvektor[i]) / nenner1
            b = 0 if nenner2 == 0 else (knotenvektor[i + p + 1] - t) / nenner2
            return a * basisfunktion(t, knotenvektor, i, p - 1) + b * basisfunktion(t, knotenvektor, i + 1, p - 1)

def ableitungsbasisfunktion(t, knotenvektor, i, p):
    if p == 0:
        return 0
    else:
        denom1 = knotenvektor[i + p] - knotenvektor[i]
        denom2 = knotenvektor[i + p + 1] - knotenvektor[i + 1]
        a = 0 if denom1 == 0 else p / denom1
        b = 0 if denom2 == 0 else -p / denom2
        return a * basisfunktion(t, knotenvektor, i, p - 1) + b * basisfunktion(t, knotenvektor, i + 1, p - 1)        
        
def ableitung_kter_ordnung(t, knotenvektor, i, p, k):
    if k == 0:
        return ableitungsbasisfunktion(t, knotenvektor, i, p)
    else:
        nenner1 = knotenvektor[i + p] - knotenvektor[i]
        nenner2 = knotenvektor[i + p + 1] - knotenvektor[i + 1]
        a = 0 if nenner1 == 0 else p / nenner1
        b = 0 if nenner2 == 0 else -p / nenner2
        return a * ableitung_kter_ordnung(t, knotenvektor, i, p - 1, k - 1) + b * ableitung_kter_ordnung(t, knotenvektor, i + 1, p - 1, k - 1)
