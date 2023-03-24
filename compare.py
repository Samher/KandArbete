import numpy as np


def find_diffs(l1, l2):
    diffs = []
    for e1 in l1:
        closest_e = min(l2, key=lambda e2:abs(e2-e1))
        if np.abs(closest_e-e1) >= 1e-2*e1:
            diffs.append(e1)
    for e2 in l2:
        closest_e = min(l1, key=lambda e1:abs(e1-e2))
        if np.abs(closest_e-e2) >= 1e-2*e2:
            diffs.append(e2)
    return diffs


a = [1, 2, 3, 4, 5, 6, 7, 8]
b = [1, 2, 2.001, 3.5, 4, 5.001, 7.99, 9]

print(find_diffs(a,b))
