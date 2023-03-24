import numpy as np
import csv
import matplotlib.pyplot as plt
import matplotlib
import argparse


font = {'family' : 'serif',
        'size' : 22}

matplotlib.rc('font', **font)

def partition_function(eLevels, J, T):
    Z = []
    h = 4.135667696e-15  # in eV*s
    c = 2.99792458e10  # in cm/s

    for i in range(len(eLevels)):
        e_eV = eLevels[i] * h*c
        Z.append((2*J[i]+1)*np.exp(-e_eV / T))

    return np.sum(Z), Z

parser = argparse.ArgumentParser()
parser.add_argument("-Z", type=int, help="atomic number")
parser.add_argument("-I", type=int, help="level of ionization")
args = parser.parse_args()
Zn = args.Z
I = args.I
FILENAME = f"Output\\{Zn}0{I}\\energy_levels_{Zn}0{I}.txt"
energies = []
with open(FILENAME, newline='') as file:
    reader = csv.reader(file, delimiter=',')
    for row in reader:
        tag = row[0].split(';')
        for i in row[1:]:
            energies.append([float(i), float(tag[2])])

energies = np.array(energies)
T = np.linspace(0.001, 2, 100)
Z = np.zeros(100)
for i in range(len(Z)):
    Z[i] = partition_function(energies[:,0], energies[:,1], T[i])[0]

NIST_Zs_02600 = np.array([21.71, 30.89, 47.99, 78.54, 127.43, 196.65, 284.88, 388.89])
NIST_Zs_05800 = np.array([71.69, 248.81, 539.72, 905.06, 1296.81, 1683.94, 2050.80, 2391.02])
NIST_Zs_05801 = np.array([77.41, 237.64, 412.46, 584.43, 750.13, 907.72, 1055.98, 1194.34])
NIST_Ts = np.array([0.25, 0.5, 0.75, 1, 1.25, 1.5, 1.75, 2])

plt.figure(figsize=(10,8))
plt.plot(T, Z)
plt.xlabel("T (eV)")
plt.ylabel("Z")
plt.title(f"Partition function for {Zn}0{I}")
#plt.scatter(NIST_Ts, NIST_Zs_05801)
plt.show()


#print(f"Z={Z[0]}")

# x = np.arange(0, len(Z[1]))
# plt.scatter(x, Z[1])
# plt.show()
