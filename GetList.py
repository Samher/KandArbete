import numpy as np
import csv
import matplotlib.pyplot as plt
import matplotlib
import argparse
import requests
import os
import fractions


def spec2num(l):
    match(l.upper()):
        case "S":
            lnum = 0
        case "P":
            lnum = 1
        case "D": 
            lnum = 2
        case _:
            lnum = ord(l)-67
    return lnum


class LevelList:
    def __init__(self):
        self._list = {}
        self._warnings = []

    @property
    def list(self):
        return self._list

    @property
    def warnings(self):
        return "\n".join(self._warnings)

    # ----------- ORDER OF QNUMS -----------
    # LS: J, L, S, P
    # JK: J, K, Jc, S2, P
    # JJ: J, J1, J2, P
    # LK: J, L, S2, K, P
    # UC: J, P
    def add_level(self, coupling, label, energy, qnums):
        J = qnums[0]
        P = qnums[-1]
        match coupling:
            case "LS":
                L = qnums[1]
                S = qnums[2]
                tag = f"{label};{coupling};{J};{L};{S};{P}"

            case "JK":
                K = qnums[1]
                if len(qnums) == 4:
                    Jc = None
                    S2 = qnums[2]
                else:
                    Jc = qnums[2]
                    S2 = qnums[3]
                tag = f"{label};{coupling};{J};{K};{Jc};{S2};{P}"

            case "JJ":
                J1 = qnums[1]
                J2 = qnums[2]
                tag = f"{label};{coupling};{J};{J1};{J2};{P}"

            case "LK":
                if len(qnums) == 4:
                    L = None
                    S2 = qnums[2]
                    K = qnums[3]
                else:
                    L = qnums[2]
                    S2 = qnums[3]
                    K = qnums[4]
                tag = f"{label};{coupling};{J};{L};{S2};{K};{P}"

            case "UC":
                tag = f"{label};{coupling};{J};{P}"


        if tag not in self._list:  # If tag does not exist, add it and associate the energy
            self._list[tag] = [energy]

        else:
            closest_energy = min(self._list[tag], key=lambda E:abs(E-energy))
            closest_energy_index = np.argmin(np.abs(np.array(self._list[tag]) - energy))

            if energy == closest_energy:  # If energy already exists for given tag, do nothing
                pass

            elif abs(closest_energy - energy) < 1:  # If energy is very close to another energy for given tag, give warning
                warning = f"Tag: {tag}, Index: {closest_energy_index}, Energies: {energy}, {closest_energy}"
                if warning not in self._warnings:
                    self._warnings.append(warning)

            else:  # If energy is unique and sufficiently distant from others, add it as a new index
                self._list[tag].append(energy)
        
        return f"{tag},{len(self._list[tag])-1}"

    def __str__(self):
        s = []
        for tag in self._list:
            s.append([tag, self._list[tag]])
        s.sort(key=lambda E:E[1][0])
        for i in range(len(s)):
            s[i] = s[i][0] + "," + ",".join(map(str, s[i][1]))
        s = "\n".join(s)
        return s
    

# Parse arguments from command line
parser = argparse.ArgumentParser()
parser.add_argument("-Z", type=int, help="atomic number")
parser.add_argument("-I", type=int, help="level of ionization")
args = parser.parse_args()
Z = args.Z
I = args.I
FOLDER = f"Output\\{Z}0{I}\\"
FILENAME = f"parse_output_{Z}0{I}.txt"
FILENAME_NIST = f"{FOLDER}NIST_data_{Z}0{I}.txt"

if not os.path.exists(FOLDER):
    os.makedirs(FOLDER)

# Read and parse VALD data
levels = LevelList()
trans_levels = []
transitions = []
with open(FILENAME, newline='') as file:
    reader = csv.reader(file, delimiter=',')
    for row in reader:
        if row[0] in ["Autoionization", "Allowed_transition", "Forbidden_transition"]:
            transitions.append([row[0], trans_levels[1], trans_levels[0], row[1], row[2], row[3], row[4]])
            trans_levels = []
            continue

        label = row[0]
        energy = float(row[1])
        coupling_scheme = row[2]
        quantum_numbers = row[3:]

        current_tag = levels.add_level(coupling_scheme, label, energy, quantum_numbers)
        trans_levels.append(current_tag)


# Write energy levels to file
file = open(f"{FOLDER}energy_levels_{Z}0{I}.txt", "w")
file.write(levels.__str__())
file.close()

# Write warnings to file
file = open(f"{FOLDER}warnings_{Z}0{I}.txt", "w")
file.write(levels.warnings)
file.close()

# Write transitions to file
for i in range(len(transitions)):
    transitions[i] = ",".join(transitions[i])
transitions_text = "\n".join(transitions)
file = open(f"{FOLDER}transitions_{Z}0{I}.txt", "w")
file.write(transitions_text)
file.close()

# If NIST data is not available, download it
if not os.path.isfile(FILENAME_NIST):
    nist_data = requests.get(f'https://physics.nist.gov/cgi-bin/ASD/energy1.pl?de=0&spectrum=Z%3D{Z}+{I}&submit=Retrieve+Data&units=0&format=2&output=0&page_size=15&multiplet_ordered=1&conf_out=on&term_out=on&level_out=on&j_out=on&temp=')
    content = nist_data.text.replace("\"","")
    content = content.replace("=","")
    lines = content.splitlines()

    file = open(FILENAME_NIST, "w")
    for l in lines[1:]:
        if "Limit" in l:
            break
        file.write(l+"\n")
    file.close()

# Parse NIST data
nist_levels = LevelList()
with open(FILENAME_NIST, newline='') as file:
    reader = csv.reader(file, delimiter=',')

    # ----------- NIST order -----------
    # LS: label, letter (2S+1)L*, J, , energy, ,
    # JK: label <Jc>, (2S_2+1)[K]*, J, , energy, , 
    # JJ: label (no period before <J2>), (J1, J2)*, J, , energy, , 
    # LK?: label, L (2S_2+1)[K]*, J, , energy, , 
    # UC: ev. label, * or 12345e, J, , energy, , 
    for row in reader:
        label = row[0].replace("?","")
        term = row[1].replace("?","")
        if "(" in term:  # JJ case
            J = float(fractions.Fraction(row[3]))
            J1 = float(fractions.Fraction(row[1][1:]))
            J2 = float(fractions.Fraction(row[2][:-2]))
            energy = float(row[5])
            label = label[:label.rfind("<")] + "." + label[label.rfind("<"):]
            if row[2][-1] == "*":
                P = 1
            else: 
                P = 2   
            nist_levels.add_level("JJ", label, energy, [J, J1, J2, P])         
            

        elif "[" in term:
            J = float(fractions.Fraction(row[2]))
            energy = float(row[4])
            if term[-1] == "*":
                P = 1
            else: 
                P = 2

            if term[0].isalpha():  # LK case, still needs to be tested
                L = spec2num(term[0])
                S2 = (float(term[2]) - 1)/2
                K = float(fractions.Fraction(term[term.find("[")+1:term.find("]")]))
                nist_levels.add_level("LK", label, energy, [J, L, S2, K, P])

            else:  # JK case
                K = float(fractions.Fraction(term[term.find("[")+1:term.find("]")]))
                Jc = float(fractions.Fraction(label[label.find("<")+1:label.find(">")]))
                S2 = (float(term[0]) - 1)/2
                nist_levels.add_level("JK", label, energy, [J, K, Jc, S2, P])


        elif len(term) > 1 and term[-1] != "e":  # LS case
            J = float(fractions.Fraction(row[2]))
            energy = float(row[4])
            if term[0].islower():
                S = (float(term[2]) - 1)/2
            else:
                S = (float(term[0]) - 1)/2
            if term[-1] == "*":
                P = 1
                L = spec2num(term[-2])
            else: 
                P = 2
                L = spec2num(term[-1])
            nist_levels.add_level("LS", label, energy, [J, L, S, P])

        else:  # UC case
            pass
             

# Plot energy levels as horizontal lines
fig, axs = plt.subplots(1,2, sharey=True, figsize = (13,8))
for key in levels.list:
    for i in range(len(levels.list[key])):
        axs[0].axhline(float(levels.list[key][i]), color='k', linewidth='0.5')

for key in nist_levels.list:
    for i in range(len(nist_levels.list[key])):
        axs[1].axhline(float(nist_levels.list[key][i]), color='r', linewidth='0.5')

# Configure the plot
font = {'family' : 'serif',
        'size' : 22}

matplotlib.rc('font', **font)
axs[0].set_ylabel("Energy (cm^-1)")
axs[0].set_xticks([])
axs[1].set_xticks([])
axs[0].set_title("VALD")
axs[1].set_title("NIST")
# plt.savefig("diagram.pdf")
plt.show()
