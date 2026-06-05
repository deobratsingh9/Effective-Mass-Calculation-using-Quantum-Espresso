#!/usr/bin/env python3

import numpy as np
import matplotlib.pyplot as plt
import re
import os

# ==========================================================
# FILES
# ==========================================================

BAND_FILE    = "band_files/output_band.dat"
GNUPLOT_FILE = "gnuplot_files/gnuplot.tmp_band"

FIT_POINTS = 5

# ==========================================================
# CONSTANTS
# ==========================================================

hbar = 1.054571817e-34
eV_to_J = 1.602176634e-19
m0 = 9.1093837015e-31

bohr_to_ang = 0.52917721092
ang_to_m = 1e-10

celldm1 = 10.26
a_m = celldm1 * bohr_to_ang * ang_to_m
G = 2*np.pi/a_m

# ==========================================================
# READ BAND.DAT
# ==========================================================

with open(BAND_FILE) as f:
    lines = f.readlines()

nbnd = int(lines[0].split("nbnd=")[1].split(",")[0])
nks  = int(lines[0].split("nks=")[1].split("/")[0])

klist = []
bands = []

i = 1
while i < len(lines):

    kxyz = np.array(list(map(float, lines[i].split())))
    i += 1

    energies = []
    while len(energies) < nbnd:
        energies.extend(map(float, lines[i].split()))
        i += 1

    klist.append(kxyz)
    bands.append(np.array(energies[:nbnd]))

klist = np.array(klist)
bands = np.array(bands)

# ==========================================================
# GNUPLOT PARSE
# ==========================================================

with open(GNUPLOT_FILE) as f:
    txt = f.read()

xscale = float(re.search(r'xscale=\s*([0-9Ee.+\-]+)', txt).group(1))
eref   = float(re.findall(r'eref=\s*([0-9Ee.+\-]+)', txt)[-1])

bands -= eref

# ==========================================================
# FIND VBM / CBM
# ==========================================================

VBM = -1e20
CBM =  1e20

for ib in range(nbnd):

    E = bands[:, ib]

    occ = np.where(E <= 0)[0]
    if len(occ):
        vmax = E[occ].max()
        if vmax > VBM:
            VBM = vmax
            VBM_band = ib
            VBM_kidx = occ[np.argmax(E[occ])]

    emp = np.where(E > 0)[0]
    if len(emp):
        cmin = E[emp].min()
        if cmin < CBM:
            CBM = cmin
            CBM_band = ib
            CBM_kidx = emp[np.argmin(E[emp])]

# ==========================================================
# EFFECTIVE MASS
# ==========================================================

def effective_mass(band_idx, k_idx):

    E = bands[:, band_idx]
    idxs = np.arange(max(0,k_idx-3), min(nks,k_idx+4))

    kk = klist[idxs]
    EE = E[idxs]

    dk = kk - klist[k_idx]

    for v in dk:
        if np.linalg.norm(v) > 1e-12:
            direction = v/np.linalg.norm(v)
            break

    kproj = np.dot(dk, direction)

    coef = np.polyfit(kproj, EE, 2)

    d2 = 2*coef[0]
    d2SI = d2 * eV_to_J / (G**2)

    mstar = hbar**2 / abs(d2SI)

    return mstar/m0, idxs

hole_mass, hole_idxs = effective_mass(VBM_band, VBM_kidx)
electron_mass, electron_idxs = effective_mass(CBM_band, CBM_kidx)

# ==========================================================
# LOAD GNUPLOT BAND FILES (FIXED PLOTTING)
# ==========================================================

files = re.findall(r'"([^"]*output_pband\.dat[^"]*)"', txt)

band_data = []
kref = None

for f in files:
    if os.path.isfile(f):
        d = np.loadtxt(f)
        k = d[:,0]*xscale
        E = d[:,1]-eref

        # ✅ BREAK SEGMENTS HERE (important fix)
        kplot = [k[0]]
        Eplot = [E[0]]

        dk = np.diff(k)
        threshold = 5*np.median(np.abs(dk))

        for i in range(1,len(k)):
            if abs(k[i] - k[i-1]) > threshold:
                kplot.append(np.nan)
                Eplot.append(np.nan)

            kplot.append(k[i])
            Eplot.append(E[i])

        band_data.append((np.array(kplot), np.array(Eplot)))

        if kref is None:
            kref = k

# ==========================================================
# SYMMETRY LABELS
# ==========================================================

sym_x = []
sym_lab = []

for line in txt.splitlines():
    if "set label" not in line:
        continue

    m1 = re.search(r'"(.*?)"', line)
    m2 = re.search(r'at\s+([0-9.]+)', line)

    if m1 and m2:
        label = m1.group(1).replace("{/Symbol G}", "Γ")
        sym_lab.append(label)
        sym_x.append(float(m2.group(1))*xscale)

# ==========================================================
# FINAL PLOT (CORRECT)
# ==========================================================

plt.figure(figsize=(10,7))

# bands
for k,E in band_data:
    plt.plot(k, E, lw=2)

# markers (correct mapping)
plt.scatter(kref[VBM_kidx], VBM, color='red', s=120, label="VBM")
plt.scatter(kref[CBM_kidx], CBM, color='blue', s=120, label="CBM")

plt.scatter(kref[hole_idxs],
            bands[hole_idxs, VBM_band],
            color='green', s=70, label="hole fit")

plt.scatter(kref[electron_idxs],
            bands[electron_idxs, CBM_band],
            color='orange', s=70, label="electron fit")

plt.xticks(sym_x, sym_lab)

plt.axhline(0, ls='--')

plt.xlabel("k-path")
plt.ylabel("Energy (eV)")
plt.title("Band Structure (FIXED)")

plt.legend()
plt.tight_layout()

plt.savefig("band_structure_fixed.png", dpi=300)
plt.close()

# ==========================================================
# OUTPUT
# ==========================================================

print("Gap =", CBM-VBM, "eV")
print("Hole mass =", hole_mass, "m0")
print("Electron mass =", electron_mass, "m0")
