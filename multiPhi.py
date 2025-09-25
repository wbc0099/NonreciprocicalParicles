#!/usr/bin/env python3

import os
import numpy as np
import matplotlib.pyplot as plt
import tools  # your own helpers
from matplotlib import font_manager as fm

# ------------------------ Config ------------------------
path = "../conf-data/"
# Optional: image size & output
figsize = (6.0, 5.0)
out_png = "./PS.png"
dpi = 300

## Set font globally
#fontPath="/usr/share/fonts/times.ttf"
#times=fm.FontProperties(fname=fontPath)
#plt.rcParams["font.family"]=times.get_name()
plt.rcParams["font.family"]="Times New Roman"
plt.rcParams["font.size"]=20
plt.rc('text',usetex=True)

# -------------------- Collect folders -------------------
items = [item for item in os.listdir(path) if item.startswith("N")]
# Sort by (proportion, kba)
def parse_prop(folder):
    # expects ... "proportion<val>--"
    return float(folder.split("proportion")[1].split("--")[0])
def parse_kba(folder):
    # expects ... "kba-<val>" at the end
    return float(folder.split("kba-")[1])

items = sorted(items, key=lambda x: (parse_prop(x), parse_kba(x)))
print("Folders:", items)

# -------------------- Scan & compute --------------------
# We’ll store one scalar per (kba, proportion)
records = []  # list of dicts: {prop, kba, value}

for item in items:
    proportion = parse_prop(item)
    kba = parse_kba(item)
    print(f"{item}  -> proportion={proportion}, kba={kba}")

    # Files
    data_path = os.path.join(path, item, "conf_1000.dat")
    rc_path   = os.path.join(path, item, "run_cpu.pl")

    # Read data & parameters
    data = tools.readFile(data_path)
    Na = tools.readParameter(rc_path, ["$Na"])[0]
    Nb = tools.readParameter(rc_path, ["$Nb"])[0]
    redParticleNumber  = int(Na)  # strip trailing ';'
    blueParticleNumber = int(Nb)
    totalN = redParticleNumber + blueParticleNumber

    # Compute order parameters (no images drawn)
    phiRed = tools.phi(
        data[:redParticleNumber],
        50, 50, threshold=2,
        filePath=os.path.join(path, item, "phi6.png"),
        order=6, drawPciture=False, returnSumPhi=True
    )
    phiBlue = tools.phi(
        data[redParticleNumber:],
        50, 50, threshold=1.2,
        filePath=os.path.join(path, item, "phi3.png"),
        order=3, drawPciture=False, returnSumPhi=True
    )

    print(f"  redParticleNumber={redParticleNumber}")
    print(f"  phiRed={phiRed}, phiBlue={phiBlue}, "
          f"phiSumDN={phiRed*redParticleNumber + phiBlue*blueParticleNumber}")

    # Your scalar per state (averaged per particle)
    value = (phiRed + phiBlue) / totalN
    records.append({"prop": proportion, "kba": kba, "val": value})

# ------------------ Build 2D grid (Ny x Nx) --------------
# Unique sorted coordinate values
prop_unique = sorted({r["prop"] for r in records})  # x-axis
kba_unique  = sorted({r["kba"]  for r in records})  # y-axis

Nx = len(prop_unique)
Ny = len(kba_unique)
print(f"Grid size: Ny={Ny} (kba) x Nx={Nx} (proportion)")

# Map value -> [row (kba index), col (prop index)]
prop_to_ix = {p:i for i, p in enumerate(prop_unique)}
kba_to_iy  = {k:i for i, k in enumerate(kba_unique)}

N = np.full((Ny, Nx), np.nan, dtype=float)
for r in records:
    iy = kba_to_iy[r["kba"]]
    ix = prop_to_ix[r["prop"]]
    N[iy, ix] = r["val"]

print("N.shape =", N.shape)

# ------------------------- Plot --------------------------
# IMPORTANT: we let imshow use *index coordinates*,
# then place ticks at those indices and set labels to your actual values.
fig, ax = plt.subplots(figsize=figsize)
im = ax.imshow(
    N, cmap="inferno", origin="lower", interpolation="nearest"  # ticks at pixel centers by default
)

# Centered ticks: use integer indices (0..Nx-1, 0..Ny-1)
ax.set_xticks(np.arange(1,Nx-1,2))
ax.set_yticks(np.arange(1,Ny-1,2))

# Label ticks with actual parameter values
# (format as you like; here keep raw or use f"{v:.3g}")
xticks=[0.4,0.8,1.2,1.6]
yticks=[0.03,0.07,0.11,0.15]
ax.set_xticklabels([f"{v:g}" for v in xticks], rotation=45, ha="right")
ax.set_yticklabels([f"{v:g}" for v in yticks])
# ax.set_xticks([0.4,0.8,1.2,1.6])
# ax.set_yticks([0.03,0.07,0.11,0.15])

ax.set_xlabel(r"$\eta$")
ax.set_ylabel(r"$k_{ba}$")

cbar = plt.colorbar(im, ax=ax)
cbar.set_label(r"$\bar{\Psi}$")  # rename if you have a better name

#ax.set_title("Heatmap of value vs (proportion, kba)")
fig.tight_layout()
plt.savefig(out_png, dpi=dpi)
plt.close()
print(f"Saved: {out_png}")

os.system(f"eog {out_png}")