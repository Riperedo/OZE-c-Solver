#!/usr/bin/env python3
"""
Generate publication-quality 2D Heatmaps of Structure Factor RMSE
across Packing Fraction phi and Temperature T*.
"""

import os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm

BASE_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DATA_DIR = os.path.join(BASE_DIR, "data")
PLOTS_DIR = os.path.join(BASE_DIR, "plots")
os.makedirs(PLOTS_DIR, exist_ok=True)

# Grid parameters
temperatures = [10.0, 1.0, 0.10, 0.01]
phi_values = [0.1, 0.2, 0.3, 0.4, 0.5]

# Metric matrices: rows=T, cols=phi
rmse_s000 = np.zeros((len(temperatures), len(phi_values)))
rmse_s110 = np.zeros((len(temperatures), len(phi_values)))
rmse_s112 = np.zeros((len(temperatures), len(phi_values)))
rmse_s1   = np.zeros((len(temperatures), len(phi_values)))

for i, T in enumerate(temperatures):
    for j, phi in enumerate(phi_values):
        ana_path = os.path.join(DATA_DIR, f"ana_phi_{phi:.1f}_T_{T:.2f}.dat")
        num_path = os.path.join(DATA_DIR, f"num_phi_{phi:.1f}_T_{T:.2f}.dat")

        if not os.path.exists(ana_path) or not os.path.exists(num_path):
            print(f"Warning: Missing files for T={T}, phi={phi}")
            continue

        ana = np.loadtxt(ana_path)
        num = np.loadtxt(num_path)

        k_ana = ana[:, 0]
        k_num = num[:, 0]

        mask = (k_ana >= 0.5) & (k_ana <= 20.0)
        k_eval = k_ana[mask]

        # S000 (col 1)
        ana_s000 = ana[mask, 1]
        num_s000 = np.interp(k_eval, k_num, num[:, 1])
        rmse_s000[i, j] = np.sqrt(np.mean((ana_s000 - num_s000)**2))

        # S110 (col 2)
        ana_s110 = ana[mask, 2]
        num_s110 = np.interp(k_eval, k_num, num[:, 2])
        rmse_s110[i, j] = np.sqrt(np.mean((ana_s110 - num_s110)**2))

        # S112 (col 3)
        ana_s112 = ana[mask, 3]
        num_s112 = np.interp(k_eval, k_num, num[:, 3])
        rmse_s112[i, j] = np.sqrt(np.mean((ana_s112 - num_s112)**2))

        # S1 (col 5)
        ana_s1 = ana[mask, 5]
        num_s1 = np.interp(k_eval, k_num, num[:, 5])
        rmse_s1[i, j] = np.sqrt(np.mean((ana_s1 - num_s1)**2))

# Plot styling
plt.rcParams.update({
    "font.family": "sans-serif",
    "font.sans-serif": ["DejaVu Sans", "Helvetica", "Arial"],
    "font.size": 11,
    "axes.labelsize": 13,
    "axes.titlesize": 13,
    "xtick.labelsize": 11,
    "ytick.labelsize": 11,
})

fig, axes = plt.subplots(2, 2, figsize=(11.5, 8.5), constrained_layout=True)

panels = [
    (axes[0, 0], rmse_s110, r"(a) Dipolar Projection $\mathrm{RMSE}(S^{110})$", 1e-6, 1e1),
    (axes[0, 1], rmse_s112, r"(b) Dipolar Anisotropy $\mathrm{RMSE}(S^{112})$", 1e-6, 1e1),
    (axes[1, 0], rmse_s1,   r"(c) Decoupled Mode $\mathrm{RMSE}(S^1)$",      1e-6, 1e1),
    (axes[1, 1], rmse_s000, r"(d) Hard-Core Base $\mathrm{RMSE}(S^{000})$",  1e-3, 1e-1),
]

y_labels = [f"$T^* = {t:.2f}$" for t in temperatures]
x_labels = [f"$\\phi = {p:.1f}$" for p in phi_values]

for ax, mat, title, vmin, vmax in panels:
    im = ax.imshow(mat, cmap="viridis_r", norm=LogNorm(vmin=vmin, vmax=vmax), aspect="auto")
    ax.set_title(title, fontweight="bold", pad=10)
    ax.set_xticks(range(len(phi_values)))
    ax.set_xticklabels(x_labels)
    ax.set_yticks(range(len(temperatures)))
    ax.set_yticklabels(y_labels)

    # Annotate exact numbers in each cell
    for i in range(len(temperatures)):
        for j in range(len(phi_values)):
            val = mat[i, j]
            log_val = np.log10(val)
            log_vmin = np.log10(vmin)
            log_vmax = np.log10(vmax)
            norm_val = (log_val - log_vmin) / (log_vmax - log_vmin)
            text_color = "white" if norm_val > 0.65 else "black"

            exp = int(np.floor(np.log10(val)))
            coeff = val / (10**exp)
            text = f"{coeff:.1f}e{exp:+d}"
            ax.text(j, i, text, ha="center", va="center", color=text_color, fontsize=10, fontweight="semibold")

    cbar = fig.colorbar(im, ax=ax, shrink=0.85, pad=0.03)
    cbar.set_label("RMSE (log scale)", rotation=270, labelpad=15)

out_pdf = os.path.join(PLOTS_DIR, "fig6_rmse_heatmap.pdf")
out_png = os.path.join(PLOTS_DIR, "fig6_rmse_heatmap.png")

plt.savefig(out_pdf, dpi=300)
plt.savefig(out_png, dpi=200)
plt.close()

print(f"Heatmaps successfully generated:\n  - {out_pdf}\n  - {out_png}")
