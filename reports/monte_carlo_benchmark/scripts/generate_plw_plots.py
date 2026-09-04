#!/usr/bin/env python3
"""
Publication Figure Generator for Patey, Levesque & Weis (1979) Monte Carlo Comparisons.
Generates vector PDFs and high-resolution PNGs matching academic journal aesthetics.
"""

import os
import numpy as np
import matplotlib.pyplot as plt

# Styling setup
plt.rcParams.update({
    'font.family': 'sans-serif',
    'font.size': 11,
    'axes.labelsize': 12,
    'axes.titlesize': 12,
    'xtick.labelsize': 10,
    'ytick.labelsize': 10,
    'legend.fontsize': 10,
    'figure.titlesize': 14,
    'lines.linewidth': 1.8,
    'grid.alpha': 0.3,
    'grid.linestyle': '--'
})

WORKSPACE_DIR = "/home/jinzo/Desktop/codigos/OZE_c_solver"
DATA_DIR = os.path.join(WORKSPACE_DIR, "reports", "monte_carlo_benchmark", "data")
PLW_MC_DIR = os.path.join(WORKSPACE_DIR, "test", "data_Monte_Carlo_Simulations", "PateyLevesqueWeis_1979")
PLOTS_DIR = os.path.join(WORKSPACE_DIR, "reports", "monte_carlo_benchmark", "plots")

os.makedirs(PLOTS_DIR, exist_ok=True)

# Color scheme
C_MC = '#111111'     # Black circles for Monte Carlo
C_MSA = '#D9534F'    # Crimson Red (dashed)
C_LHNC = '#F0AD4E'   # Amber / Orange (dash-dot)
C_QHNC = '#0275D8'   # Cobalt Blue (solid)
C_RHNC = '#9467BD'   # Purple (dotted)

def load_sol(state_key, closure):
    fpath = os.path.join(DATA_DIR, f"solution_plw_{state_key}_{closure}.dat")
    if not os.path.exists(fpath):
        fpath = os.path.join(DATA_DIR, f"solution_{state_key}_{closure}.dat")
    if os.path.exists(fpath):
        data = np.loadtxt(fpath)
        return {"r": data[:, 0], "g000": data[:, 1] + 1.0, "h110": data[:, 2], "h112": data[:, 3]}
    return None

def plot_single_g000(state_key, mc_filename, fig_num_plw, rho_val, mu2_val, out_base, xlim=(0.9, 3.5), ylim=(0.0, 5.5)):
    fig, ax = plt.subplots(figsize=(6.5, 4.8), dpi=300)
    
    # Load MC
    mc_file = os.path.join(PLW_MC_DIR, mc_filename)
    if os.path.exists(mc_file):
        mc = np.loadtxt(mc_file)
        ax.scatter(mc[:, 0], mc[:, 1], color=C_MC, s=36, facecolors='none', edgecolors='black', linewidth=1.5, zorder=5, label='Monte Carlo (PLW 1979)')

    # Load Closures
    sol_msa = load_sol(state_key, "MSA")
    sol_lhnc = load_sol(state_key, "LHNC")
    sol_qhnc = load_sol(state_key, "QHNC")
    sol_rhnc = load_sol(state_key, "RHNC")

    if sol_msa:
        mask = (sol_msa["r"] >= 1.0) & (sol_msa["r"] <= xlim[1] + 0.5)
        ax.plot(sol_msa["r"][mask], sol_msa["g000"][mask], color=C_MSA, linestyle='--', label='MSA')

    if sol_lhnc:
        mask = (sol_lhnc["r"] >= 1.0) & (sol_lhnc["r"] <= xlim[1] + 0.5)
        ax.plot(sol_lhnc["r"][mask], sol_lhnc["g000"][mask], color=C_LHNC, linestyle='-.', label='LHNC')

    if sol_qhnc:
        mask = (sol_qhnc["r"] >= 1.0) & (sol_qhnc["r"] <= xlim[1] + 0.5)
        ax.plot(sol_qhnc["r"][mask], sol_qhnc["g000"][mask], color=C_QHNC, linestyle='-', linewidth=2.0, label='QHNC')

    if sol_rhnc:
        mask = (sol_rhnc["r"] >= 1.0) & (sol_rhnc["r"] <= xlim[1] + 0.5)
        ax.plot(sol_rhnc["r"][mask], sol_rhnc["g000"][mask], color=C_RHNC, linestyle=':', linewidth=2.2, label='RHNC')

    ax.set_xlabel(r'$r / \sigma$', fontweight='bold')
    ax.set_ylabel(r'$g^{000}(r)$', fontweight='bold')
    ax.set_title(f"Patey et al. (1979) Fig. {fig_num_plw}: Radial Distribution $g^{{000}}(r)$\n$\\rho^* = {rho_val}$, $\\mu^{{*2}} = {mu2_val}$, $T^* = 1.0$", fontsize=11, fontweight='bold')
    ax.set_xlim(xlim)
    ax.set_ylim(ylim)
    ax.grid(True)
    ax.legend(frameon=True, loc='upper right', framealpha=0.9)
    plt.tight_layout()

    pdf_path = os.path.join(PLOTS_DIR, f"{out_base}.pdf")
    png_path = os.path.join(PLOTS_DIR, f"{out_base}.png")
    plt.savefig(pdf_path)
    plt.savefig(png_path)
    plt.close()
    print(f"✓ Generated {out_base}")

def plot_angular_projections_plw():
    """Generates a 2x3 panel figure for h110 and h112 across rho = 0.15, 0.40, 0.60"""
    fig, axes = plt.subplots(2, 3, figsize=(14, 8), dpi=300, sharex='col')
    
    densities = [
        {"key": "rho_0.15_mu2_2.0", "rho": 0.15, "mu2": 2.0, "h110_mc": "Fig6_rho_0.15_mu2_2.0.dat", "h112_mc": "Fig9_rho_0.15_mu2_2.0.dat", "fig_h110": "6", "fig_h112": "9", "ylim_h110": (-0.2, 1.2), "ylim_h112": (-0.2, 3.8)},
        {"key": "rho_0.4_mu2_2.75", "rho": 0.40, "mu2": 2.75, "h110_mc": "Fig7_rho_0.4_mu2_2.75.dat", "h112_mc": "Fig10_rho_0.4_mu_2.75.dat", "fig_h110": "7", "fig_h112": "10", "ylim_h110": (-0.3, 3.2), "ylim_h112": (-0.2, 4.5)},
        {"key": "rho_0.6_mu2_2.75", "rho": 0.60, "mu2": 2.75, "h110_mc": "Fig8_rho_0.6_mu2_2.75.dat", "h112_mc": "Fig11_rho_0.6_mu2_2.75.dat", "fig_h110": "8", "fig_h112": "11", "ylim_h110": (-0.3, 2.8), "ylim_h112": (-0.2, 3.8)}
    ]

    for col_idx, d in enumerate(densities):
        ax_top = axes[0, col_idx]
        ax_bot = axes[1, col_idx]

        # Top row: h110(r)
        mc_110 = os.path.join(PLW_MC_DIR, d["h110_mc"])
        if os.path.exists(mc_110):
            data_mc = np.loadtxt(mc_110)
            ax_top.scatter(data_mc[:, 0], data_mc[:, 1], color=C_MC, s=28, facecolors='none', edgecolors='black', linewidth=1.4, zorder=5, label='MC (PLW 1979)')

        sol_msa = load_sol(d["key"], "MSA")
        sol_lhnc = load_sol(d["key"], "LHNC")
        sol_qhnc = load_sol(d["key"], "QHNC")
        sol_rhnc = load_sol(d["key"], "RHNC")

        if sol_msa:
            mask = (sol_msa["r"] >= 1.0) & (sol_msa["r"] <= 4.0)
            ax_top.plot(sol_msa["r"][mask], sol_msa["h110"][mask], color=C_MSA, linestyle='--', label='MSA')
            ax_bot.plot(sol_msa["r"][mask], sol_msa["h112"][mask], color=C_MSA, linestyle='--', label='MSA')

        if sol_lhnc:
            mask = (sol_lhnc["r"] >= 1.0) & (sol_lhnc["r"] <= 4.0)
            ax_top.plot(sol_lhnc["r"][mask], sol_lhnc["h110"][mask], color=C_LHNC, linestyle='-.', label='LHNC')
            ax_bot.plot(sol_lhnc["r"][mask], sol_lhnc["h112"][mask], color=C_LHNC, linestyle='-.', label='LHNC')

        if sol_qhnc:
            mask = (sol_qhnc["r"] >= 1.0) & (sol_qhnc["r"] <= 4.0)
            ax_top.plot(sol_qhnc["r"][mask], sol_qhnc["h110"][mask], color=C_QHNC, linestyle='-', linewidth=2.0, label='QHNC')
            ax_bot.plot(sol_qhnc["r"][mask], sol_qhnc["h112"][mask], color=C_QHNC, linestyle='-', linewidth=2.0, label='QHNC')

        if sol_rhnc:
            mask = (sol_rhnc["r"] >= 1.0) & (sol_rhnc["r"] <= 4.0)
            ax_top.plot(sol_rhnc["r"][mask], sol_rhnc["h110"][mask], color=C_RHNC, linestyle=':', linewidth=2.2, label='RHNC')
            ax_bot.plot(sol_rhnc["r"][mask], sol_rhnc["h112"][mask], color=C_RHNC, linestyle=':', linewidth=2.2, label='RHNC')

        # Bottom row: h112(r)
        mc_112 = os.path.join(PLW_MC_DIR, d["h112_mc"])
        if os.path.exists(mc_112):
            data_mc = np.loadtxt(mc_112)
            ax_bot.scatter(data_mc[:, 0], data_mc[:, 1], color=C_MC, s=28, facecolors='none', edgecolors='black', linewidth=1.4, zorder=5, label='MC (PLW 1979)')

        ax_top.set_title(f"$\\rho^* = {d['rho']}$, $\\mu^{{*2}} = {d['mu2']}$ (Fig. {d['fig_h110']})", fontsize=11, fontweight='bold')
        ax_top.set_ylabel(r'$h^{110}(r)$', fontweight='bold')
        ax_top.set_xlim(0.9, 4.0)
        ax_top.set_ylim(d["ylim_h110"])
        ax_top.grid(True)
        if col_idx == 0:
            ax_top.legend(frameon=True, loc='upper right')

        ax_bot.set_title(f"$\\rho^* = {d['rho']}$, $\\mu^{{*2}} = {d['mu2']}$ (Fig. {d['fig_h112']})", fontsize=11, fontweight='bold')
        ax_bot.set_xlabel(r'$r / \sigma$', fontweight='bold')
        ax_bot.set_ylabel(r'$h^{112}(r)$', fontweight='bold')
        ax_bot.set_xlim(0.9, 4.0)
        ax_bot.set_ylim(d["ylim_h112"])
        ax_bot.grid(True)
        if col_idx == 0:
            ax_bot.legend(frameon=True, loc='upper right')

    plt.suptitle("Angular Correlation Projections $h^{110}(r)$ and $h^{112}(r)$ across Density Regimes\nPatey, Levesque & Weis (1979) vs. Integral Equation Closures", fontsize=13, fontweight='bold', y=0.98)
    plt.tight_layout(rect=[0, 0, 1, 0.95])

    pdf_path = os.path.join(PLOTS_DIR, "fig_plw_angular_grid.pdf")
    png_path = os.path.join(PLOTS_DIR, "fig_plw_angular_grid.png")
    plt.savefig(pdf_path)
    plt.savefig(png_path)
    plt.close()
    print("✓ Generated fig_plw_angular_grid")

def plot_g000_density_progression():
    """Generates 4-panel g000 density progression (Figs 1, 2, 3, 4 of PLW 1979)"""
    fig, axes = plt.subplots(2, 2, figsize=(11, 9), dpi=300)
    
    panels = [
        {"ax": axes[0, 0], "key": "rho_0.15_mu2_2.0", "rho": 0.15, "mu2": 2.0, "mc": "Fig1_rho_0.15_mu2_2.0.dat", "fig_num": "1", "ylim": (0.5, 2.8)},
        {"ax": axes[0, 1], "key": "rho_0.4_mu2_2.75", "rho": 0.40, "mu2": 2.75, "mc": "Fig2_rho_0.4_mu2_2.75.dat", "fig_num": "2", "ylim": (0.5, 4.5)},
        {"ax": axes[1, 0], "key": "rho_0.6_mu2_2.75", "rho": 0.60, "mu2": 2.75, "mc": "Fig3_rho_0.6_mu2_2.75.dat", "fig_num": "3", "ylim": (0.5, 4.5)},
        {"ax": axes[1, 1], "key": "rho_0.8_mu2_2.75", "rho": 0.80, "mu2": 2.75, "mc": "Fig4_rho_0.8_mu2_2.75.dat", "fig_num": "4", "ylim": (0.5, 6.0)},
    ]

    for p in panels:
        ax = p["ax"]
        mc_f = os.path.join(PLW_MC_DIR, p["mc"])
        if os.path.exists(mc_f):
            data_mc = np.loadtxt(mc_f)
            ax.scatter(data_mc[:, 0], data_mc[:, 1], color=C_MC, s=32, facecolors='none', edgecolors='black', linewidth=1.4, zorder=5, label='MC (PLW 1979)')

        sol_msa = load_sol(p["key"], "MSA")
        sol_lhnc = load_sol(p["key"], "LHNC")
        sol_qhnc = load_sol(p["key"], "QHNC")
        sol_rhnc = load_sol(p["key"], "RHNC")

        if sol_msa:
            mask = (sol_msa["r"] >= 1.0) & (sol_msa["r"] <= 3.5)
            ax.plot(sol_msa["r"][mask], sol_msa["g000"][mask], color=C_MSA, linestyle='--', label='MSA')
        if sol_lhnc:
            mask = (sol_lhnc["r"] >= 1.0) & (sol_lhnc["r"] <= 3.5)
            ax.plot(sol_lhnc["r"][mask], sol_lhnc["g000"][mask], color=C_LHNC, linestyle='-.', label='LHNC')
        if sol_qhnc:
            mask = (sol_qhnc["r"] >= 1.0) & (sol_qhnc["r"] <= 3.5)
            ax.plot(sol_qhnc["r"][mask], sol_qhnc["g000"][mask], color=C_QHNC, linestyle='-', linewidth=2.0, label='QHNC')
        if sol_rhnc:
            mask = (sol_rhnc["r"] >= 1.0) & (sol_rhnc["r"] <= 3.5)
            ax.plot(sol_rhnc["r"][mask], sol_rhnc["g000"][mask], color=C_RHNC, linestyle=':', linewidth=2.2, label='RHNC')

        ax.set_title(f"Panel ({chr(97 + panels.index(p))}): $\\rho^* = {p['rho']}$, $\\mu^{{*2}} = {p['mu2']}$ (PLW Fig. {p['fig_num']})", fontsize=11, fontweight='bold')
        ax.set_xlabel(r'$r / \sigma$', fontweight='bold')
        ax.set_ylabel(r'$g^{000}(r)$', fontweight='bold')
        ax.set_xlim(0.9, 3.5)
        ax.set_ylim(p["ylim"])
        ax.grid(True)
        ax.legend(frameon=True, loc='upper right')

    plt.suptitle("Radial Distribution Function $g^{000}(r)$ Across Density Regimes $\\rho^* \\in [0.15, 0.80]$\nPatey, Levesque & Weis (1979) vs. Integral Equation Closures", fontsize=13, fontweight='bold', y=0.98)
    plt.tight_layout(rect=[0, 0, 1, 0.95])

    pdf_path = os.path.join(PLOTS_DIR, "fig_plw_g000_grid.pdf")
    png_path = os.path.join(PLOTS_DIR, "fig_plw_g000_grid.png")
    plt.savefig(pdf_path)
    plt.savefig(png_path)
    plt.close()
    print("✓ Generated fig_plw_g000_grid")

def plot_plw_error_bars():
    """Generates bar chart comparing RMSE across all 10 PLW datasets for MSA, LHNC, QHNC, RHNC"""
    summary_file = os.path.join(DATA_DIR, "plw_benchmark_summary.dat")
    if not os.path.exists(summary_file):
        return

    # Parse summary
    data = []
    with open(summary_file) as f:
        for line in f:
            if line.startswith("#"):
                continue
            parts = line.strip().split()
            if len(parts) >= 8:
                data.append({
                    "state": parts[0],
                    "closure": parts[1],
                    "g000": float(parts[5]) if parts[5] != "---" else np.nan,
                    "h110": float(parts[6]) if parts[6] != "---" else np.nan,
                    "h112": float(parts[7]) if parts[7] != "---" else np.nan
                })

    states = ["rho_0.15_mu2_2.0", "rho_0.4_mu2_2.75", "rho_0.6_mu2_2.75", "rho_0.8_mu2_2.75"]
    state_labels = [r"$\rho^*=0.15$", r"$\rho^*=0.40$", r"$\rho^*=0.60$", r"$\rho^*=0.80$"]
    closures = ["MSA", "LHNC", "QHNC", "RHNC"]
    colors = [C_MSA, C_LHNC, C_QHNC, C_RHNC]

    fig, axes = plt.subplots(1, 3, figsize=(14, 4.8), dpi=300)
    n_cl = len(closures)
    width = 0.8 / n_cl

    # Subplot 1: RMSE g000
    ax = axes[0]
    x = np.arange(len(states))

    for i, cl in enumerate(closures):
        vals = [next((d["g000"] for d in data if d["state"] == s and d["closure"] == cl), np.nan) for s in states]
        pos = x + (i - (n_cl - 1) / 2.0) * width
        ax.bar(pos, vals, width, label=cl, color=colors[i], alpha=0.85, edgecolor='black', linewidth=0.8)

    ax.set_ylabel(r'$\text{RMSE}(g^{000})$', fontweight='bold')
    ax.set_xticks(x)
    ax.set_xticklabels(state_labels)
    ax.set_title("Radial Distribution $g^{000}(r)$ Error", fontweight='bold')
    ax.grid(True, axis='y')
    ax.legend()

    # Subplot 2: RMSE h110
    ax = axes[1]
    x_sub = np.arange(3)
    for i, cl in enumerate(closures):
        vals = [next((d["h110"] for d in data if d["state"] == s and d["closure"] == cl), np.nan) for s in states[:3]]
        pos = x_sub + (i - (n_cl - 1) / 2.0) * width
        ax.bar(pos, vals, width, label=cl, color=colors[i], alpha=0.85, edgecolor='black', linewidth=0.8)

    ax.set_ylabel(r'$\text{RMSE}(h^{110})$', fontweight='bold')
    ax.set_xticks(x_sub)
    ax.set_xticklabels(state_labels[:3])
    ax.set_title("Dipolar Projection $h^{110}(r)$ Error", fontweight='bold')
    ax.grid(True, axis='y')

    # Subplot 3: RMSE h112
    ax = axes[2]
    for i, cl in enumerate(closures):
        vals = [next((d["h112"] for d in data if d["state"] == s and d["closure"] == cl), np.nan) for s in states[:3]]
        pos = x_sub + (i - (n_cl - 1) / 2.0) * width
        ax.bar(pos, vals, width, label=cl, color=colors[i], alpha=0.85, edgecolor='black', linewidth=0.8)

    ax.set_ylabel(r'$\text{RMSE}(h^{112})$', fontweight='bold')
    ax.set_xticks(x_sub)
    ax.set_xticklabels(state_labels[:3])
    ax.set_title("Anisotropic Projection $h^{112}(r)$ Error", fontweight='bold')
    ax.grid(True, axis='y')

    plt.suptitle("Root Mean Square Error Comparison across Densities (Patey et al. 1979)", fontsize=13, fontweight='bold', y=0.98)
    plt.tight_layout(rect=[0, 0, 1, 0.95])

    pdf_path = os.path.join(PLOTS_DIR, "fig_plw_error_comparison.pdf")
    png_path = os.path.join(PLOTS_DIR, "fig_plw_error_comparison.png")
    plt.savefig(pdf_path)
    plt.savefig(png_path)
    plt.close()
    print("✓ Generated fig_plw_error_comparison")

def main():
    print("Generating Publication Plots for Patey, Levesque & Weis (1979)...")
    
    # Individual g000 figures
    plot_single_g000("rho_0.15_mu2_2.0", "Fig1_rho_0.15_mu2_2.0.dat", "1", 0.15, 2.0, "fig_plw_1_g000_rho015", ylim=(0.5, 2.8))
    plot_single_g000("rho_0.4_mu2_2.75", "Fig2_rho_0.4_mu2_2.75.dat", "2", 0.40, 2.75, "fig_plw_2_g000_rho04", ylim=(0.5, 4.5))
    plot_single_g000("rho_0.6_mu2_2.75", "Fig3_rho_0.6_mu2_2.75.dat", "3", 0.60, 2.75, "fig_plw_3_g000_rho06", ylim=(0.5, 4.5))
    plot_single_g000("rho_0.8_mu2_2.75", "Fig4_rho_0.8_mu2_2.75.dat", "4", 0.80, 2.75, "fig_plw_4_g000_rho08", ylim=(0.5, 6.0))

    # Grid progression figures
    plot_g000_density_progression()
    plot_angular_projections_plw()
    plot_plw_error_bars()

    print("All PLW plots successfully generated in reports/monte_carlo_benchmark/plots/")

if __name__ == "__main__":
    main()
