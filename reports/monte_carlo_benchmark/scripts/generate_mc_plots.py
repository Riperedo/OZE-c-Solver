#!/usr/bin/env python3
"""
Publication Figure Generator for Monte Carlo Benchmark Report.
Generates Figures 1 to 6 comparing MSA, LHNC, and QHNC against Monte Carlo simulations
(Fries & Patey, J. Chem. Phys. 82, 429, 1985).
"""

import os
import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d

# Style configuration
plt.rcParams.update({
    'font.family': 'serif',
    'font.size': 11,
    'axes.labelsize': 13,
    'axes.titlesize': 13,
    'xtick.labelsize': 11,
    'ytick.labelsize': 11,
    'legend.fontsize': 10,
    'lines.linewidth': 2.0,
    'lines.markersize': 6.0,
    'figure.dpi': 300,
    'savefig.dpi': 300,
    'savefig.bbox': 'tight'
})

BASE_DIR = os.path.abspath(os.path.join(os.path.dirname(__file__), "../../.."))
MC_DIR = os.path.join(BASE_DIR, "test/data_Monte_Carlo_Simulations/FriesPatey_JCP_824291985")
DATA_DIR = os.path.join(BASE_DIR, "reports/monte_carlo_benchmark/data")
PLOTS_DIR = os.path.join(BASE_DIR, "reports/monte_carlo_benchmark/plots")
os.makedirs(PLOTS_DIR, exist_ok=True)

# Color and style palette
STYLES = {
    "MC": {"color": "black", "marker": "o", "ls": "none", "label": "Monte Carlo (Fries & Patey 1985)"},
    "MSA": {"color": "#1f77b4", "ls": "--", "label": "MSA"},
    "LHNC": {"color": "#2ca02c", "ls": "-.", "label": "LHNC"},
    "QHNC": {"color": "#d62728", "ls": "-", "label": "QHNC"}
}

def load_theory_data(state_tag, closure):
    path = os.path.join(DATA_DIR, f"solution_{state_tag}_{closure}.dat")
    if not os.path.exists(path):
        # Fallback to test/ directory if available
        fb_path = os.path.join(MC_DIR, f"hdr_{state_tag}_{closure}.dat")
        if os.path.exists(fb_path):
            path = fb_path
    if not os.path.exists(path):
        return None
    data = np.loadtxt(path)
    return {
        "r": data[:,0],
        "h000": data[:,1],
        "h110": data[:,2],
        "h112": data[:,3],
        "g000": data[:,1] + 1.0
    }

def plot_fig1():
    """Fig 1: g000(r) at rho*=0.8, mu*2=2.75 (Contact and Medium range)."""
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5))
    
    # Left: Contact region
    mc_1a = np.loadtxt(os.path.join(MC_DIR, "fig1a_rho_0.8_mu2_2.75.dat"))
    ax1.plot(mc_1a[:,0], mc_1a[:,1], STYLES["MC"]["marker"], color=STYLES["MC"]["color"], 
             label=STYLES["MC"]["label"], markeredgecolor='black', markerfacecolor='black', markersize=6)
    
    # Right: Medium range
    mc_1b = np.loadtxt(os.path.join(MC_DIR, "fig1b_rho_0.8_mu2_2.75.dat"))
    ax2.plot(mc_1b[:,0], mc_1b[:,1], STYLES["MC"]["marker"], color=STYLES["MC"]["color"], 
             label=STYLES["MC"]["label"], markeredgecolor='black', markerfacecolor='black', markersize=6)
    
    for cl in ["MSA", "LHNC", "QHNC"]:
        th = load_theory_data("rho_0.8_mu2_2.75", cl)
        if th is not None:
            mask = th["r"] >= 1.0
            ax1.plot(th["r"][mask], th["g000"][mask], ls=STYLES[cl]["ls"], color=STYLES[cl]["color"], label=STYLES[cl]["label"])
            ax2.plot(th["r"][mask], th["g000"][mask], ls=STYLES[cl]["ls"], color=STYLES[cl]["color"], label=STYLES[cl]["label"])
            
    ax1.set_xlim(1.0, 1.6)
    ax1.set_ylim(0.0, 5.5)
    ax1.set_xlabel(r"$r / \sigma$")
    ax1.set_ylabel(r"$g^{000}(r)$")
    ax1.set_title(r"(a) Contact Region ($r \in [1.0, 1.6]\sigma$)")
    ax1.grid(True, linestyle=":", alpha=0.6)
    ax1.legend(loc="upper right", framealpha=0.9)

    ax2.set_xlim(1.0, 4.0)
    ax2.set_ylim(0.5, 1.6)
    ax2.set_xlabel(r"$r / \sigma$")
    ax2.set_ylabel(r"$g^{000}(r)$")
    ax2.set_title(r"(b) Medium-Range Structure ($r \in [1.0, 4.0]\sigma$)")
    ax2.grid(True, linestyle=":", alpha=0.6)
    ax2.legend(loc="upper right", framealpha=0.9)

    plt.suptitle(r"Radial Distribution Function $g^{000}(r)$ at $\rho^* = 0.8, \mu^{*2} = 2.75$ ($T^* = 0.3636$)", fontsize=14, y=1.02)
    plt.tight_layout()
    plt.savefig(os.path.join(PLOTS_DIR, "fig1_g000_mu2_2.75.pdf"))
    plt.savefig(os.path.join(PLOTS_DIR, "fig1_g000_mu2_2.75.png"))
    plt.close()
    print("✓ Generated Fig 1: fig1_g000_mu2_2.75")

def plot_fig2():
    """Fig 2: g000(r) at rho*=0.8, mu*2=2.0 (Contact and Medium range)."""
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5))
    
    # Left: Contact region
    mc_2a = np.loadtxt(os.path.join(MC_DIR, "fig2a_rho_0.8_mu2_2.0.dat"))
    ax1.plot(mc_2a[:,0], mc_2a[:,1], STYLES["MC"]["marker"], color=STYLES["MC"]["color"], 
             label=STYLES["MC"]["label"], markeredgecolor='black', markerfacecolor='black', markersize=6)
    
    # Right: Medium range
    mc_2b = np.loadtxt(os.path.join(MC_DIR, "fig2b_rho_0.8_mu2_2.0.dat"))
    ax2.plot(mc_2b[:,0], mc_2b[:,1], STYLES["MC"]["marker"], color=STYLES["MC"]["color"], 
             label=STYLES["MC"]["label"], markeredgecolor='black', markerfacecolor='black', markersize=6)
    
    for cl in ["MSA", "LHNC", "QHNC"]:
        th = load_theory_data("rho_0.8_mu2_2.0", cl)
        if th is not None:
            mask = th["r"] >= 1.0
            ax1.plot(th["r"][mask], th["g000"][mask], ls=STYLES[cl]["ls"], color=STYLES[cl]["color"], label=STYLES[cl]["label"])
            ax2.plot(th["r"][mask], th["g000"][mask], ls=STYLES[cl]["ls"], color=STYLES[cl]["color"], label=STYLES[cl]["label"])
            
    ax1.set_xlim(1.0, 1.6)
    ax1.set_ylim(0.0, 5.0)
    ax1.set_xlabel(r"$r / \sigma$")
    ax1.set_ylabel(r"$g^{000}(r)$")
    ax1.set_title(r"(a) Contact Region ($r \in [1.0, 1.6]\sigma$)")
    ax1.grid(True, linestyle=":", alpha=0.6)
    ax1.legend(loc="upper right", framealpha=0.9)

    ax2.set_xlim(1.0, 4.0)
    ax2.set_ylim(0.5, 1.6)
    ax2.set_xlabel(r"$r / \sigma$")
    ax2.set_ylabel(r"$g^{000}(r)$")
    ax2.set_title(r"(b) Medium-Range Structure ($r \in [1.0, 4.0]\sigma$)")
    ax2.grid(True, linestyle=":", alpha=0.6)
    ax2.legend(loc="upper right", framealpha=0.9)

    plt.suptitle(r"Radial Distribution Function $g^{000}(r)$ at $\rho^* = 0.8, \mu^{*2} = 2.00$ ($T^* = 0.5000$)", fontsize=14, y=1.02)
    plt.tight_layout()
    plt.savefig(os.path.join(PLOTS_DIR, "fig2_g000_mu2_2.0.pdf"))
    plt.savefig(os.path.join(PLOTS_DIR, "fig2_g000_mu2_2.0.png"))
    plt.close()
    print("✓ Generated Fig 2: fig2_g000_mu2_2.0")

def plot_fig3():
    """Fig 3: Dipolar Projection h110(r) at rho*=0.8, mu*2=2.0."""
    plt.figure(figsize=(7.5, 5.5))
    mc_3 = np.loadtxt(os.path.join(MC_DIR, "fig3_rho_0.8_mu2_2.0.dat"))
    plt.plot(mc_3[:,0], mc_3[:,1], STYLES["MC"]["marker"], color=STYLES["MC"]["color"], 
             label=STYLES["MC"]["label"], markeredgecolor='black', markerfacecolor='black', markersize=6)
    
    for cl in ["MSA", "LHNC", "QHNC"]:
        th = load_theory_data("rho_0.8_mu2_2.0", cl)
        if th is not None:
            mask = th["r"] >= 1.0
            plt.plot(th["r"][mask], th["h110"][mask], ls=STYLES[cl]["ls"], color=STYLES[cl]["color"], label=STYLES[cl]["label"])

    plt.xlim(1.0, 4.0)
    plt.ylim(-0.5, 3.5)
    plt.xlabel(r"$r / \sigma$")
    plt.ylabel(r"$h^{110}(r)$")
    plt.title(r"Dipolar Angular Projection $h^{110}(r)$ at $\rho^* = 0.8, \mu^{*2} = 2.00$ ($T^* = 0.5000$)")
    plt.grid(True, linestyle=":", alpha=0.6)
    plt.legend(loc="upper right", framealpha=0.9)
    plt.tight_layout()
    plt.savefig(os.path.join(PLOTS_DIR, "fig3_h110_mu2_2.0.pdf"))
    plt.savefig(os.path.join(PLOTS_DIR, "fig3_h110_mu2_2.0.png"))
    plt.close()
    print("✓ Generated Fig 3: fig3_h110_mu2_2.0")

def plot_fig4():
    """Fig 4: Anisotropic Projection h112(r) at rho*=0.8, mu*2=2.0."""
    plt.figure(figsize=(7.5, 5.5))
    mc_4 = np.loadtxt(os.path.join(MC_DIR, "fig4_rho_0.8_mu2_2.0.dat"))
    plt.plot(mc_4[:,0], mc_4[:,1], STYLES["MC"]["marker"], color=STYLES["MC"]["color"], 
             label=STYLES["MC"]["label"], markeredgecolor='black', markerfacecolor='black', markersize=6)
    
    for cl in ["MSA", "LHNC", "QHNC"]:
        th = load_theory_data("rho_0.8_mu2_2.0", cl)
        if th is not None:
            mask = th["r"] >= 1.0
            plt.plot(th["r"][mask], th["h112"][mask], ls=STYLES[cl]["ls"], color=STYLES[cl]["color"], label=STYLES[cl]["label"])

    plt.xlim(1.0, 4.0)
    plt.ylim(-0.5, 4.5)
    plt.xlabel(r"$r / \sigma$")
    plt.ylabel(r"$h^{112}(r)$")
    plt.title(r"Anisotropic Projection $h^{112}(r)$ at $\rho^* = 0.8, \mu^{*2} = 2.00$ ($T^* = 0.5000$)")
    plt.grid(True, linestyle=":", alpha=0.6)
    plt.legend(loc="upper right", framealpha=0.9)
    plt.tight_layout()
    plt.savefig(os.path.join(PLOTS_DIR, "fig4_h112_mu2_2.0.pdf"))
    plt.savefig(os.path.join(PLOTS_DIR, "fig4_h112_mu2_2.0.png"))
    plt.close()
    print("✓ Generated Fig 4: fig4_h112_mu2_2.0")

def plot_fig5():
    """Fig 5: Comprehensive Error Metric Bar Chart across Closures."""
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(13, 5.5))
    
    closures = ["MSA", "LHNC", "QHNC"]
    colors = [STYLES[c]["color"] for c in closures]
    
    # State 1: mu2 = 2.75
    rmse_1_contact = [0.8121, 0.0987, 0.0576]
    rmse_1_med = [0.0367, 0.0246, 0.0238]
    
    x = np.arange(2)
    width = 0.25
    
    ax1.bar(x - width, [rmse_1_contact[0], rmse_1_med[0]], width, label="MSA", color=colors[0], alpha=0.85, edgecolor='black')
    ax1.bar(x, [rmse_1_contact[1], rmse_1_med[1]], width, label="LHNC", color=colors[1], alpha=0.85, edgecolor='black')
    ax1.bar(x + width, [rmse_1_contact[2], rmse_1_med[2]], width, label="QHNC", color=colors[2], alpha=0.85, edgecolor='black')
    
    ax1.set_xticks(x)
    ax1.set_xticklabels([r"$g^{000}_{\mathrm{contact}}$ ($r \leq 1.6\sigma$)", r"$g^{000}_{\mathrm{medium}}$ ($r \leq 4.0\sigma$)"])
    ax1.set_ylabel("RMSE vs. Monte Carlo")
    ax1.set_title(r"(a) State 1: $\rho^* = 0.8, \mu^{*2} = 2.75$ ($T^* = 0.3636$)")
    ax1.set_yscale("log")
    ax1.grid(True, which="both", linestyle=":", alpha=0.6)
    ax1.legend(loc="upper right")

    # State 2: mu2 = 2.0
    rmse_2_contact = [0.5235, 0.0876, 0.0586]
    rmse_2_med = [0.0715, 0.0248, 0.0236]
    rmse_2_110 = [0.3693, 0.1051, 0.1042]
    rmse_2_112 = [0.3800, 0.1709, 0.1731]

    x2 = np.arange(4)
    ax2.bar(x2 - width, [rmse_2_contact[0], rmse_2_med[0], rmse_2_110[0], rmse_2_112[0]], width, label="MSA", color=colors[0], alpha=0.85, edgecolor='black')
    ax2.bar(x2, [rmse_2_contact[1], rmse_2_med[1], rmse_2_110[1], rmse_2_112[1]], width, label="LHNC", color=colors[1], alpha=0.85, edgecolor='black')
    ax2.bar(x2 + width, [rmse_2_contact[2], rmse_2_med[2], rmse_2_110[2], rmse_2_112[2]], width, label="QHNC", color=colors[2], alpha=0.85, edgecolor='black')
    
    ax2.set_xticks(x2)
    ax2.set_xticklabels([r"$g^{000}_{\mathrm{cont}}$", r"$g^{000}_{\mathrm{med}}$", r"$h^{110}$", r"$h^{112}$"])
    ax2.set_ylabel("RMSE vs. Monte Carlo")
    ax2.set_title(r"(b) State 2: $\rho^* = 0.8, \mu^{*2} = 2.00$ ($T^* = 0.5000$)")
    ax2.set_yscale("log")
    ax2.grid(True, which="both", linestyle=":", alpha=0.6)
    ax2.legend(loc="upper right")

    plt.suptitle("Comparative Accuracy of MSA, LHNC, and QHNC against Monte Carlo Benchmarks", fontsize=14, y=1.02)
    plt.tight_layout()
    plt.savefig(os.path.join(PLOTS_DIR, "fig5_error_comparison.pdf"))
    plt.savefig(os.path.join(PLOTS_DIR, "fig5_error_comparison.png"))
    plt.close()
    print("✓ Generated Fig 5: fig5_error_comparison")

def plot_fig6():
    """Fig 6: Convergence Iteration Progression across Phases."""
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5))
    
    phases = ["Phase 1\n(Cold Start)", "Phase 2\n(Continuation)", "Phase 3\n(Warm-Start)"]
    x = np.arange(len(phases))
    width = 0.25

    # State 1: mu2 = 2.75 iters
    iters_1_msa = [65, 19, 22]
    iters_1_lhnc = [157, 67, 1]
    iters_1_qhnc = [210, 89, 110] # 210 represents stalled

    ax1.bar(x - width, iters_1_msa, width, label="MSA", color=STYLES["MSA"]["color"], alpha=0.85, edgecolor='black')
    ax1.bar(x, iters_1_lhnc, width, label="LHNC", color=STYLES["LHNC"]["color"], alpha=0.85, edgecolor='black')
    ax1.bar(x + width, iters_1_qhnc, width, label="QHNC", color=STYLES["QHNC"]["color"], alpha=0.85, edgecolor='black')

    ax1.set_xticks(x)
    ax1.set_xticklabels(phases)
    ax1.set_ylabel("Total Iterations to Convergence")
    ax1.set_title(r"(a) State 1: $\rho^* = 0.8, \mu^{*2} = 2.75$ ($T^* = 0.3636$)")
    ax1.grid(True, linestyle=":", alpha=0.6)
    ax1.legend(loc="upper right")

    # State 2: mu2 = 2.0 iters
    iters_2_msa = [101, 52, 15]
    iters_2_lhnc = [130, 45, 1]
    iters_2_qhnc = [210, 61, 75]

    ax2.bar(x - width, iters_2_msa, width, label="MSA", color=STYLES["MSA"]["color"], alpha=0.85, edgecolor='black')
    ax2.bar(x, iters_2_lhnc, width, label="LHNC", color=STYLES["LHNC"]["color"], alpha=0.85, edgecolor='black')
    ax2.bar(x + width, iters_2_qhnc, width, label="QHNC", color=STYLES["QHNC"]["color"], alpha=0.85, edgecolor='black')

    ax2.set_xticks(x)
    ax2.set_xticklabels(phases)
    ax2.set_ylabel("Total Iterations to Convergence")
    ax2.set_title(r"(b) State 2: $\rho^* = 0.8, \mu^{*2} = 2.00$ ($T^* = 0.5000$)")
    ax2.grid(True, linestyle=":", alpha=0.6)
    ax2.legend(loc="upper right")

    plt.suptitle("Computational Cost and Convergence Robustness across Evaluation Phases", fontsize=14, y=1.02)
    plt.tight_layout()
    plt.savefig(os.path.join(PLOTS_DIR, "fig6_convergence_phases.pdf"))
    plt.savefig(os.path.join(PLOTS_DIR, "fig6_convergence_phases.png"))
    plt.close()
    print("✓ Generated Fig 6: fig6_convergence_phases")

def main():
    print("Generating Publication Plots for Monte Carlo Benchmark Report...")
    plot_fig1()
    plot_fig2()
    plot_fig3()
    plot_fig4()
    plot_fig5()
    plot_fig6()
    print("All plots successfully generated in reports/monte_carlo_benchmark/plots/")

if __name__ == "__main__":
    main()
