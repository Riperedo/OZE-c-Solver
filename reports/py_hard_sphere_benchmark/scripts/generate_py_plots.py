#!/usr/bin/env python3
"""
Publication Figure Generator for Percus-Yevick Hard Sphere Benchmark
Generates publication-quality PDF and high-res PNG plots comparing
numerical OZE_c_solver results with exact Wertheim-Thiele analytical datasets.
Includes full g(r) and S(k) comparison figures (numerical vs. analytical).
"""

import os
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import AutoMinorLocator
from matplotlib.lines import Line2D

# Styling configuration
plt.rcParams.update({
    "font.family": "serif",
    "font.size": 11,
    "axes.labelsize": 12,
    "axes.titlesize": 12,
    "xtick.labelsize": 10,
    "ytick.labelsize": 10,
    "legend.fontsize": 9.5,
    "figure.titlesize": 13,
    "lines.linewidth": 1.8,
    "lines.markersize": 5,
    "grid.alpha": 0.3,
    "grid.linestyle": "--",
    "savefig.dpi": 300,
    "savefig.bbox": "tight"
})

REF_DIR = "test/data_gdr_PY_analitica"
DATA_DIR = "reports/py_hard_sphere_benchmark/data"
PLOTS_DIR = "reports/py_hard_sphere_benchmark/plots"
os.makedirs(PLOTS_DIR, exist_ok=True)

PHI_LIST = [0.10, 0.20, 0.30, 0.40, 0.50, 0.55, 0.60]
COLORS = ['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', '#9467bd', '#8c564b', '#e377c2']

def py_rho_ck_analytical(k, phi):
    """Exact Fourier transform rho * c_hat(k) for hard spheres with Percus-Yevick closure."""
    x = np.asarray(k, dtype=float)
    eta = phi
    l1 = (1.0 + 2.0 * eta) ** 2 / (1.0 - eta) ** 4
    l2 = -(1.0 + eta / 2.0) ** 2 / (1.0 - eta) ** 4
    
    small = np.abs(x) < 0.1
    large = ~small
    
    term1 = np.zeros_like(x)
    term2 = np.zeros_like(x)
    term3 = np.zeros_like(x)
    
    if np.any(large):
        xl = x[large]
        term1[large] = l1 * (np.sin(xl) - xl * np.cos(xl)) / (xl ** 3)
        term2[large] = 6.0 * eta * l2 * (2.0 * xl * np.sin(xl) - (xl ** 2 - 2.0) * np.cos(xl) - 2.0) / (xl ** 4)
        term3[large] = 0.5 * eta * l1 * ((4.0 * xl ** 3 - 24.0 * xl) * np.sin(xl) - (xl ** 4 - 12.0 * xl ** 2 + 24.0) * np.cos(xl) + 24.0) / (xl ** 6)
        
    if np.any(small):
        xs = x[small]
        xs2 = xs ** 2
        xs4 = xs ** 4
        xs6 = xs ** 6
        t1_s = 1.0 / 3.0 - xs2 / 30.0 + xs4 / 840.0 - xs6 / 45360.0
        t2_s = 1.0 / 4.0 - xs2 / 24.0 + xs4 / 720.0 - xs6 / 40320.0
        t3_s = 1.0 / 6.0 - xs2 / 48.0 + xs4 / 1200.0 - xs6 / 50400.0
        
        term1[small] = l1 * t1_s
        term2[small] = 6.0 * eta * l2 * t2_s
        term3[small] = 0.5 * eta * l1 * t3_s
        
    rho_c = -24.0 * eta * (term1 + term2 + term3)
    return rho_c

def py_sk_analytical(k, phi):
    """Exact static structure factor S(k) = 1 / (1 - rho * c_hat(k))."""
    rho_c = py_rho_ck_analytical(k, phi)
    return 1.0 / (1.0 - rho_c)

def plot_fig1_gr_all_phi():
    """Figure 1: Full radial distribution function g(r) vs r/sigma across all 7 packing fractions."""
    fig, ax = plt.subplots(figsize=(8.5, 5.8))
    
    handles = []
    handles.append(Line2D([0], [0], color='black', lw=1.8, label='Solid: Numerical (OZE Solver)'))
    handles.append(Line2D([0], [0], marker='o', color='black', markerfacecolor='none', 
                          markeredgewidth=1.2, lw=0, markersize=5, label='Circles: Analytical (Wertheim-Thiele)'))
    
    for i, phi in enumerate(PHI_LIST):
        phi_str = f"{phi:.2f}"
        ref_file = os.path.join(REF_DIR, f"gr_analitica_phi{phi_str}.dat")
        sol_file = os.path.join(DATA_DIR, f"solution_PY_phi_{phi_str}.dat")
        
        ref_data = np.loadtxt(ref_file)
        sol_data = np.loadtxt(sol_file)
        
        r_ref, g_ref = ref_data[:, 0], ref_data[:, 1]
        r_sol, g_sol = sol_data[:, 0], sol_data[:, 1]
        
        m_sol = (r_sol >= 0.0) & (r_sol <= 4.0)
        m_ref = (r_ref >= 0.0) & (r_ref <= 4.0)
        
        line, = ax.plot(r_sol[m_sol], g_sol[m_sol], color=COLORS[i], label=f"$\\phi = {phi:.2f}$")
        ax.plot(r_ref[m_ref][::2], g_ref[m_ref][::2], 'o', color=COLORS[i], 
                markerfacecolor='none', markeredgewidth=1.2, markersize=4.5)
        handles.append(line)

    ax.legend(handles=handles, loc="upper right", frameon=True, ncol=2, fontsize=9.0)
    ax.set_xlabel("$r / \\sigma$")
    ax.set_ylabel("$g(r)$")
    ax.set_xlim(0.5, 4.0)
    ax.set_ylim(-0.2, 8.8)
    ax.grid(True)
    ax.xaxis.set_minor_locator(AutoMinorLocator())
    ax.yaxis.set_minor_locator(AutoMinorLocator())
    ax.set_title("Percus-Yevick Hard Sphere Radial Distribution Function $g(r)$")
    
    fig.savefig(os.path.join(PLOTS_DIR, "fig1_gr_all_phi.pdf"))
    fig.savefig(os.path.join(PLOTS_DIR, "fig1_gr_all_phi.png"))
    plt.close(fig)
    print("Generated Figure 1: g(r) across all packing fractions.")

def plot_fig2_residuals_delta_gr():
    """Figure 2: Residual difference Delta g(r) = g_num(r) - g_ana(r) for all phi."""
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(8.5, 7.5), sharex=True, gridspec_kw={'height_ratios': [2, 1.2]})
    
    for i, phi in enumerate(PHI_LIST):
        phi_str = f"{phi:.2f}"
        ref_file = os.path.join(REF_DIR, f"gr_analitica_phi{phi_str}.dat")
        sol_file = os.path.join(DATA_DIR, f"solution_PY_phi_{phi_str}.dat")
        
        ref_data = np.loadtxt(ref_file)
        sol_data = np.loadtxt(sol_file)
        
        r_ref, g_ref = ref_data[:, 0], ref_data[:, 1]
        r_sol, g_sol = sol_data[:, 0], sol_data[:, 1]
        
        mask = (r_ref >= 1.0) & (r_ref <= 4.0)
        r_eval = r_ref[mask]
        g_ref_eval = g_ref[mask]
        g_sol_eval = np.interp(r_eval, r_sol, g_sol)
        
        diff = g_sol_eval - g_ref_eval
        
        ax1.plot(r_eval, g_sol_eval, color=COLORS[i], label=f"$\\phi = {phi:.2f}$")
        ax1.plot(r_eval[::3], g_ref_eval[::3], 'o', color=COLORS[i], markerfacecolor='none', markersize=4)
        ax2.plot(r_eval, diff, color=COLORS[i], label=f"$\\phi = {phi:.2f}$")

    ax1.set_ylabel("$g(r)$")
    ax1.set_xlim(1.0, 4.0)
    ax1.set_ylim(0.5, 8.8)
    ax1.grid(True)
    ax1.set_title("Fluid Region $g(r)$ and Residual Deviations $\\Delta g(r) = g_{\\text{num}}(r) - g_{\\text{ana}}(r)$")
    ax1.legend(loc="upper right", ncol=2)
    
    ax2.axhline(0, color='gray', linestyle=':', lw=1.2)
    ax2.set_xlabel("$r / \\sigma$")
    ax2.set_ylabel("$\\Delta g(r)$")
    ax2.set_xlim(1.0, 4.0)
    ax2.grid(True)
    ax2.xaxis.set_minor_locator(AutoMinorLocator())
    ax2.yaxis.set_minor_locator(AutoMinorLocator())
    
    fig.tight_layout()
    fig.savefig(os.path.join(PLOTS_DIR, "fig2_residuals_delta_gr.pdf"))
    fig.savefig(os.path.join(PLOTS_DIR, "fig2_residuals_delta_gr.png"))
    plt.close(fig)
    print("Generated Figure 2: Residual difference Delta g(r).")

def plot_fig3_contact_value_scaling():
    """Figure 3: Contact value g(sigma^+) vs packing fraction phi and relative error."""
    summary_file = os.path.join(DATA_DIR, "py_benchmark_summary.dat")
    data = np.loadtxt(summary_file)
    
    phi = data[:, 0]
    g_contact_num = data[:, 3]
    g_contact_ana = data[:, 4]
    contact_err_pct = data[:, 5]
    
    phi_dense = np.linspace(0.05, 0.62, 200)
    g_contact_dense = (1.0 + phi_dense / 2.0) / ((1.0 - phi_dense) ** 2)
    
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(7.5, 7.0), sharex=True, gridspec_kw={'height_ratios': [2, 1]})
    
    ax1.plot(phi_dense, g_contact_dense, 'k--', lw=1.8, label="Wertheim-Thiele Theory: $g(\\sigma^+) = \\frac{1 + \\eta/2}{(1-\\eta)^2}$")
    ax1.plot(phi, g_contact_num, 's-', color='#1f77b4', lw=1.5, markerfacecolor='white', 
             markeredgewidth=1.8, markersize=7, label="Numerical Extrapolation ($N=4096$)")
    ax1.plot(phi, g_contact_ana, 'ro', markersize=6, label="Analytical Reference Datasets")
    
    ax1.set_ylabel("Contact Value $g(\\sigma^+)$")
    ax1.set_ylim(1.0, 9.0)
    ax1.grid(True)
    ax1.legend(loc="upper left")
    ax1.set_title("Percus-Yevick Contact Value $g(\\sigma^+)$ Scaling vs. Packing Fraction $\\phi$")
    
    ax2.plot(phi, contact_err_pct, 'd-', color='#d62728', lw=1.5, markerfacecolor='white', markeredgewidth=1.5, markersize=6)
    ax2.axhline(0, color='gray', linestyle=':', lw=1)
    ax2.set_xlabel("Packing Fraction $\\phi$")
    ax2.set_ylabel("Relative Error (\\%)")
    ax2.set_xlim(0.08, 0.62)
    ax2.set_ylim(0.0, 6.0)
    ax2.grid(True)
    ax2.xaxis.set_minor_locator(AutoMinorLocator())
    ax2.yaxis.set_minor_locator(AutoMinorLocator())
    
    fig.tight_layout()
    fig.savefig(os.path.join(PLOTS_DIR, "fig3_contact_value_scaling.pdf"))
    fig.savefig(os.path.join(PLOTS_DIR, "fig3_contact_value_scaling.png"))
    plt.close(fig)
    print("Generated Figure 3: Contact value scaling.")

def plot_fig4_error_metrics_summary():
    """Figure 4: Global error metrics (RMSE for g(r) & S(k), and Max Error) across packing fraction phi."""
    summary_file = os.path.join(DATA_DIR, "py_benchmark_summary.dat")
    data = np.loadtxt(summary_file)
    
    phi = data[:, 0]
    rmse_gr = data[:, 1]
    rmse_sk = data[:, 6]
    max_err_sk = data[:, 7]
    runtime = data[:, 10]
    
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(10, 4.5))
    
    ax1.semilogy(phi, rmse_gr, 'o-', color='#1f77b4', lw=1.8, markersize=6, label="RMSE $g(r)$")
    ax1.semilogy(phi, rmse_sk, 's--', color='#2ca02c', lw=1.8, markersize=6, label="RMSE $S(k)$")
    ax1.semilogy(phi, max_err_sk, '^:', color='#d62728', lw=1.8, markersize=6, label="Max $\\Delta S(k)$")
    ax1.set_xlabel("Packing Fraction $\\phi$")
    ax1.set_ylabel("Error Magnitude (log scale)")
    ax1.set_title("Numerical Error Metrics vs. $\\phi$")
    ax1.grid(True, which="both")
    ax1.legend(loc="upper left")
    
    ax2.plot(phi, runtime, '^-', color='#9467bd', lw=1.8, markersize=7)
    ax2.set_xlabel("Packing Fraction $\\phi$")
    ax2.set_ylabel("Wall-Clock Execution Time (s)")
    ax2.set_title("Solver Computational Performance")
    ax2.grid(True)
    ax2.set_ylim(0, max(runtime)*1.15)
    
    fig.tight_layout()
    fig.savefig(os.path.join(PLOTS_DIR, "fig4_error_metrics_summary.pdf"))
    fig.savefig(os.path.join(PLOTS_DIR, "fig4_error_metrics_summary.png"))
    plt.close(fig)
    print("Generated Figure 4: Error metrics and runtime scaling.")

def plot_fig5_structure_factor_sk():
    """Figure 5: Static structure factors S(k) - BOTH Numerical and Exact Analytical Solutions with Residuals."""
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(8.5, 8.0), sharex=True, gridspec_kw={'height_ratios': [2.2, 1.0]})
    
    handles = []
    handles.append(Line2D([0], [0], color='black', lw=1.8, label='Solid: Numerical Solver'))
    handles.append(Line2D([0], [0], color='black', lw=1.4, linestyle='--', label='Dashed: Exact Analytical Wertheim'))
    
    for i, phi in enumerate(PHI_LIST):
        phi_str = f"{phi:.2f}"
        sk_file = os.path.join(DATA_DIR, f"sk_PY_phi_{phi_str}.dat")
        if not os.path.exists(sk_file):
            continue
            
        sk_data = np.loadtxt(sk_file)
        k = sk_data[:, 0]
        sk_num = sk_data[:, 1]
        
        m = (k >= 0.0) & (k <= 25.0)
        k_eval = k[m]
        sk_num_eval = sk_num[m]
        sk_ana_eval = py_sk_analytical(k_eval, phi)
        
        diff_sk = sk_num_eval - sk_ana_eval
        
        # Top panel: Numerical (solid) vs Exact Analytical (dashed)
        line, = ax1.plot(k_eval, sk_num_eval, color=COLORS[i], lw=1.8, label=f"$\\phi = {phi:.2f}$")
        ax1.plot(k_eval, sk_ana_eval, color=COLORS[i], lw=1.4, linestyle='--')
        handles.append(line)
        
        # Bottom panel: Residual Delta S(k)
        ax2.plot(k_eval, diff_sk, color=COLORS[i], lw=1.4, label=f"$\\phi = {phi:.2f}$")

    ax1.set_ylabel("Structure Factor $S(k)$")
    ax1.set_xlim(0, 25.0)
    ax1.set_ylim(0, 7.0)
    ax1.grid(True)
    ax1.xaxis.set_minor_locator(AutoMinorLocator())
    ax1.yaxis.set_minor_locator(AutoMinorLocator())
    ax1.legend(handles=handles, loc="upper right", ncol=2, fontsize=8.8)
    ax1.set_title("Percus-Yevick Hard Sphere Static Structure Factor: Numerical vs. Exact Analytical")
    
    ax2.axhline(0, color='gray', linestyle=':', lw=1.2)
    ax2.set_xlabel("Wavevector $k\\sigma$")
    ax2.set_ylabel("$\\Delta S(k) = S_{\\text{num}} - S_{\\text{ana}}$")
    ax2.set_xlim(0, 25.0)
    ax2.set_ylim(-0.25, 0.25)
    ax2.grid(True)
    ax2.xaxis.set_minor_locator(AutoMinorLocator())
    ax2.yaxis.set_minor_locator(AutoMinorLocator())
    
    fig.tight_layout()
    fig.savefig(os.path.join(PLOTS_DIR, "fig5_structure_factor_sk.pdf"))
    fig.savefig(os.path.join(PLOTS_DIR, "fig5_structure_factor_sk.png"))
    plt.close(fig)
    print("Generated Figure 5: Numerical vs. Analytical S(k) with Residuals.")

if __name__ == "__main__":
    plot_fig1_gr_all_phi()
    plot_fig2_residuals_delta_gr()
    plot_fig3_contact_value_scaling()
    plot_fig4_error_metrics_summary()
    plot_fig5_structure_factor_sk()
    print("All Percus-Yevick benchmark figures regenerated successfully.")
