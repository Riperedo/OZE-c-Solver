#!/usr/bin/env python3
"""
Percus-Yevick Hard Sphere Benchmark Suite
Evaluates OZE_c_solver against the exact analytical Wertheim-Thiele Percus-Yevick solution
across 7 volume fractions: phi in {0.10, 0.20, 0.30, 0.40, 0.50, 0.55, 0.60}.
Computes both g(r) and S(k) numerical vs exact analytical error metrics.
"""

import os
import sys
import time
import subprocess
import numpy as np

# Packing fractions present in test/data_gdr_PY_analitica
PHI_VALUES = [0.10, 0.20, 0.30, 0.40, 0.50, 0.55, 0.60]
EXECUTABLE = "./build/facdes_solver"
REF_DIR = "test/data_gdr_PY_analitica"
OUT_DATA_DIR = "reports/py_hard_sphere_benchmark/data"
DEFAULT_NODES = 4096
DEFAULT_KNODES = 1024

def get_py_contact_analytical(phi):
    """Percus-Yevick contact value from Wertheim-Thiele theory: g(1+) = (1 + phi/2) / (1 - phi)^2"""
    return (1.0 + phi / 2.0) / ((1.0 - phi) ** 2)

def get_py_sk0_analytical(phi):
    """Percus-Yevick isothermal compressibility S(0) = (1 - phi)^4 / (1 + 2*phi)^2"""
    return ((1.0 - phi) ** 4) / ((1.0 + 2.0 * phi) ** 2)

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

def extrapolate_contact_value(r, g):
    """Extrapolate g(r) to r=1.0+ using the first few fluid grid points (r >= 1.0)"""
    idx_fluid = np.where(r >= 1.0)[0]
    if len(idx_fluid) < 4:
        return 0.0
    idx_fit = idx_fluid[:4]
    p = np.polyfit(r[idx_fit], g[idx_fit], 2)
    return np.polyval(p, 1.0)

def run_benchmark():
    os.makedirs(OUT_DATA_DIR, exist_ok=True)
    
    summary_results = []
    
    print("=" * 105)
    print("      PERCUS-YEVICK HARD SPHERE NUMERICAL BENCHMARK EVALUATION (g(r) and S(k))")
    print("=" * 105)
    print(f"{'phi':>6} | {'g(r) RMSE':>12} | {'g_num(1+)':>10} | {'g_ana(1+)':>10} | {'Cont. Err%':>10} | {'S(k) RMSE':>12} | {'S(k) MaxErr':>12} | {'Time (s)':>8}")
    print("-" * 105)
    
    for phi in PHI_VALUES:
        phi_str = f"{phi:.2f}"
        ref_file = os.path.join(REF_DIR, f"gr_analitica_phi{phi_str}.dat")
        if not os.path.exists(ref_file):
            print(f"Warning: Reference file {ref_file} not found. Skipping.")
            continue
            
        ref_data = np.loadtxt(ref_file)
        r_ref = ref_data[:, 0]
        g_ref = ref_data[:, 1]
        
        t0 = time.time()
        cmd = [
            EXECUTABLE,
            "--closure", "PY",
            "--potential", "7",
            "--volfactor", str(phi),
            "--temp", "1.0",
            "--nodes", str(DEFAULT_NODES),
            "--knodes", str(DEFAULT_KNODES)
        ]
        
        proc = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
        runtime = time.time() - t0
        
        if proc.returncode != 0:
            print(f"Error running solver for phi={phi}: {proc.stderr}")
            continue
            
        # Read numerical output
        num_gr_file = "output/PY_GdeR.dat"
        num_sk_file = "output/PY_SdeK.dat"
        
        if not os.path.exists(num_gr_file):
            print(f"Error: Output file {num_gr_file} not generated.")
            continue
            
        num_data = np.loadtxt(num_gr_file)
        r_num = num_data[:, 0]
        g_num = num_data[:, 1]
        
        sk_num_data = np.loadtxt(num_sk_file) if os.path.exists(num_sk_file) else None
        
        # Save solution copy in benchmark data dir
        dest_gr_file = os.path.join(OUT_DATA_DIR, f"solution_PY_phi_{phi_str}.dat")
        dest_sk_file = os.path.join(OUT_DATA_DIR, f"sk_PY_phi_{phi_str}.dat")
        dest_sk_ana_file = os.path.join(OUT_DATA_DIR, f"sk_ana_PY_phi_{phi_str}.dat")
        
        np.savetxt(dest_gr_file, num_data, fmt="%.10f", header="r/sigma g(r)")
        if sk_num_data is not None:
            np.savetxt(dest_sk_file, sk_num_data, fmt="%.10f", header="k*sigma S(k)")
            # Generate analytical S(k) on identical k-grid
            k_grid = sk_num_data[:, 0]
            sk_ana_grid = py_sk_analytical(k_grid, phi)
            np.savetxt(dest_sk_ana_file, np.column_stack((k_grid, sk_ana_grid)), fmt="%.10f", header="k*sigma S_ana(k)")
            
        # Interpolate numerical solution onto analytical reference grid (fluid phase r >= 1.0)
        mask_fluid = r_ref >= 1.0
        r_eval = r_ref[mask_fluid]
        g_ref_eval = g_ref[mask_fluid]
        
        g_num_eval = np.interp(r_eval, r_num, g_num)
        
        # g(r) error metrics
        diff_gr = g_num_eval - g_ref_eval
        rmse_gr = np.sqrt(np.mean(diff_gr ** 2))
        max_err_gr = np.max(np.abs(diff_gr))
        
        # S(k) error metrics (for k in [0, 25.0])
        if sk_num_data is not None:
            k_mask = (sk_num_data[:, 0] >= 0.0) & (sk_num_data[:, 0] <= 25.0)
            k_eval = sk_num_data[k_mask, 0]
            sk_num_eval = sk_num_data[k_mask, 1]
            sk_ana_eval = py_sk_analytical(k_eval, phi)
            diff_sk = sk_num_eval - sk_ana_eval
            rmse_sk = np.sqrt(np.mean(diff_sk ** 2))
            max_err_sk = np.max(np.abs(diff_sk))
        else:
            rmse_sk = 0.0
            max_err_sk = 0.0
            
        # Contact value comparison
        g_contact_ana = get_py_contact_analytical(phi)
        g_contact_num = extrapolate_contact_value(r_num, g_num)
        contact_rel_err = abs(g_contact_num - g_contact_ana) / g_contact_ana * 100.0
        
        # Structure factor S(0) comparison
        s0_ana = get_py_sk0_analytical(phi)
        s0_num = sk_num_data[0, 1] if sk_num_data is not None else 0.0
        
        summary_results.append({
            "phi": phi,
            "rmse_gr": rmse_gr,
            "max_err_gr": max_err_gr,
            "g_contact_num": g_contact_num,
            "g_contact_ana": g_contact_ana,
            "contact_rel_err": contact_rel_err,
            "rmse_sk": rmse_sk,
            "max_err_sk": max_err_sk,
            "s0_num": s0_num,
            "s0_ana": s0_ana,
            "runtime": runtime
        })
        
        print(f"{phi:6.2f} | {rmse_gr:12.4e} | {g_contact_num:10.4f} | {g_contact_ana:10.4f} | {contact_rel_err:9.2f}% | {rmse_sk:12.4e} | {max_err_sk:12.4e} | {runtime:8.2f}s")
        
    print("-" * 105)
    
    # Save benchmark summary
    summary_file = os.path.join(OUT_DATA_DIR, "py_benchmark_summary.dat")
    with open(summary_file, "w") as f:
        f.write("# phi  RMSE_gr  MaxErr_gr  g_contact_num  g_contact_ana  ContactErrPct  RMSE_sk  MaxErr_sk  S0_num  S0_ana  Runtime_s\n")
        for res in summary_results:
            f.write(f"{res['phi']:.4f}  {res['rmse_gr']:.6e}  {res['max_err_gr']:.6e}  "
                    f"{res['g_contact_num']:.6f}  {res['g_contact_ana']:.6f}  {res['contact_rel_err']:.4f}  "
                    f"{res['rmse_sk']:.6e}  {res['max_err_sk']:.6e}  "
                    f"{res['s0_num']:.6f}  {res['s0_ana']:.6f}  {res['runtime']:.4f}\n")
                    
    # Generate LaTeX Table
    tex_table_file = "reports/py_hard_sphere_benchmark/table_py_summary.tex"
    with open(tex_table_file, "w") as f:
        f.write("\\begin{table}[htbp]\n")
        f.write("\\centering\n")
        f.write("\\small\n")
        f.write("\\caption{Numerical error metrics for $g(r)$ and $S(k)$, contact values, and execution times for hard spheres under Percus-Yevick closure across seven packing fractions $\\phi$.}\n")
        f.write("\\label{tab:py_summary}\n")
        f.write("\\resizebox{\\textwidth}{!}{\n")
        f.write("\\begin{tabular}{c c c c c c c c}\n")
        f.write("\\hline\\hline\n")
        f.write("$\\phi$ & RMSE $g(r)$ & $g_{\\text{num}}(\\sigma^+)$ & $g_{\\text{ana}}(\\sigma^+)$ & Error $g(\\sigma^+)$ (\\%) & RMSE $S(k)$ & Max $\\Delta S(k)$ & Time (s) \\\\\n")
        f.write("\\hline\n")
        for res in summary_results:
            f.write(f"{res['phi']:.2f} & {res['rmse_gr']:.3e} & "
                    f"{res['g_contact_num']:.4f} & {res['g_contact_ana']:.4f} & {res['contact_rel_err']:.2f}\\% & "
                    f"{res['rmse_sk']:.3e} & {res['max_err_sk']:.3e} & {res['runtime']:.2f} \\\\\n")
        f.write("\\hline\\hline\n")
        f.write("\\end{tabular}\n")
        f.write("}\n")
        f.write("\\end{table}\n")
        
    print(f"Summary written to {summary_file}")
    print(f"LaTeX table written to {tex_table_file}")

if __name__ == "__main__":
    run_benchmark()
