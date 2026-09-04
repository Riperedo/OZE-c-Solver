#!/usr/bin/env python3
"""
Monte Carlo Benchmark Runner for Dipolar Hard Sphere Solver.
Evaluates MSA, LHNC, and QHNC closures against Monte Carlo simulation data
from Fries & Patey (J. Chem. Phys. 82, 429, 1985).

Phases:
- Phase 1: Cold Start (c^(0) = 0, no ramps)
- Phase 2: Geometric Temperature Continuation (Annealing ramps)
- Phase 3: Direct Analytical / Warm-Start Initialization
"""

import os
import sys
import time
import subprocess
import numpy as np
from scipy.interpolate import interp1d

# Paths
BASE_DIR = os.path.abspath(os.path.join(os.path.dirname(__file__), "../../.."))
MC_DIR = os.path.join(BASE_DIR, "test/data_Monte_Carlo_Simulations/FriesPatey_JCP_824291985")
OUTPUT_DIR = os.path.join(BASE_DIR, "reports/monte_carlo_benchmark/data")
os.makedirs(OUTPUT_DIR, exist_ok=True)

EXECUTABLE = os.path.join(BASE_DIR, "build/facdes_solver")

# Thermodynamic States (Fries & Patey 1985)
# rho* = rho * sigma^3 = 0.8 -> packing fraction phi = (pi/6) * 0.8
PHI = 0.8 * np.pi / 6.0

STATES = {
    "state_2.0": {
        "name": r"$\rho^*=0.8, \mu^{*2}=2.0$",
        "mu2": 2.0,
        "temp": 0.5,
        "phi": PHI,
        "tag": "rho_0.8_mu2_2.0"
    },
    "state_2.75": {
        "name": r"$\rho^*=0.8, \mu^{*2}=2.75$",
        "mu2": 2.75,
        "temp": 1.0 / 2.75,
        "phi": PHI,
        "tag": "rho_0.8_mu2_2.75"
    }
}

CLOSURES = ["MSA", "LHNC", "QHNC", "RHNC"]

def run_solver(closure, temp, phi, use_ramp=False, temp_start=10.0, temp_steps=15, init_sk=None):
    """Executes the C solver with specified arguments."""
    cmd = [
        EXECUTABLE,
        "--closure", closure,
        "--potential", "14",
        "--volfactor", str(phi),
        "--temp", str(temp),
        "--dipole", "1.0",
        "--nodes", "1024",
        "--knodes", "1024",
        "--rmax", "30.0"
    ]
    if use_ramp:
        cmd.extend(["--ramp", "--temp-start", str(temp_start), "--temp-steps", str(temp_steps)])
    if init_sk:
        cmd.extend(["--init-sk", init_sk])

    t0 = time.time()
    res = subprocess.run(cmd, capture_output=True, text=True, cwd=BASE_DIR)
    t1 = time.time()

    elapsed = t1 - t0
    stdout = res.stdout

    # Parse iterations and final error
    iter_lines = [l for l in stdout.split("\n") if "Iter" in l]
    total_iters = len(iter_lines)
    final_error = 1e9
    status = "Stalled"

    if iter_lines:
        last_line = iter_lines[-1]
        try:
            parts = last_line.split(":")
            it_num = int(parts[0].replace("Iter", "").strip())
            err_val = float(parts[1].replace("Error =", "").strip())
            final_error = err_val
            total_iters = it_num
            if final_error < 1e-5:
                status = "Converged"
        except Exception:
            pass

    return {
        "status": status,
        "iterations": total_iters,
        "final_error": final_error,
        "elapsed_sec": elapsed,
        "stdout": stdout
    }

def compute_mc_metrics(state_key, closure_tag, r_th, g000_th, h110_th, h112_th):
    """Interpolates numerical results onto Monte Carlo points and calculates RMSE and L_inf."""
    metrics = {}
    
    interp_g000 = interp1d(r_th, g000_th, kind="cubic", bounds_error=False, fill_value="extrapolate")
    interp_h110 = interp1d(r_th, h110_th, kind="cubic", bounds_error=False, fill_value="extrapolate")
    interp_h112 = interp1d(r_th, h112_th, kind="cubic", bounds_error=False, fill_value="extrapolate")

    if state_key == "state_2.75":
        # Fig 1a: g000 contact
        f1a = np.loadtxt(os.path.join(MC_DIR, "fig1a_rho_0.8_mu2_2.75.dat"))
        diff_1a = interp_g000(f1a[:,0]) - f1a[:,1]
        metrics["rmse_g000_contact"] = np.sqrt(np.mean(diff_1a**2))
        metrics["linf_g000_contact"] = np.max(np.abs(diff_1a))

        # Fig 1b: g000 medium
        f1b = np.loadtxt(os.path.join(MC_DIR, "fig1b_rho_0.8_mu2_2.75.dat"))
        diff_1b = interp_g000(f1b[:,0]) - f1b[:,1]
        metrics["rmse_g000_medium"] = np.sqrt(np.mean(diff_1b**2))
        metrics["linf_g000_medium"] = np.max(np.abs(diff_1b))

    elif state_key == "state_2.0":
        # Fig 2a: g000 contact
        f2a = np.loadtxt(os.path.join(MC_DIR, "fig2a_rho_0.8_mu2_2.0.dat"))
        diff_2a = interp_g000(f2a[:,0]) - f2a[:,1]
        metrics["rmse_g000_contact"] = np.sqrt(np.mean(diff_2a**2))
        metrics["linf_g000_contact"] = np.max(np.abs(diff_2a))

        # Fig 2b: g000 medium
        f2b = np.loadtxt(os.path.join(MC_DIR, "fig2b_rho_0.8_mu2_2.0.dat"))
        diff_2b = interp_g000(f2b[:,0]) - f2b[:,1]
        metrics["rmse_g000_medium"] = np.sqrt(np.mean(diff_2b**2))
        metrics["linf_g000_medium"] = np.max(np.abs(diff_2b))

        # Fig 3: h110
        f3 = np.loadtxt(os.path.join(MC_DIR, "fig3_rho_0.8_mu2_2.0.dat"))
        diff_3 = interp_h110(f3[:,0]) - f3[:,1]
        metrics["rmse_h110"] = np.sqrt(np.mean(diff_3**2))
        metrics["linf_h110"] = np.max(np.abs(diff_3))

        # Fig 4: h112
        f4 = np.loadtxt(os.path.join(MC_DIR, "fig4_rho_0.8_mu2_2.0.dat"))
        diff_4 = interp_h112(f4[:,0]) - f4[:,1]
        metrics["rmse_h112"] = np.sqrt(np.mean(diff_4**2))
        metrics["linf_h112"] = np.max(np.abs(diff_4))

    return metrics

def main():
    print("================================================================================")
    print("  MONTE CARLO BENCHMARK SUITE: MSA vs LHNC vs QHNC (Fries & Patey 1985)")
    print("================================================================================")

    all_phase_results = {"phase1": {}, "phase2": {}, "phase3": {}}

    # -------------------------------------------------------------------------
    # PHASE 1: Cold Start (c^(0) = 0)
    # -------------------------------------------------------------------------
    print("\n>>> Executing Phase 1: Standard Cold-Start Solver (c^(0) = 0)...")
    for s_key, s_data in STATES.items():
        for closure in CLOSURES:
            print(f"  Running [{closure:<4}] for {s_data['name']} (Cold Start)...", end="", flush=True)
            res = run_solver(closure, s_data["temp"], s_data["phi"], use_ramp=False)
            all_phase_results["phase1"][(s_key, closure)] = res
            print(f" Status: {res['status']}, Iters: {res['iterations']}, Error: {res['final_error']:.2e}")

    # -------------------------------------------------------------------------
    # PHASE 2: Temperature Continuation (Annealing Ramps)
    # -------------------------------------------------------------------------
    print("\n>>> Executing Phase 2: Geometric Temperature Continuation Ramps...")
    for s_key, s_data in STATES.items():
        for closure in CLOSURES:
            print(f"  Running [{closure:<4}] for {s_data['name']} (Continuation)...", end="", flush=True)
            res = run_solver(closure, s_data["temp"], s_data["phi"], use_ramp=True, temp_start=10.0, temp_steps=15)
            all_phase_results["phase2"][(s_key, closure)] = res
            
            # Save converged profile to data/
            if res["status"] == "Converged":
                out_path = os.path.join(OUTPUT_DIR, f"solution_{s_data['tag']}_{closure}.dat")
                subprocess.run(["cp", "output/output_dipolar.dat", out_path], cwd=BASE_DIR)
                sk_path = os.path.join(OUTPUT_DIR, f"sk_{s_data['tag']}_{closure}.dat")
                subprocess.run(["cp", "output/output_dipolar_sk.dat", sk_path], cwd=BASE_DIR)
            print(f" Status: {res['status']}, Iters: {res['iterations']}, Error: {res['final_error']:.2e}, Time: {res['elapsed_sec']:.2f}s")

    # -------------------------------------------------------------------------
    # PHASE 3: Analytical / Warm-Start Direct Ingestion
    # -------------------------------------------------------------------------
    print("\n>>> Executing Phase 3: Analytical / Warm-Start Direct Ingestion...")
    for s_key, s_data in STATES.items():
        # Use LHNC S(k) as warm-start provider
        sk_init_file = os.path.join(OUTPUT_DIR, f"sk_{s_data['tag']}_LHNC.dat")
        for closure in CLOSURES:
            print(f"  Running [{closure:<4}] for {s_data['name']} (Warm-Start)...", end="", flush=True)
            res = run_solver(closure, s_data["temp"], s_data["phi"], use_ramp=False, init_sk=sk_init_file)
            all_phase_results["phase3"][(s_key, closure)] = res
            print(f" Status: {res['status']}, Iters: {res['iterations']}, Error: {res['final_error']:.2e}")

    # -------------------------------------------------------------------------
    # Error Metrics Computation against Monte Carlo
    # -------------------------------------------------------------------------
    print("\n>>> Computing Error Metrics against Monte Carlo Datasets...")
    summary_lines = []
    summary_lines.append("# State Closure Status Phase1_Iters Phase2_Iters Phase3_Iters RMSE_g000_contact RMSE_g000_med RMSE_h110 RMSE_h112")
    
    mc_table_tex = []
    
    for s_key, s_data in STATES.items():
        for closure in CLOSURES:
            sol_path = os.path.join(OUTPUT_DIR, f"solution_{s_data['tag']}_{closure}.dat")
            if not os.path.exists(sol_path):
                continue
            data = np.loadtxt(sol_path)
            r = data[:,0]
            h000 = data[:,1]
            h110 = data[:,2]
            h112 = data[:,3]
            g000 = h000 + 1.0

            metrics = compute_mc_metrics(s_key, closure, r, g000, h110, h112)
            
            p1_it = all_phase_results["phase1"][(s_key, closure)]["iterations"]
            p2_it = all_phase_results["phase2"][(s_key, closure)]["iterations"]
            p3_it = all_phase_results["phase3"][(s_key, closure)]["iterations"]
            
            rmse_c = metrics.get("rmse_g000_contact", 0.0)
            rmse_m = metrics.get("rmse_g000_medium", 0.0)
            rmse_110 = metrics.get("rmse_h110", 0.0)
            rmse_112 = metrics.get("rmse_h112", 0.0)
            
            summary_lines.append(
                f"{s_data['tag']} {closure} Converged {p1_it} {p2_it} {p3_it} "
                f"{rmse_c:.6e} {rmse_m:.6e} {rmse_110:.6e} {rmse_112:.6e}"
            )
            
            # Form LaTeX table row
            if s_key == "state_2.0":
                row_tex = f"$\\mu^{{*2}}=2.00$ & {closure} & {p1_it} & {p2_it} & {p3_it} & {rmse_c:.4f} & {rmse_m:.4f} & {rmse_110:.4f} & {rmse_112:.4f} \\\\"
            else:
                row_tex = f"$\\mu^{{*2}}=2.75$ & {closure} & {p1_it} & {p2_it} & {p3_it} & {rmse_c:.4f} & {rmse_m:.4f} & --- & --- \\\\"
            mc_table_tex.append(row_tex)

    # Save summary dat
    summary_path = os.path.join(OUTPUT_DIR, "mc_benchmark_summary.dat")
    with open(summary_path, "w") as f:
        f.write("\n".join(summary_lines) + "\n")
    print(f"Saved: {summary_path}")

    # Save LaTeX table fragment
    table_tex_path = os.path.join(BASE_DIR, "reports/monte_carlo_benchmark/table_mc_summary.tex")
    with open(table_tex_path, "w") as f:
        f.write("% Auto-generated Monte Carlo error comparison table\n")
        f.write("\\begin{table}[htbp]\n")
        f.write("    \\centering\n")
        f.write("    \\caption{Comparative convergence and error performance of MSA, LHNC, QHNC, and RHNC against Monte Carlo simulations (Fries \\& Patey 1985) across all three evaluation phases at $\\rho^* = 0.8$.}\n")
        f.write("    \\label{tab:mc_summary}\n")
        f.write("    \\vspace{2mm}\n")
        f.write("    \\begin{tabular}{ccccccccc}\n")
        f.write("        \\hline\\hline\n")
        f.write("        State & Closure & Phase 1 & Phase 2 & Phase 3 & $\\text{RMSE}(g_{\\text{cont}}^{000})$ & $\\text{RMSE}(g_{\\text{med}}^{000})$ & $\\text{RMSE}(h^{110})$ & $\\text{RMSE}(h^{112})$ \\\\\n")
        f.write("        \\hline\n")
        n_cl = len(CLOSURES)
        f.write("        \\multicolumn{9}{l}{\\textbf{State 1: } $\\rho^* = 0.8, \\mu^{*2} = 2.75, T^* = 0.3636$} \\\\\n")
        for row in mc_table_tex[n_cl:]:
            f.write("        " + row + "\n")
        f.write("        \\hline\n")
        f.write("        \\multicolumn{9}{l}{\\textbf{State 2: } $\\rho^* = 0.8, \\mu^{*2} = 2.00, T^* = 0.5000$} \\\\\n")
        for row in mc_table_tex[:n_cl]:
            f.write("        " + row + "\n")
        f.write("        \\hline\\hline\n")
        f.write("    \\end{tabular}\n")
        f.write("\\end{table}\n")
    print(f"Saved: {table_tex_path}")

if __name__ == "__main__":
    main()
