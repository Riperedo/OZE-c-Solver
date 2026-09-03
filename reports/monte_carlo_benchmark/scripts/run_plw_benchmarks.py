#!/usr/bin/env python3
"""
Automated Benchmark Suite for Patey, Levesque & Weis (1979) Monte Carlo Data.
Evaluates MSA, LHNC, and QHNC across Three Evaluation Phases:
  - Phase 1: Cold Start (c(0) = 0)
  - Phase 2: Temperature Continuation (Annealing Ramps)
  - Phase 3: Warm-Start Ingestion

Computes quantitative RMSE and L_inf errors against all 10 digitized MC figures:
  - Fig 1: g000 (rho=0.15, mu2=2.0)
  - Fig 2: g000 (rho=0.40, mu2=2.75)
  - Fig 3: g000 (rho=0.60, mu2=2.75)
  - Fig 4: g000 (rho=0.80, mu2=2.75)
  - Fig 6: h110 (rho=0.15, mu2=2.0)
  - Fig 7: h110 (rho=0.40, mu2=2.75)
  - Fig 8: h110 (rho=0.60, mu2=2.75)
  - Fig 9: h112 (rho=0.15, mu2=2.0)
  - Fig 10: h112 (rho=0.40, mu2=2.75)
  - Fig 11: h112 (rho=0.60, mu2=2.75)
"""

import os
import sys
import math
import subprocess
import shutil
import numpy as np

WORKSPACE_DIR = "/home/jinzo/Desktop/codigos/OZE_c_solver"
BUILD_DIR = os.path.join(WORKSPACE_DIR, "build")
EXECUTABLE = os.path.join(BUILD_DIR, "facdes_solver")
OUTPUT_DIR = os.path.join(WORKSPACE_DIR, "output")
REPORT_DATA_DIR = os.path.join(WORKSPACE_DIR, "reports", "monte_carlo_benchmark", "data")
PLW_MC_DIR = os.path.join(WORKSPACE_DIR, "test", "data_Monte_Carlo_Simulations", "PateyLevesqueWeis_1979")
TABLE_TEX_PATH = os.path.join(WORKSPACE_DIR, "reports", "monte_carlo_benchmark", "table_plw_summary.tex")

os.makedirs(REPORT_DATA_DIR, exist_ok=True)

# 4 Thermodynamic State Points in Patey, Levesque & Weis (1979)
PLW_STATES = {
    "rho_0.15_mu2_2.0": {
        "rho": 0.15, "mu2": 2.0, "temp": 1.0, "dipole": math.sqrt(2.0),
        "phi": math.pi * 0.15 / 6.0,
        "figures": {"g000": "Fig1_rho_0.15_mu2_2.0.dat", "h110": "Fig6_rho_0.15_mu2_2.0.dat", "h112": "Fig9_rho_0.15_mu2_2.0.dat"}
    },
    "rho_0.4_mu2_2.75": {
        "rho": 0.40, "mu2": 2.75, "temp": 1.0, "dipole": math.sqrt(2.75),
        "phi": math.pi * 0.40 / 6.0,
        "figures": {"g000": "Fig2_rho_0.4_mu2_2.75.dat", "h110": "Fig7_rho_0.4_mu2_2.75.dat", "h112": "Fig10_rho_0.4_mu_2.75.dat"}
    },
    "rho_0.6_mu2_2.75": {
        "rho": 0.60, "mu2": 2.75, "temp": 1.0, "dipole": math.sqrt(2.75),
        "phi": math.pi * 0.60 / 6.0,
        "figures": {"g000": "Fig3_rho_0.6_mu2_2.75.dat", "h110": "Fig8_rho_0.6_mu2_2.75.dat", "h112": "Fig11_rho_0.6_mu2_2.75.dat"}
    },
    "rho_0.8_mu2_2.75": {
        "rho": 0.80, "mu2": 2.75, "temp": 1.0, "dipole": math.sqrt(2.75),
        "phi": math.pi * 0.80 / 6.0,
        "figures": {"g000": "Fig4_rho_0.8_mu2_2.75.dat"}
    }
}

CLOSURES = [
    {"name": "MSA", "id": 0},
    {"name": "LHNC", "id": 1},
    {"name": "QHNC", "id": 2}
]

def run_solver(closure_name, phi, temp, dipole, rmax=30.0, nodes=4096, knodes=1024, ramp=False, temp_start=10.0, temp_steps=12, init_sk=None):
    cmd = [
        EXECUTABLE,
        "--potential", "14",
        "--closure", str(closure_name),
        "--volfactor", f"{phi:.8f}",
        "--temp", f"{temp:.8f}",
        "--dipole", f"{dipole:.8f}",
        "--nodes", str(nodes),
        "--knodes", str(knodes),
        "--rmax", f"{rmax:.1f}"
    ]
    if ramp:
        cmd.extend(["--ramp", "--temp-start", f"{temp_start:.1f}", "--temp-steps", str(temp_steps)])
    if init_sk and os.path.exists(init_sk):
        cmd.extend(["--init-sk", init_sk])

    res = subprocess.run(cmd, capture_output=True, text=True)
    
    # Parse iterations and error
    iters = 0
    final_err = 1.0
    converged = False
    for line in res.stdout.split("\n"):
        if "Iter" in line and "Error =" in line:
            parts = line.split()
            try:
                iters = int(parts[1].replace(":", ""))
                final_err = float(parts[-1])
            except (ValueError, IndexError):
                pass
        if "Written output/output_dipolar.dat" in line:
            converged = True

    if not converged and iters >= 200:
        return {"converged": False, "iters": "Stalled", "error": final_err, "raw": res.stdout}

    return {"converged": converged, "iters": iters, "error": final_err, "raw": res.stdout}

def compute_rmse(num_r, num_val, mc_file):
    if not os.path.exists(mc_file):
        return None
    mc_data = np.loadtxt(mc_file)
    r_mc = mc_data[:, 0]
    y_mc = mc_data[:, 1]
    
    # Interpolate numerical solution onto MC grid
    y_num_interp = np.interp(r_mc, num_r, num_val)
    
    rmse = np.sqrt(np.mean((y_num_interp - y_mc)**2))
    linf = np.max(np.abs(y_num_interp - y_mc))
    return {"rmse": rmse, "linf": linf}

def main():
    print("=" * 80)
    print("  RUNNING MULTI-CLOSURE BENCHMARK: PATEY, LEVESQUE & WEIS (1979)")
    print("=" * 80)

    results = []

    for state_key, state in PLW_STATES.items():
        print(f"\n--- State: {state_key} (rho={state['rho']}, mu2={state['mu2']}, T=1.0) ---")

        for cl in CLOSURES:
            cl_name = cl["name"]
            cl_id = cl["id"]
            print(f"  Evaluating Closure: {cl_name} (ID: {cl_id})")

            # Phase 1: Cold Start
            res_p1 = run_solver(cl_name, state["phi"], state["temp"], state["dipole"])
            p1_iters = res_p1["iters"] if res_p1["converged"] else "Stalled"
            print(f"    Phase 1 (Cold Start):    {p1_iters} iters")

            # Phase 2: Temperature Continuation
            res_p2 = run_solver(cl_name, state["phi"], state["temp"], state["dipole"], ramp=True, temp_start=10.0, temp_steps=12)
            p2_iters = res_p2["iters"] if res_p2["converged"] else "Stalled"
            print(f"    Phase 2 (Continuation):  {p2_iters} iters")

            # Save converged solution from Phase 2
            sol_out = os.path.join(OUTPUT_DIR, "output_dipolar.dat")
            sk_out = os.path.join(OUTPUT_DIR, "output_dipolar_sk.dat")

            target_sol = os.path.join(REPORT_DATA_DIR, f"solution_plw_{state_key}_{cl_name}.dat")
            target_sk = os.path.join(REPORT_DATA_DIR, f"sk_plw_{state_key}_{cl_name}.dat")
            if os.path.exists(sol_out):
                shutil.copyfile(sol_out, target_sol)
            if os.path.exists(sk_out):
                shutil.copyfile(sk_out, target_sk)

            # Phase 3: Warm-Start Ingestion from prior S(k)
            res_p3 = run_solver(cl_name, state["phi"], state["temp"], state["dipole"], init_sk=target_sk)
            p3_iters = res_p3["iters"] if res_p3["converged"] else "Stalled"
            print(f"    Phase 3 (Warm-Start):    {p3_iters} iters")

            # Load solution to compute errors against MC figures
            errors = {}
            if os.path.exists(target_sol):
                sol_data = np.loadtxt(target_sol)
                r_arr = sol_data[:, 0]
                g000_arr = sol_data[:, 1] + 1.0
                h110_arr = sol_data[:, 2]
                h112_arr = sol_data[:, 3]

                # Compare with Fig g000
                if "g000" in state["figures"]:
                    mc_f = os.path.join(PLW_MC_DIR, state["figures"]["g000"])
                    e = compute_rmse(r_arr, g000_arr, mc_f)
                    errors["g000"] = e["rmse"] if e else None

                # Compare with Fig h110
                if "h110" in state["figures"]:
                    mc_f = os.path.join(PLW_MC_DIR, state["figures"]["h110"])
                    e = compute_rmse(r_arr, h110_arr, mc_f)
                    errors["h110"] = e["rmse"] if e else None

                # Compare with Fig h112
                if "h112" in state["figures"]:
                    mc_f = os.path.join(PLW_MC_DIR, state["figures"]["h112"])
                    e = compute_rmse(r_arr, h112_arr, mc_f)
                    errors["h112"] = e["rmse"] if e else None

            print(f"    Errors (RMSE): g000={errors.get('g000')}, h110={errors.get('h110')}, h112={errors.get('h112')}")

            results.append({
                "state": state_key,
                "rho": state["rho"],
                "mu2": state["mu2"],
                "closure": cl_name,
                "p1_iters": p1_iters,
                "p2_iters": p2_iters,
                "p3_iters": p3_iters,
                "rmse_g000": errors.get("g000"),
                "rmse_h110": errors.get("h110"),
                "rmse_h112": errors.get("h112")
            })

    # Save summary data
    summary_txt = os.path.join(REPORT_DATA_DIR, "plw_benchmark_summary.dat")
    with open(summary_txt, "w") as f:
        f.write("# State Closure Phase1 Phase2 Phase3 RMSE_g000 RMSE_h110 RMSE_h112\n")
        for r in results:
            g_str = f"{r['rmse_g000']:.5f}" if r['rmse_g000'] is not None else "---"
            h110_str = f"{r['rmse_h110']:.5f}" if r['rmse_h110'] is not None else "---"
            h112_str = f"{r['rmse_h112']:.5f}" if r['rmse_h112'] is not None else "---"
            f.write(f"{r['state']} {r['closure']} {r['p1_iters']} {r['p2_iters']} {r['p3_iters']} {g_str} {h110_str} {h112_str}\n")
    print(f"\n✓ Saved summary dataset to {summary_txt}")

    # Generate LaTeX Table Fragment
    with open(TABLE_TEX_PATH, "w") as f:
        f.write("% Auto-generated Patey, Levesque & Weis (1979) error comparison table\n")
        f.write("\\begin{table}[htbp]\n")
        f.write("    \\centering\n")
        f.write("    \\caption{Comparative convergence and error performance of MSA, LHNC, and QHNC against Monte Carlo simulations of Patey, Levesque \\& Weis (1979) across all density regimes at $T^* = 1.0$.}\n")
        f.write("    \\label{tab:plw_summary}\n")
        f.write("    \\vspace{2mm}\n")
        f.write("    \\resizebox{\\textwidth}{!}{\n")
        f.write("    \\begin{tabular}{cccccccc}\n")
        f.write("        \\toprule\n")
        f.write("        \\textbf{State Point} & \\textbf{Closure} & \\textbf{Phase 1 (Cold)} & \\textbf{Phase 2 (Cont.)} & \\textbf{Phase 3 (Warm)} & $\\text{RMSE}(g^{000})$ & $\\text{RMSE}(h^{110})$ & $\\text{RMSE}(h^{112})$ \\\\\n")
        f.write("        \\midrule\n")

        current_st = None
        for r in results:
            if r["state"] != current_st:
                current_st = r["state"]
                f.write(f"        \\multicolumn{{8}}{{l}}{{\\textit{{{current_st.replace('_', ' ')}}}: $\\rho^* = {r['rho']}, \\mu^{{*2}} = {r['mu2']}$}} \\\\\n")

            g_str = f"{r['rmse_g000']:.4f}" if r['rmse_g000'] is not None else "---"
            h110_str = f"{r['rmse_h110']:.4f}" if r['rmse_h110'] is not None else "---"
            h112_str = f"{r['rmse_h112']:.4f}" if r['rmse_h112'] is not None else "---"

            f.write(f"        $\\rho^*={r['rho']}$ & {r['closure']} & {r['p1_iters']} & {r['p2_iters']} & {r['p3_iters']} & {g_str} & {h110_str} & {h112_str} \\\\\n")

        f.write("        \\bottomrule\n")
        f.write("    \\end{tabular}\n")
        f.write("    }\n")
        f.write("\\end{table}\n")

    print(f"✓ Generated LaTeX table fragment to {TABLE_TEX_PATH}")

if __name__ == "__main__":
    main()
