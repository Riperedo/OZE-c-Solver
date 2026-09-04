#!/usr/bin/env python3
"""
Phase 3 Benchmark: Direct Analytical Structure Factor Initialization (Warm-Start)
Evaluates dipolar hard spheres with direct analytical S(k) ingestion (no continuation ramps).
"""

import os
import subprocess
import numpy as np

PHIS = [0.1, 0.2, 0.3, 0.4, 0.5]
TEMPS = [0.1, 0.01]
OUTPUT_DIR = "reports/msa_benchmark/data/phase3"
ANA_DIR = "data/input_sk"

os.makedirs(OUTPUT_DIR, exist_ok=True)

def compute_rmse(num_file, ana_file):
    if not os.path.exists(num_file) or not os.path.exists(ana_file):
        return None
    try:
        num = np.loadtxt(num_file)
        ana = np.loadtxt(ana_file)
        
        k_num = num[:, 0]
        s000_num = num[:, 1]
        s110_num = num[:, 2]
        s112_num = num[:, 3]
        
        k_ana = ana[:, 0]
        s000_ana = ana[:, 1]
        s110_ana = ana[:, 2]
        s112_ana = ana[:, 3]
        
        s000_num_interp = np.interp(k_ana, k_num, s000_num)
        s110_num_interp = np.interp(k_ana, k_num, s110_num)
        s112_num_interp = np.interp(k_ana, k_num, s112_num)
        
        rmse_000 = np.sqrt(np.mean((s000_num_interp - s000_ana)**2))
        rmse_110 = np.sqrt(np.mean((s110_num_interp - s110_ana)**2))
        rmse_112 = np.sqrt(np.mean((s112_num_interp - s112_ana)**2))
        
        s1_num = s110_num_interp - s112_num_interp
        s1_ana = s110_ana - s112_ana
        rmse_1 = np.sqrt(np.mean((s1_num - s1_ana)**2))
        
        return rmse_000, rmse_110, rmse_112, rmse_1
    except Exception as e:
        return None

results = []
print("=" * 80)
print("=== Phase 3: Direct Analytical S(k) Warm-Start Benchmark ===")
print("=" * 80)
print(f"{'phi':<5} {'T*':<6} {'Coupling':<10} {'Iter':<6} {'L2 Error':<14} {'RMSE(S110)':<12} {'RMSE(S112)':<12} {'Status'}")
print("-" * 80)

for T in TEMPS:
    beta_mu2 = 1.0 / T
    for phi in PHIS:
        ana_cand = f"{ANA_DIR}/sk_phi_{phi:.1f}_T_{T:g}.dat"
        if not os.path.exists(ana_cand):
            ana_cand = f"{ANA_DIR}/sk_phi_{phi:.1f}_T_{T:.2f}.dat"
            
        cmd = [
            "./build/facdes_solver",
            "--closure", "MSA",
            "--potential", "14",
            "--volfactor", f"{phi:.1f}",
            "--temp", f"{T:g}",
            "--dipole", "1.0",
            "--nodes", "1024",
            "--knodes", "1024",
            "--rmax", "30.0",
            "--init-sk", ana_cand
        ]
        proc = subprocess.run(cmd, capture_output=True, text=True)
        
        lines = [l for l in proc.stdout.split("\n") if "Iter" in l]
        last_iter_str = lines[-1] if lines else "N/A"
        
        iter_num = 0
        l2_err = 1e9
        status = "Failed"
        if "Iter" in last_iter_str:
            parts = last_iter_str.replace(":", "").split()
            try:
                iter_num = int(parts[1])
                l2_err = float(parts[3])
                if l2_err < 1e-5 and iter_num < 590:
                    status = "Converged"
                else:
                    status = "Stalled"
            except:
                pass
        
        out_sk = f"{OUTPUT_DIR}/sk_phi_{phi:.1f}_T_{T:g}.dat"
        if os.path.exists("output/output_dipolar_sk.dat"):
            subprocess.run(["cp", "output/output_dipolar_sk.dat", out_sk])
            
        rmses = compute_rmse(out_sk, ana_cand)
        rmse_110 = rmses[1] if rmses else 999.0
        rmse_112 = rmses[2] if rmses else 999.0
        
        results.append((phi, T, beta_mu2, iter_num, l2_err, rmse_110, rmse_112, status))
        print(f"{phi:<5.1f} {T:<6.2f} {beta_mu2:<10.2f} {iter_num:<6d} {l2_err:<14.2e} {rmse_110:<12.4e} {rmse_112:<12.4e} {status}")

print("=" * 80)
# Save summary table
with open(f"{OUTPUT_DIR}/summary_phase3.dat", "w") as f:
    f.write("# phi T beta_mu2 iters l2_err rmse_110 rmse_112 status\n")
    for r in results:
        f.write(f"{r[0]:.2f} {r[1]:.4f} {r[2]:.4f} {r[3]} {r[4]:.6e} {r[5]:.6e} {r[6]:.6e} {r[7]}\n")
print(f"Summary written to {OUTPUT_DIR}/summary_phase3.dat")
