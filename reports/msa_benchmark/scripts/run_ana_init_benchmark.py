#!/usr/bin/env python3
"""
Benchmark OZE dipolar solver with direct analytical structure factor initialization (Warm-Start).
Evaluates T* in {0.10, 0.01} across phi in {0.1, 0.2, 0.3, 0.4, 0.5}.
"""

import os, subprocess, shutil
import numpy as np

phis = [0.1, 0.2, 0.3, 0.4, 0.5]
temps = [0.10, 0.01]

results = []

for T in temps:
    for phi in phis:
        # Match input filename (check if T=0.1 or 0.10)
        cand1 = f"data/input_sk/sk_phi_{phi:.1f}_T_{T:.2f}.dat"
        cand2 = f"data/input_sk/sk_phi_{phi:.1f}_T_{T:g}.dat"
        init_file = cand1 if os.path.exists(cand1) else cand2
        
        cmd = [
            "./build/facdes_solver",
            "--closure", "MSA",
            "--potential", "14",
            "--volfactor", str(phi),
            "--temp", str(T),
            "--dipole", "1.0",
            "--nodes", "1024",
            "--knodes", "1024",
            "--rmax", "30.0",
            "--init-sk", init_file
        ]
        
        print(f"\n==========================================")
        print(f"Running: phi={phi:.1f}, T*={T:.2f}")
        print(f"Init File: {init_file}")
        print(f"==========================================")
        
        res = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
        
        # Save output
        out_saved = f"reports/msa_benchmark/data/num_ana_init_phi_{phi:.1f}_T_{T:.2f}.dat"
        if os.path.exists("output/output_dipolar_sk.dat"):
            shutil.copy("output/output_dipolar_sk.dat", out_saved)
            
            # Compute RMSE against analytical data
            ana_data = np.loadtxt(init_file)
            num_data = np.loadtxt(out_saved)
            
            k_ana, s000_a, s110_a, s112_a = ana_data[:,0], ana_data[:,1], ana_data[:,2], ana_data[:,3]
            k_num, s000_n, s110_n, s112_n = num_data[:,0], num_data[:,1], num_data[:,2], num_data[:,3]
            s1_n = num_data[:, 5]
            s1_a = (s110_a + 1.0) - s112_a
            
            mask = (k_ana >= 0.5) & (k_ana <= 20.0)
            k_eval = k_ana[mask]
            
            s000_interp = np.interp(k_eval, k_num, s000_n)
            s110_interp = np.interp(k_eval, k_num, s110_n)
            s112_interp = np.interp(k_eval, k_num, s112_n)
            s1_interp = np.interp(k_eval, k_num, s1_n)
            
            rmse_000 = np.sqrt(np.mean((s000_interp - s000_a[mask])**2))
            rmse_110 = np.sqrt(np.mean((s110_interp - s110_a[mask])**2))
            rmse_112 = np.sqrt(np.mean((s112_interp - s112_a[mask])**2))
            rmse_s1  = np.sqrt(np.mean((s1_interp - s1_a[mask])**2))
            
            # Parse iterations from stdout
            lines = res.stdout.strip().split("\n")
            last_iter_line = [l for l in lines if l.startswith("Iter") or "completed in" in l]
            summary = last_iter_line[-1] if last_iter_line else "Unknown"
            
            results.append({
                "phi": phi, "T": T,
                "rmse_000": rmse_000, "rmse_110": rmse_110,
                "rmse_112": rmse_112, "rmse_s1": rmse_s1,
                "status": summary
            })
            print(f"Result: {summary}")
            print(f"RMSE: S000={rmse_000:.4e}, S110={rmse_110:.4e}, S112={rmse_112:.4e}, S1={rmse_s1:.4e}")
        else:
            print(f"FAILED to produce output.")

# Summary Table
print("\n" + "="*80)
print(f"{'phi':<6} {'T*':<6} {'RMSE S000':<12} {'RMSE S110':<12} {'RMSE S112':<12} {'RMSE S1':<12} {'Status'}")
print("="*80)
for r in results:
    print(f"{r['phi']:<6.1f} {r['T']:<6.2f} {r['rmse_000']:<12.4e} {r['rmse_110']:<12.4e} {r['rmse_112']:<12.4e} {r['rmse_s1']:<12.4e} {r['status']}")
print("="*80)
