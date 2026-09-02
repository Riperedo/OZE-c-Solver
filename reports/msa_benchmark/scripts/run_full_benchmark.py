import os
import subprocess
import numpy as np

# Base paths
ROOT_DIR = "/home/jinzo/Desktop/codigos/OZE_c_solver"
DATA_ANA_DIR = os.path.join(ROOT_DIR, "test/data_MSA_analitica_HD_dipolar")
BENCH_DIR = os.path.join(ROOT_DIR, "reports/msa_benchmark")
DATA_OUT_DIR = os.path.join(BENCH_DIR, "data")
PLOTS_DIR = os.path.join(BENCH_DIR, "plots")
SCRIPTS_DIR = os.path.join(BENCH_DIR, "scripts")

os.makedirs(DATA_OUT_DIR, exist_ok=True)
os.makedirs(PLOTS_DIR, exist_ok=True)

temps = [10.0, 1.0, 0.1, 0.01]
phis = [0.1, 0.2, 0.3, 0.4, 0.5]

benchmark_records = []

print("=" * 80)
print("STARTING FULL MSA NUMERICAL VS ANALYTICAL BENCHMARK")
print("=" * 80)

for T in temps:
    for phi in phis:
        tag = f"phi_{phi:.1f}_T_{T:.2f}"
        ana_file = os.path.join(DATA_ANA_DIR, f"phi_{phi}_T_{T}.dat")
        if not os.path.exists(ana_file):
            print(f"[WARN] Missing analytical file: {ana_file}")
            continue

        # Run numerical solver
        cmd = [
            os.path.join(ROOT_DIR, "build/facdes_solver"),
            "--closure", "MSA",
            "--potential", "14",
            "--volfactor", str(phi),
            "--temp", str(T),
            "--dipole", "1.0",
            "--nodes", "1024",
            "--knodes", "1024",
            "--rmax", "30.0"
        ]
        
        res = subprocess.run(cmd, cwd=ROOT_DIR, capture_output=True, text=True)
        if res.returncode != 0:
            print(f"[ERROR] Solver failed for {tag}: {res.stderr}")
            continue

        # Load datasets
        ana = np.loadtxt(ana_file)
        num = np.loadtxt(os.path.join(ROOT_DIR, "output/output_dipolar_sk.dat"))

        # Save to benchmark data dir
        np.savetxt(os.path.join(DATA_OUT_DIR, f"num_{tag}.dat"), num,
                   header="k S000 S110 S112 S10 S11", fmt="%.8e")
        np.savetxt(os.path.join(DATA_OUT_DIR, f"ana_{tag}.dat"), ana,
                   header="k S000 S110 S112 S10 S11", fmt="%.8e")

        k_ana = ana[:, 0]
        k_num = num[:, 0]

        # Interpolation evaluation range k in [0.5, 20.0]
        mask = (k_ana >= 0.5) & (k_ana <= 20.0)
        k_eval = k_ana[mask]

        rec = {
            "T": T,
            "phi": phi,
            "tag": tag
        }

        col_names = [(1, "S000"), (2, "S110"), (3, "S112"), (4, "S10"), (5, "S11")]
        for col_idx, name in col_names:
            ana_val = ana[mask, col_idx]
            num_val = np.interp(k_eval, k_num, num[:, col_idx])

            rmse = np.sqrt(np.mean((ana_val - num_val)**2))
            max_err = np.max(np.abs(ana_val - num_val))
            
            # MARE with safe denominator
            denom = np.maximum(np.abs(ana_val), 1e-3)
            mare = np.mean(np.abs(ana_val - num_val) / denom) * 100.0

            rec[f"{name}_rmse"] = rmse
            rec[f"{name}_max"] = max_err
            rec[f"{name}_mare"] = mare

        benchmark_records.append(rec)
        print(f"Done: T={T:5.2f}, phi={phi:.1f} | RMSE S000={rec['S000_rmse']:.3e}, S110={rec['S110_rmse']:.3e}, S112={rec['S112_rmse']:.3e}, S10={rec['S10_rmse']:.3e}, S11={rec['S11_rmse']:.3e}")

print("=" * 80)
print("BENCHMARK SIMULATION COMPLETE")
print("=" * 80)

# Generate LaTeX table code
table_tex_path = os.path.join(BENCH_DIR, "table_benchmark_summary.tex")
with open(table_tex_path, "w") as f:
    f.write("% Auto-generated Benchmark Summary Table\n")
    f.write("\\begin{table}[htbp]\n")
    f.write("\\centering\n")
    f.write("\\small\n")
    f.write("\\caption{Summary of numerical vs. analytical MSA structure factor error metrics ($k \\in [0.5, 20.0]$).}\n")
    f.write("\\label{tab:msa_benchmark_errors}\n")
    f.write("\\begin{tabular}{cc|cc|cc|cc|cc}\n")
    f.write("\\hline\\hline\n")
    f.write("$T^*$ & $\\phi$ & $\\text{RMSE}_{S^{000}}$ & $L^\\infty_{S^{000}}$ & $\\text{RMSE}_{S^{110}}$ & $L^\\infty_{S^{110}}$ & $\\text{RMSE}_{S^{112}}$ & $L^\\infty_{S^{112}}$ & $\\text{RMSE}_{S^0}$ & $\\text{RMSE}_{S^1}$ \\\\\n")
    f.write("\\hline\n")
    current_T = None
    for r in benchmark_records:
        if current_T is not None and r['T'] != current_T:
            f.write("\\hline\n")
        current_T = r['T']
        f.write(f"{r['T']:.2f} & {r['phi']:.1f} & "
                f"{r['S000_rmse']:.2e} & {r['S000_max']:.2e} & "
                f"{r['S110_rmse']:.2e} & {r['S110_max']:.2e} & "
                f"{r['S112_rmse']:.2e} & {r['S112_max']:.2e} & "
                f"{r['S10_rmse']:.2e} & {r['S11_rmse']:.2e} \\\\\n")
    f.write("\\hline\\hline\n")
    f.write("\\end{tabular}\n")
    f.write("\\end{table}\n")

print(f"Generated LaTeX table at: {table_tex_path}")
