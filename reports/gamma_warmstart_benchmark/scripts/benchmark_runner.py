#!/usr/bin/env python3
"""
benchmark_runner.py
Automated benchmark comparing cold-start vs. warm-start performance in OZE-c-Solver.
Evaluates Percus-Yevick closure for Hard Spheres across volume fractions phi in {0.50, 0.55, 0.60, 0.64}.
Pure Python implementation (no external dependencies required).
"""

import os
import sys
import time
import subprocess
import math
import statistics

SOLVER_BIN = os.path.abspath(os.path.join(os.path.dirname(__file__), "../../../build/facdes_solver"))
DATA_DIR = os.path.abspath(os.path.join(os.path.dirname(__file__), "../data"))
OUTPUT_DIR = os.path.abspath(os.path.join(os.path.dirname(__file__), "../../../output"))

PHIS = [0.50, 0.55, 0.60, 0.64]
NODES = 4096
KNODES = 1024
TIMEOUT_SEC = 35

def run_solver(phi, init_gamma=None, timeout=TIMEOUT_SEC):
    cmd = [
        SOLVER_BIN,
        "--closure", "PY",
        "--potential", "7",
        "--volfactor", f"{phi:.4f}",
        "--temp", "1.0",
        "--nodes", str(NODES),
        "--knodes", str(KNODES),
        "--gamma"
    ]
    if init_gamma is not None and os.path.exists(init_gamma):
        cmd.extend(["--init-gamma", init_gamma])

    t0 = time.perf_counter()
    try:
        proc = subprocess.run(
            cmd,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True,
            timeout=timeout
        )
        elapsed = time.perf_counter() - t0
        converged = (proc.returncode == 0)
        return elapsed, converged, proc.stdout, proc.stderr
    except subprocess.TimeoutExpired:
        return timeout, False, "", "TIMEOUT"

def analytical_sk_scalar(k_val, phi, sigma=1.0):
    """Exact Wertheim-Thiele solution for Hard Sphere Percus-Yevick S(k) at scalar k."""
    eta = phi
    lambda1 = (1.0 + 2.0 * eta)**2 / (1.0 - eta)**4
    lambda2 = -1.5 * eta * (1.0 + eta / 2.0)**2 / (1.0 - eta)**4
    
    k_sig = max(abs(k_val * sigma), 1e-9)
    s = math.sin(k_sig)
    c = math.cos(k_sig)
    
    alpha = lambda1
    beta = 6.0 * eta * lambda2
    gamma_c = 0.5 * eta * lambda1
    
    term1 = alpha * (s - k_sig * c) / (k_sig**3)
    term2 = beta * (2.0 * k_sig * s - (k_sig**2 - 2.0) * c - 2.0) / (k_sig**4)
    term3 = gamma_c * ((4.0 * k_sig**3 - 24.0 * k_sig) * s - (k_sig**4 - 12.0 * k_sig**2 + 24.0) * c + 24.0) / (k_sig**6)
    
    rho = 6.0 * eta / (math.pi * sigma**3)
    c_k = -4.0 * math.pi * sigma**3 * (term1 + term2 + term3)
    
    return 1.0 / (1.0 - rho * c_k)

def main():
    print("=" * 70)
    print("OZE-c-Solver: Cold-Start vs Warm-Start Benchmark Harness")
    print(f"Discretization: {NODES} real-space nodes, {KNODES} Fourier nodes")
    print("=" * 70)

    os.makedirs(DATA_DIR, exist_ok=True)
    
    seed_050 = os.path.join(DATA_DIR, "gamma_phi_0.50.dat")
    if not os.path.exists(seed_050):
        print("[Setup] Generating baseline gamma for phi=0.50...")
        run_solver(0.50)
        subprocess.run(["cp", f"{OUTPUT_DIR}/PY_GammaDeR.dat", seed_050], check=True)

    results = []
    prior_seed = seed_050

    for phi in PHIS:
        print(f"\n---> Benchmarking phi = {phi:.2f} <---")
        
        current_seed = seed_050 if phi == 0.50 else prior_seed

        # 1. Benchmark Cold-Start (2 iterations)
        cold_times = []
        cold_converged = True
        print("  Running Cold-Start runs...")
        for it in range(2):
            t, conv, out, err = run_solver(phi, init_gamma=None, timeout=TIMEOUT_SEC)
            if conv:
                cold_times.append(t)
                print(f"    Cold run {it+1}: {t:.3f} s (converged)")
            else:
                cold_converged = False
                cold_times.append(TIMEOUT_SEC)
                print(f"    Cold run {it+1}: FAILED / TIMED OUT ({TIMEOUT_SEC} s)")
                break

        if cold_converged:
            cold_mean = statistics.mean(cold_times)
            cold_std = statistics.stdev(cold_times) if len(cold_times) > 1 else 0.0
        else:
            cold_mean = float("nan")
            cold_std = float("nan")

        # Save cold outputs if converged
        if cold_converged and phi < 0.64:
            subprocess.run(["cp", f"{OUTPUT_DIR}/PY_GdeR.dat", f"{DATA_DIR}/gr_cold_phi_{phi:.2f}.dat"], check=True)
            subprocess.run(["cp", f"{OUTPUT_DIR}/PY_SdeK.dat", f"{DATA_DIR}/sk_cold_phi_{phi:.2f}.dat"], check=True)

        # 2. Benchmark Warm-Start (2 iterations)
        warm_times = []
        warm_converged = True
        print(f"  Running Warm-Start runs (seeded with {os.path.basename(current_seed)})...")
        for it in range(2):
            t, conv, out, err = run_solver(phi, init_gamma=current_seed, timeout=TIMEOUT_SEC)
            if conv:
                warm_times.append(t)
                print(f"    Warm run {it+1}: {t:.3f} s (converged)")
            else:
                warm_converged = False
                warm_times.append(TIMEOUT_SEC)
                print(f"    Warm run {it+1}: FAILED / TIMED OUT")
                break

        if warm_converged:
            warm_mean = statistics.mean(warm_times)
            warm_std = statistics.stdev(warm_times) if len(warm_times) > 1 else 0.0
        else:
            warm_mean = float("nan")
            warm_std = float("nan")

        # Save warm outputs and prepare next seed
        if warm_converged:
            next_seed = os.path.join(DATA_DIR, f"gamma_phi_{phi:.2f}.dat")
            subprocess.run(["cp", f"{OUTPUT_DIR}/PY_GammaDeR.dat", next_seed], check=True)
            subprocess.run(["cp", f"{OUTPUT_DIR}/PY_GdeR.dat", f"{DATA_DIR}/gr_warm_phi_{phi:.2f}.dat"], check=True)
            subprocess.run(["cp", f"{OUTPUT_DIR}/PY_SdeK.dat", f"{DATA_DIR}/sk_warm_phi_{phi:.2f}.dat"], check=True)
            prior_seed = next_seed

        speedup = (cold_mean / warm_mean) if (cold_converged and warm_converged) else float("inf")

        results.append({
            "phi": phi,
            "cold_mean": cold_mean,
            "cold_std": cold_std,
            "cold_conv": cold_converged,
            "warm_mean": warm_mean,
            "warm_std": warm_std,
            "warm_conv": warm_converged,
            "speedup": speedup
        })

    # Write summary table to file for gnuplot and LaTeX
    summary_path = os.path.join(DATA_DIR, "timing_summary.dat")
    with open(summary_path, "w") as f:
        f.write("# phi cold_mean cold_std warm_mean warm_std speedup cold_status warm_status\n")
        for r in results:
            c_m = f"{r['cold_mean']:.4f}" if r['cold_conv'] else "35.0000"
            c_s = f"{r['cold_std']:.4f}" if r['cold_conv'] else "0.0000"
            w_m = f"{r['warm_mean']:.4f}" if r['warm_conv'] else "35.0000"
            w_s = f"{r['warm_std']:.4f}" if r['warm_conv'] else "0.0000"
            sp = f"{r['speedup']:.2f}" if r['cold_conv'] else "inf"
            c_stat = "CONV" if r['cold_conv'] else "DIVERGED"
            w_stat = "CONV" if r['warm_conv'] else "FAILED"
            f.write(f"{r['phi']:.2f} {c_m} {c_s} {w_m} {w_s} {sp} {c_stat} {w_stat}\n")

    # Generate analytical Wertheim S(k) curves for comparison
    for phi in PHIS:
        sk_warm_file = os.path.join(DATA_DIR, f"sk_warm_phi_{phi:.2f}.dat")
        if os.path.exists(sk_warm_file):
            k_pts = []
            with open(sk_warm_file, "r") as f_in:
                for line in f_in:
                    parts = line.strip().split()
                    if len(parts) >= 2:
                        try:
                            k_pts.append(float(parts[0]))
                        except ValueError:
                            pass
            
            ana_file = os.path.join(DATA_DIR, f"sk_analytical_phi_{phi:.2f}.dat")
            with open(ana_file, "w") as f_out:
                f_out.write("# k S_analytical(k)\n")
                for k in k_pts:
                    val = analytical_sk_scalar(k, phi, sigma=1.0)
                    f_out.write(f"{k:.8f} {val:.8f}\n")

    print("\n" + "=" * 70)
    print("BENCHMARK SUMMARY")
    print("=" * 70)
    print(f"{'phi':<6} | {'Cold-Start (s)':<18} | {'Warm-Start (s)':<18} | {'Speedup':<10} | {'Status'}")
    print("-" * 70)
    for r in results:
        c_str = f"{r['cold_mean']:.2f} ± {r['cold_std']:.2f}" if r['cold_conv'] else "DIVERGED (>35s)"
        w_str = f"{r['warm_mean']:.2f} ± {r['warm_std']:.2f}" if r['warm_conv'] else "FAILED"
        sp_str = f"{r['speedup']:.2f}x" if r['cold_conv'] else "N/A (Critical)"
        status = "Success" if r['warm_conv'] else "Failed"
        print(f"{r['phi']:<6.2f} | {c_str:<18} | {w_str:<18} | {sp_str:<10} | {status}")
    print("=" * 70)
    print(f"Summary written to: {summary_path}")

if __name__ == "__main__":
    main()
