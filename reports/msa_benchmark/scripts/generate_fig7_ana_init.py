#!/usr/bin/env python3
"""
Generate Figure 7: Structure factor profiles under direct analytical initialization at T*=0.10 and T*=0.01.
"""

import os
import matplotlib.pyplot as plt
import numpy as np

os.makedirs("reports/msa_benchmark/plots", exist_ok=True)

plt.style.use('seaborn-v0_8-paper' if 'seaborn-v0_8-paper' in plt.style.available else 'default')
plt.rcParams.update({
    'font.size': 10,
    'axes.labelsize': 11,
    'axes.titlesize': 11,
    'legend.fontsize': 9,
    'xtick.labelsize': 9,
    'ytick.labelsize': 9,
    'figure.autolayout': True,
    'lines.linewidth': 1.6
})

fig, axs = plt.subplots(2, 2, figsize=(10, 8), dpi=300)

phis_to_plot = [0.1, 0.3, 0.5]
colors = {0.1: '#1f77b4', 0.3: '#2ca02c', 0.5: '#d62728'}

# Panel (a): S110 at T* = 0.10
ax = axs[0, 0]
for phi in phis_to_plot:
    ana_path = f"data/input_sk/sk_phi_{phi:.1f}_T_0.1.dat"
    num_path = f"reports/msa_benchmark/data/num_ana_init_phi_{phi:.1f}_T_0.10.dat"
    ana = np.loadtxt(ana_path)
    num = np.loadtxt(num_path)
    c = colors[phi]
    ax.plot(ana[:, 0], ana[:, 2], '--', color=c, label=f'Ana $\\phi={phi:.1f}$')
    ax.plot(num[:, 0], num[:, 2], '-', color=c, label=f'Num $\\phi={phi:.1f}$')
ax.set_xlim(0, 15)
ax.set_ylim(-0.5, 5.0)
ax.set_xlabel('$k\\sigma$')
ax.set_ylabel('$S^{110}(k)$')
ax.set_title('(a) $S^{110}(k)$ at $T^* = 0.10$ (Direct Warm-Start)')
ax.legend(ncol=2, frameon=True, fontsize=8)
ax.grid(True, linestyle=':', alpha=0.6)

# Panel (b): S112 at T* = 0.10
ax = axs[0, 1]
for phi in phis_to_plot:
    ana_path = f"data/input_sk/sk_phi_{phi:.1f}_T_0.1.dat"
    num_path = f"reports/msa_benchmark/data/num_ana_init_phi_{phi:.1f}_T_0.10.dat"
    ana = np.loadtxt(ana_path)
    num = np.loadtxt(num_path)
    c = colors[phi]
    ax.plot(ana[:, 0], ana[:, 3], '--', color=c, label=f'Ana $\\phi={phi:.1f}$')
    ax.plot(num[:, 0], num[:, 3], '-', color=c, label=f'Num $\\phi={phi:.1f}$')
ax.set_xlim(0, 15)
ax.set_ylim(-2.5, 0.5)
ax.set_xlabel('$k\\sigma$')
ax.set_ylabel('$S^{112}(k)$')
ax.set_title('(b) $S^{112}(k)$ at $T^* = 0.10$ (Direct Warm-Start)')
ax.legend(ncol=2, frameon=True, fontsize=8)
ax.grid(True, linestyle=':', alpha=0.6)

# Panel (c): S110 at T* = 0.01
ax = axs[1, 0]
for phi in phis_to_plot:
    ana_path = f"data/input_sk/sk_phi_{phi:.1f}_T_0.01.dat"
    num_path = f"reports/msa_benchmark/data/num_ana_init_phi_{phi:.1f}_T_0.01.dat"
    ana = np.loadtxt(ana_path)
    num = np.loadtxt(num_path)
    c = colors[phi]
    ax.plot(ana[:, 0], ana[:, 2], '--', color=c, label=f'Ana $\\phi={phi:.1f}$')
    ax.plot(num[:, 0], num[:, 2], '-', color=c, label=f'Num $\\phi={phi:.1f}$')
ax.set_xlim(0, 15)
ax.set_xlabel('$k\\sigma$')
ax.set_ylabel('$S^{110}(k)$')
ax.set_title('(c) $S^{110}(k)$ at $T^* = 0.01$ (Direct Warm-Start)')
ax.legend(ncol=2, frameon=True, fontsize=8)
ax.grid(True, linestyle=':', alpha=0.6)

# Panel (d): S112 at T* = 0.01
ax = axs[1, 1]
for phi in phis_to_plot:
    ana_path = f"data/input_sk/sk_phi_{phi:.1f}_T_0.01.dat"
    num_path = f"reports/msa_benchmark/data/num_ana_init_phi_{phi:.1f}_T_0.01.dat"
    ana = np.loadtxt(ana_path)
    num = np.loadtxt(num_path)
    c = colors[phi]
    ax.plot(ana[:, 0], ana[:, 3], '--', color=c, label=f'Ana $\\phi={phi:.1f}$')
    ax.plot(num[:, 0], num[:, 3], '-', color=c, label=f'Num $\\phi={phi:.1f}$')
ax.set_xlim(0, 15)
ax.set_xlabel('$k\\sigma$')
ax.set_ylabel('$S^{112}(k)$')
ax.set_title('(d) $S^{112}(k)$ at $T^* = 0.01$ (Direct Warm-Start)')
ax.legend(ncol=2, frameon=True, fontsize=8)
ax.grid(True, linestyle=':', alpha=0.6)

plt.tight_layout()
pdf_path = "reports/msa_benchmark/plots/fig7_ana_init_benchmark.pdf"
png_path = "reports/msa_benchmark/plots/fig7_ana_init_benchmark.png"
plt.savefig(pdf_path, format='pdf')
plt.savefig(png_path, format='png', dpi=300)
print(f"Saved: {pdf_path} and {png_path}")
