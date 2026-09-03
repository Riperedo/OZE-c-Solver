import matplotlib.pyplot as plt
import numpy as np

mc = np.loadtxt("fig2a_rho_0.8_mu2_2.0.dat")
rhnc = np.loadtxt("hdr_rho_0.8_mu2_2.0_RHNC.dat")
lhnc = np.loadtxt("hdr_rho_0.8_mu2_2.0_LHNC.dat")

plt.figure(figsize=(8,6))
plt.plot(mc[:,0], mc[:,1], 'ko', label="Monte Carlo")
plt.plot(rhnc[:,0], rhnc[:,1]+1.0, 'r-', linewidth=2, label="Exact RHNC")
plt.plot(lhnc[:,0], lhnc[:,1]+1.0, 'b--', linewidth=1.5, label="LHNC")
plt.xlim(1.0, 1.6)
plt.ylim(0, 5)
plt.xlabel("r/sigma")
plt.ylabel("g000(r)")
plt.legend()
plt.title("Fries & Patey Validation (rho=0.8, mu2=2.0)")
plt.savefig("Fig2a_RHNC_Validation.png")
