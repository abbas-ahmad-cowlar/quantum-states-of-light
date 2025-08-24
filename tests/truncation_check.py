"""
Truncation Convergence Check
Verify that Hilbert space dimension N is large enough for accurate results.
"""
import numpy as np
import qutip
import sys
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parents[1]

SRC_DIR = PROJECT_ROOT / "src"
if str(SRC_DIR) not in sys.path:
    sys.path.insert(0, str(SRC_DIR))
from qo_utils import (
    photon_distribution,
    mean_photon_number,
    compute_g2_zero,
    coherent_tail,
)

alpha = 5.0  # Coherent state with <n> = 25

print("=== Truncation Convergence Test ===\n")
print(f"Coherent state |alpha={alpha}>, theoretical <n> = {alpha**2:.0f}\n")
print(f"{'N':>5} {'<n>':>10} {'P(N-1)':>12} {'tail':>12} {'g2(0)':>10} {'Converged?':>12}")
print("-" * 75)

for N in [15, 20, 25, 30, 40, 50, 60]:
    a = qutip.destroy(N)
    state = qutip.coherent(N, alpha)

    mn = mean_photon_number(state, a)
    P = photon_distribution(state)
    g2 = compute_g2_zero(state, a)
    tail = coherent_tail(N, alpha)

    converged = "YES" if abs(mn - alpha**2) < 0.1 and abs(g2 - 1.0) < 1e-3 and tail < 1e-6 else "NO"
    print(f"{N:>5} {mn:>10.4f} {P[-1]:>12.2e} {tail:>12.2e} {g2:>10.6f} {converged:>12}")

print(f"\nRule of thumb: N >= <n> + 5*sqrt(<n>) = {alpha**2 + 5*alpha:.0f}")
