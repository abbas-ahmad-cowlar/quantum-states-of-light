#!/usr/bin/env python
"""verify_notebooks.py — Verify all notebooks, figures, and physics sentinels."""
import sys
from pathlib import Path

import numpy as np
import qutip

HERE = Path(__file__).resolve()
PROJECT_ROOT = HERE.parent
if not (PROJECT_ROOT / "notebooks").exists() and (PROJECT_ROOT.parent / "notebooks").exists():
    PROJECT_ROOT = PROJECT_ROOT.parent

SRC_DIR = PROJECT_ROOT / "src"
if str(SRC_DIR) not in sys.path:
    sys.path.insert(0, str(SRC_DIR))

from qo_utils import compute_g2_zero, mean_photon_number, mandel_Q, wigner_normalization

FIG_DIR = PROJECT_ROOT / "figures"

EXPECTED_FIGURE_BASES = [
    "fock_photon_distributions", "fock_wigner_n0_n1_n3",
    "coherent_photon_distributions", "coherent_poisson_overlay",
    "coherent_wigner_phase_space", "coherent_wigner_sweep",
    "fock_vs_coherent_comparison",
    "squeezed_photon_distributions", "squeezed_quadrature_scaling",
    "squeezed_uncertainty_ellipses", "squeezed_wigner_progression",
    "squeezed_coherent_wigner",
    "all_states_photon_distributions", "g2_bar_chart_comparison",
    "all_states_wigner_comparison", "mandel_Q_vs_g2_scatter",
    "hom_effect_and_dip", "bs_output_probabilities",
    "tmsv_photon_correlations", "bell_chsh_violation",
    "bell_correlation_function", "bell_S_vs_concurrence",
    "wigner_gallery_6panel",
]

results = {}


def assert_close(name, value, expected, atol):
    if not np.isclose(value, expected, atol=atol, rtol=0):
        raise AssertionError(f"{name}: got {value}, expected {expected} +/- {atol}")


def run_physics_sentinels():
    """Independent numerical sentinels for central physics claims."""
    N = 80
    a = qutip.destroy(N)

    assert_close("g2(|1>)", compute_g2_zero(qutip.basis(N, 1), a), 0.0, 1e-10)

    alpha = np.sqrt(3)
    coherent = qutip.coherent(N, alpha)
    assert_close("coherent mean", mean_photon_number(coherent, a), 3.0, 1e-3)
    assert_close("coherent g2", compute_g2_zero(coherent, a), 1.0, 1e-3)

    thermal = qutip.thermal_dm(N, 3.0, method="operator")
    assert_close("thermal g2", compute_g2_zero(thermal, a), 2.0, 5e-3)
    assert_close("thermal Q", mandel_Q(thermal, a), 3.0, 5e-3)

    r = np.arcsinh(np.sqrt(3))
    squeezed = qutip.squeeze(N, r) * qutip.basis(N, 0)
    assert_close("squeezed mean", mean_photon_number(squeezed, a), 3.0, 1e-3)
    assert_close("squeezed g2", compute_g2_zero(squeezed, a), 3.0 + 1.0/3.0, 5e-3)

    xvec = np.linspace(-6, 6, 301)
    W = qutip.wigner(qutip.basis(N, 1), xvec, xvec)
    assert_close("Wigner norm |1>", wigner_normalization(W, xvec, xvec), 1.0, 5e-3)

    # Bell/CHSH sentinel
    sx, sz = qutip.sigmax(), qutip.sigmaz()
    def sigma_pol(phi):
        return np.cos(2*phi) * sz + np.sin(2*phi) * sx
    def bell_E(psi, a, b):
        return float(np.real(qutip.expect(qutip.tensor(sigma_pol(a), sigma_pol(b)), psi)))

    psi_minus = (qutip.tensor(qutip.basis(2,0), qutip.basis(2,1))
                 - qutip.tensor(qutip.basis(2,1), qutip.basis(2,0))).unit()
    S = abs(bell_E(psi_minus,0,np.pi/8) - bell_E(psi_minus,0,3*np.pi/8)
            + bell_E(psi_minus,np.pi/4,np.pi/8) + bell_E(psi_minus,np.pi/4,3*np.pi/8))
    assert_close("CHSH singlet", S, 2*np.sqrt(2), 1e-10)

    # HOM sentinel
    N_bs = 4
    a1 = qutip.tensor(qutip.destroy(N_bs), qutip.qeye(N_bs))
    a2 = qutip.tensor(qutip.qeye(N_bs), qutip.destroy(N_bs))
    U = (np.pi/4 * (a1.dag()*a2 - a1*a2.dag())).expm()
    psi_out = U * qutip.tensor(qutip.basis(N_bs,1), qutip.basis(N_bs,1))
    P11 = abs(qutip.tensor(qutip.basis(N_bs,1), qutip.basis(N_bs,1)).overlap(psi_out))**2
    assert_close("HOM P(1,1)", P11, 0.0, 1e-10)

    print("  [PASS] All physics sentinels passed")


# Check figures
print("=== Figure Manifest Check ===")
for fig_base in EXPECTED_FIGURE_BASES:
    for ext in [".png", ".svg"]:
        fig = f"{fig_base}{ext}"
        path = FIG_DIR / fig
        if not path.exists():
            results[f"figure:{fig}"] = "FAIL"
            print(f"  MISSING: {fig}")
        elif path.stat().st_size < 10_000:
            results[f"figure:{fig}"] = "FAIL"
            print(f"  TOO SMALL: {fig}")
        else:
            results[f"figure:{fig}"] = "PASS"

# Check required files
print("\n=== Required Files ===")
for artifact in ["requirements.txt", "requirements-lock.txt", "environment_report.txt"]:
    path = PROJECT_ROOT / artifact
    if path.exists() and path.stat().st_size > 0:
        results[f"file:{artifact}"] = "PASS"
        print(f"  OK: {artifact}")
    else:
        results[f"file:{artifact}"] = "FAIL"
        print(f"  MISSING/EMPTY: {artifact}")

# Physics sentinels
print("\n=== Physics Sentinels ===")
try:
    run_physics_sentinels()
    results["physics_sentinels"] = "PASS"
except Exception as exc:
    results["physics_sentinels"] = "FAIL"
    print(f"  FAIL: {exc}")

# Summary
print(f"\n{'='*50}")
print("SUMMARY")
print(f"{'='*50}")
passed = sum(1 for s in results.values() if s == "PASS")
failed = sum(1 for s in results.values() if s == "FAIL")
total = len(results)
print(f"  {passed} passed, {failed} failed / {total} total")
if failed > 0:
    for k, v in results.items():
        if v == "FAIL":
            print(f"  FAIL: {k}")
    sys.exit(1)
else:
    print("  ALL STATIC CHECKS PASS (figures, files, physics sentinels).")
    print("  Note: This does not execute notebooks. Run 'jupyter nbconvert --execute' separately.")
