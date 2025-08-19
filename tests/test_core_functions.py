"""
Unit tests for core physics functions in qo_utils.py.
Tests mean_photon_number, photon_variance, mandel_Q, compute_g2_zero,
and photon_distribution against known analytical results.

Uses N=40 for coherent alpha=3 (mean_n=9) to avoid truncation noise.
"""
import sys
from pathlib import Path
import math
import numpy as np
import pytest
import qutip

PROJECT_ROOT = Path(__file__).resolve().parents[1]
SRC_DIR = PROJECT_ROOT / "src"
if str(SRC_DIR) not in sys.path:
    sys.path.insert(0, str(SRC_DIR))
from qo_utils import (
    mean_photon_number, photon_variance, mandel_Q,
    compute_g2_zero, photon_distribution,
    coherent_tail, thermal_tail, squeezed_vacuum_tail,
)

N = 40
a = qutip.destroy(N)


# === mean_photon_number ===

def test_mean_photon_vacuum():
    assert np.isclose(mean_photon_number(qutip.basis(N, 0), a), 0.0)

def test_mean_photon_fock5():
    assert np.isclose(mean_photon_number(qutip.basis(N, 5), a), 5.0)

def test_mean_photon_coherent():
    assert np.isclose(mean_photon_number(qutip.coherent(N, 3.0), a), 9.0, atol=0.1)

def test_mean_photon_thermal():
    assert np.isclose(mean_photon_number(qutip.thermal_dm(N, 2.0), a), 2.0, atol=0.1)


# === photon_variance ===

def test_variance_fock():
    assert np.isclose(photon_variance(qutip.basis(N, 3), a), 0.0, atol=1e-10)

def test_variance_coherent():
    assert np.isclose(photon_variance(qutip.coherent(N, 3.0), a), 9.0, atol=0.1)


# === mandel_Q ===

def test_mandel_Q_fock():
    assert np.isclose(mandel_Q(qutip.basis(N, 3), a), -1.0, atol=1e-10)

def test_mandel_Q_coherent():
    assert np.isclose(mandel_Q(qutip.coherent(N, 3.0), a), 0.0, atol=1e-4)

def test_mandel_Q_thermal_positive():
    Q_th = mandel_Q(qutip.thermal_dm(N, 3.0), a)
    assert Q_th > 0, f"Q should be positive for thermal, got {Q_th}"


# === g2/Q cross-check ===

def test_g2_Q_crosscheck():
    state = qutip.basis(N, 3)
    g2 = compute_g2_zero(state, a)
    Q = mandel_Q(state, a)
    mn = mean_photon_number(state, a)
    assert np.isclose(g2, 1 + Q / mn, atol=1e-10), "g2(0) != 1 + Q/<n>"


# === compute_g2_zero ===

def test_g2_single_photon():
    g2 = compute_g2_zero(qutip.basis(N, 1), a)
    assert np.isclose(g2, 0.0, atol=1e-10)

def test_g2_fock3():
    g2 = compute_g2_zero(qutip.basis(N, 3), a)
    assert np.isclose(g2, 1 - 1 / 3, atol=1e-10)

def test_g2_coherent():
    g2 = compute_g2_zero(qutip.coherent(N, 3.0), a)
    assert np.isclose(g2, 1.0, atol=1e-6)

def test_g2_thermal():
    g2 = compute_g2_zero(qutip.thermal_dm(N, 3.0), a)
    assert np.isclose(g2, 2.0, atol=0.01)

def test_g2_vacuum_nan():
    g2 = compute_g2_zero(qutip.basis(N, 0), a)
    assert np.isnan(g2)


# === photon_distribution ===

def test_photon_dist_fock():
    P = photon_distribution(qutip.basis(N, 3))
    assert np.isclose(P[3], 1.0)
    assert np.isclose(P.sum(), 1.0)
    assert np.isclose(P[0], 0.0)

def test_photon_dist_coherent_poisson():
    state = qutip.coherent(N, 2.0)
    P = photon_distribution(state)
    assert np.isclose(P.sum(), 1.0, atol=1e-12)
    alpha_sq = 4.0
    for n in range(10):
        expected = np.exp(-alpha_sq) * alpha_sq**n / math.factorial(n)
        assert np.isclose(P[n], expected, atol=1e-6), f"P({n}) mismatch"

def test_photon_dist_thermal_normalized():
    P = photon_distribution(qutip.thermal_dm(N, 2.0))
    assert np.isclose(P.sum(), 1.0, atol=1e-12)


# === tail helpers (Finding 4 fix) ===

def test_coherent_tail_small():
    """Poisson tail for alpha=3, N=40 should be negligible."""
    assert coherent_tail(40, 3.0) < 1e-8

def test_coherent_tail_large():
    """Poisson tail for alpha=5, N=15 should be large."""
    assert coherent_tail(15, 5.0) > 0.5

def test_thermal_tail_basic():
    """Geometric tail P(n>=N) = (n_bar/(1+n_bar))^N."""
    n_bar = 2.0
    N_test = 20
    expected = (n_bar / (1 + n_bar)) ** N_test
    assert np.isclose(thermal_tail(N_test, n_bar), expected)

