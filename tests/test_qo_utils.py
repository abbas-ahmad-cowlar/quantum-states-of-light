"""
Integration test for qo_utils.py.
Tests the full utility module with a coherent state as the reference.
"""
import sys
from pathlib import Path
import numpy as np
import pytest
import qutip
import matplotlib
matplotlib.use('Agg')

PROJECT_ROOT = Path(__file__).resolve().parents[1]
SRC_DIR = PROJECT_ROOT / "src"
if str(SRC_DIR) not in sys.path:
    sys.path.insert(0, str(SRC_DIR))
from qo_utils import (
    photon_distribution,
    wigner_normalization,
    coherent_tail,
    cutoff_from_tail,
    squeezed_wigner_extent,
    wigner_sign_summary,
    compute_g2_zero,
    mean_photon_number,
    photon_variance,
    mandel_Q,
    plot_photon_distribution,
    plot_wigner,
)


# Shared fixtures
N = 40
a = qutip.destroy(N)
alpha = 3.0
state = qutip.coherent(N, alpha)


def test_photon_distribution_normalized():
    P = photon_distribution(state)
    assert np.isclose(P.sum(), 1.0, atol=1e-12)

def test_coherent_tail_negligible():
    assert coherent_tail(N, alpha) < 1e-8

def test_boundary_probability_tiny():
    P = photon_distribution(state)
    assert P[-1] < 1e-8

def test_mean_photon_number():
    mn = mean_photon_number(state, a)
    assert np.isclose(mn, 9.0, atol=0.1)

def test_photon_variance():
    var = photon_variance(state, a)
    assert np.isclose(var, 9.0, atol=0.2)

def test_g2_zero_coherent():
    g2 = compute_g2_zero(state, a)
    assert np.isclose(g2, 1.0, atol=1e-4)

def test_mandel_Q_coherent():
    Q = mandel_Q(state, a)
    assert np.isclose(Q, 0.0, atol=1e-3)

def test_g2_Q_crosscheck():
    g2 = compute_g2_zero(state, a)
    Q = mandel_Q(state, a)
    mn = mean_photon_number(state, a)
    assert np.isclose(g2, 1 + Q / mn, atol=1e-6)

