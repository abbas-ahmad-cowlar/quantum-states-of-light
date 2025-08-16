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

