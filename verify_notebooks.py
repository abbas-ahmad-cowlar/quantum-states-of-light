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


