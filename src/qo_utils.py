"""
Quantum Optics Utility Module
==============================

Shared functions for the Quantum States of Light Explorer project.
All functions operate on QuTiP quantum objects (Qobj).

Functions:
    1. photon_distribution(state, n_max) -- Compute P(n) = |<n|state>|^2
    2. plot_photon_distribution(...)      -- Bar chart of P(n)
    3. plot_wigner(...)                   -- 2D colormap of Wigner function
    3b. wigner_normalization(W, xvec, yvec) -- Numerically verify integral W dx dp ~ 1
    4. compute_g2_zero(state, a)          -- Compute g^(2)(0) = <a_dag a_dag a a>/<a_dag a>^2
    5. mean_photon_number(state, a)       -- Compute <n_hat> = <a_dag a>
    6. photon_variance(state, a)          -- Compute Dn^2 = <n_hat^2> - <n_hat>^2
    7. mandel_Q(state, a)                -- Compute Q = (Dn^2 - <n>)/<n>

"""

import numpy as np
import matplotlib.pyplot as plt
import qutip
import math
from scipy.stats import poisson


# =========================================================================
# Function 5: mean_photon_number
# =========================================================================

def mean_photon_number(state, a):
    """
    Compute the mean photon number <n_hat> = <a_dag a>.

    Physics: The average number of photons in the state.
    - Fock |n>: <n_hat> = n (exactly n photons)
    - Coherent |alpha>: <n_hat> = |alpha|^2 (proportional to intensity)
    - Thermal (n_bar): <n_hat> = n_bar = 1/(exp(hbar*omega/kT) - 1)

    Parameters
    ----------
    state : qutip.Qobj
        Quantum state (ket or density matrix).
    a : qutip.Qobj
        Annihilation operator.

    Returns
    -------
    float
        The mean photon number <n_hat>.
    """
    return np.real(qutip.expect(a.dag() * a, state))


# =========================================================================
# Function 6: photon_variance
# =========================================================================

def photon_variance(state, a):
    """
    Compute the photon number variance Dn^2 = <n_hat^2> - <n_hat>^2.

    Physics: Measures the spread in photon number.
    - Fock |n>: Dn^2 = 0 (exactly n photons, zero uncertainty)
    - Coherent |alpha>: Dn^2 = |alpha|^2 = <n> (Poissonian: variance = mean)
    - Thermal: Dn^2 = <n>(1 + <n>) (super-Poissonian: variance > mean)

    Parameters
    ----------
    state : qutip.Qobj
        Quantum state (ket or density matrix).
    a : qutip.Qobj
        Annihilation operator.

    Returns
    -------
    float
        The photon number variance Dn^2.
    """
    n_op = a.dag() * a
    mean_n = np.real(qutip.expect(n_op, state))
    mean_n_sq = np.real(qutip.expect(n_op * n_op, state))
    return mean_n_sq - mean_n**2


# =========================================================================
# Function 7: mandel_Q
# =========================================================================

def mandel_Q(state, a):
    """
    Compute the Mandel Q parameter: Q = (Dn^2 - <n>) / <n>.

    Physics: Measures departure from Poissonian statistics.
    - Q < 0: sub-Poissonian (non-classical, anti-bunching)
    - Q = 0: Poissonian (coherent state, boundary)
    - Q > 0: super-Poissonian (bunching, thermal-like)

    Related to g^(2)(0) by: g^(2)(0) = 1 + Q/<n>

    Parameters
    ----------
    state : qutip.Qobj
        Quantum state (ket or density matrix).
    a : qutip.Qobj
        Annihilation operator.

    Returns
    -------
    float
        The Mandel Q parameter. Returns np.nan for vacuum (<n> = 0).
    """
    mean_n = mean_photon_number(state, a)

    if mean_n < 1e-10:
        return np.nan

    var_n = photon_variance(state, a)
    return (var_n - mean_n) / mean_n


# =========================================================================
# Function 4: compute_g2_zero
# =========================================================================

def compute_g2_zero(state, a):
    """
    Compute the second-order coherence function g^(2)(0).

    Physics: g^(2)(0) = <a_dag a_dag a a> / <a_dag a>^2 = <n_hat(n_hat-1)> / <n_hat>^2

    This is a central zero-delay photon-correlation diagnostic:
    - g^(2)(0) = 2    -> nonzero thermal state (bunching, super-Poissonian)
    - g^(2)(0) = 1    -> Poissonian second-order statistics
    - g^(2)(0) < 1    -> anti-bunched/sub-Poissonian
    - g^(2)(0) = 0    -> single photon source |1>

    It is not a complete state classifier; interpret it with P(n), Wigner
    functions, quadrature variances, and truncation diagnostics.

    Parameters
    ----------
    state : qutip.Qobj
        Quantum state (ket or density matrix).
    a : qutip.Qobj
        Annihilation operator for the mode.

    Returns
    -------
    g2 : float
        The g^(2)(0) value. Returns np.nan for vacuum state (<n> = 0).

    Notes
    -----
    CRITICAL: The numerator is <a_dag a_dag a a>, NOT <(a_dag a)^2> = <n_hat^2>.
    The difference: <a_dag a_dag a a> = <n_hat^2 - n_hat> = <n_hat(n_hat-1)>.
    Using (a.dag()*a)**2 instead gives WRONG results.
    """
    # Compute <n_hat> = <a_dag a>
    mean_n = qutip.expect(a.dag() * a, state)

    # Guard: g^(2)(0) is undefined for vacuum (division by zero)
    if mean_n < 1e-10:
        return np.nan

    # Compute <a_dag a_dag a a> -- NOTE THE OPERATOR ORDERING!
    # This is a.dag() * a.dag() * a * a, NOT (a.dag() * a)**2
    numerator = qutip.expect(a.dag() * a.dag() * a * a, state)

    g2 = numerator / mean_n**2

    return np.real(g2)  # Should be real; discard floating-point imaginary noise


# =========================================================================
# Function 1: photon_distribution
# =========================================================================

