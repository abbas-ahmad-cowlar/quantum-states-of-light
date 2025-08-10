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

