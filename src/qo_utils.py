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

