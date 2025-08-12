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

def photon_distribution(state, n_max=None):
    """
    Compute the photon number distribution P(n) = |<n|psi>|^2 for a quantum state.

    Physics: P(n) is the probability of detecting exactly n photons.
    - For Fock states: P(n) is a delta function (all probability on one number)
    - For coherent states: P(n) is a Poisson distribution
    - For thermal states: P(n) is a Bose-Einstein (geometric) distribution

    Parameters
    ----------
    state : qutip.Qobj
        Quantum state (ket vector or density matrix).
    n_max : int, optional
        Number of probabilities to return. If omitted, uses the state's Hilbert
        space dimension.

    Returns
    -------
    probabilities : numpy.ndarray
        Array of length n_max with P(n) for n = 0, 1, ..., n_max-1.

    Raises
    ------
    ValueError
        If state is not a valid quantum object (ket or density matrix).
    """
    if not (state.isket or state.isoper):
        raise ValueError(f"Expected ket or density matrix, got {state.type}")
    if len(state.dims[0]) != 1:
        raise ValueError("photon_distribution is single-mode only.")

    dim = state.shape[0]
    if n_max is None:
        n_max = dim
    if n_max > dim:
        raise ValueError(f"n_max={n_max} exceeds state dimension {dim}")

    if state.isket:
        probabilities = np.abs(state.full().ravel())**2
    else:
        probabilities = np.real(np.diag(state.full()))

    probabilities = np.array(probabilities[:n_max], dtype=float, copy=True)
    probabilities[np.isclose(probabilities, 0.0, atol=1e-15)] = 0.0
    return probabilities


# =========================================================================
# Function 2: plot_photon_distribution
# =========================================================================

def plot_photon_distribution(state, n_max=None, title="", ax=None, color='steelblue',
                              show_stats=True):
    """
    Plot the photon number distribution P(n) as a bar chart.

    Physics: Visualizes the probability of finding n photons in the state.
    The shape of this distribution classifies the type of light:
    - Delta function -> Fock state
    - Poissonian (bell) -> Coherent state (laser)
    - Geometric (exponential decay) -> Thermal state
    - Even-only bars -> Squeezed vacuum

    Parameters
    ----------
    state : qutip.Qobj
        Quantum state (ket or density matrix).
    n_max : int, optional
        Number of Fock states to plot. If omitted, uses the state dimension.
    title : str
        Plot title (supports LaTeX: r"$|\\alpha=3\\rangle$").
    ax : matplotlib.axes.Axes, optional
        Axes to plot on. If None, creates a new figure.
    color : str
        Bar color. If omitted, uses 'steelblue'.
    show_stats : bool
        If True, annotate with <n> and Dn on the plot.

    Returns
    -------
    ax : matplotlib.axes.Axes
        The axes object (for further customization or multi-panel figures).
    """
    if len(state.dims[0]) != 1:
        raise ValueError("plot_photon_distribution is single-mode only.")

    dim = state.shape[0]
    if n_max is None:
        n_max = dim
    if n_max > dim:
        raise ValueError(f"n_max={n_max} exceeds state dimension {dim}")

    P = photon_distribution(state, n_max=n_max)
    n_values = np.arange(n_max)

    if ax is None:
        fig, ax = plt.subplots(figsize=(8, 4))

    ax.bar(n_values, P, color=color, alpha=0.8, edgecolor='black', linewidth=0.5)
    ax.set_xlabel(r'Photon number $n$', fontsize=12)
    ax.set_ylabel(r'$P(n)$', fontsize=12)
    ax.set_title(title, fontsize=14)
    ax.set_xlim(-0.5, n_max - 0.5)
    ax.set_ylim(0, max(P) * 1.15 if len(P) and max(P) > 0 else 1)

    if show_stats:
        a = qutip.destroy(dim)
        mn = mean_photon_number(state, a)
        var_n = photon_variance(state, a)
        delta_n = np.sqrt(max(var_n, 0))  # Guard against tiny negative from numerics
        ax.text(0.95, 0.95,
                f'$\\langle n \\rangle = {mn:.2f}$\n$\\Delta n = {delta_n:.2f}$',
                transform=ax.transAxes, fontsize=10,
                verticalalignment='top', horizontalalignment='right',
                bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.8))

    return ax


# =========================================================================
# Function 3: plot_wigner + wigner_normalization
# =========================================================================

def _trapz_compat(y, x, axis):
    """
    NumPy 1.26/2.x compatible trapezoidal integration.

    NumPy 2.x provides np.trapezoid; NumPy 1.x uses np.trapz. The project
    allows NumPy >=1.26, so the helper must support both.
    """
    if hasattr(np, "trapezoid"):
        return np.trapezoid(y, x, axis=axis)
    return np.trapz(y, x, axis=axis)


def wigner_normalization(W, xvec, yvec):
    """
    Numerically integrate a Wigner function over the plotted phase-space grid.

    Returns integral W(x,p) dx dp. For a sufficiently large grid this should
    be close to 1.0. Uses a compatibility helper so NumPy 1.26 and 2.x both work.
    """
    W = np.asarray(W)
    if W.shape == (len(yvec), len(xvec)):
        # Standard plotting convention: rows follow y/p, columns follow x.
        return _trapz_compat(_trapz_compat(W, xvec, axis=1), yvec, axis=0)
    if W.shape == (len(xvec), len(yvec)):
        # Fallback for transposed arrays on non-square grids.
        return _trapz_compat(_trapz_compat(W, yvec, axis=1), xvec, axis=0)
    raise ValueError(
        f"W shape {W.shape} is incompatible with xvec={len(xvec)}, yvec={len(yvec)}"
    )


def plot_wigner(state, xvec=None, yvec=None, title="", ax=None,
                colormap='RdBu_r', show_colorbar=True):
    """
    Plot the Wigner function W(x, p) of a quantum state as a 2D colormap.
    This project uses QuTiP's standard g=sqrt(2) convention, so a coherent
    state |alpha> is centered at (sqrt(2) Re alpha, sqrt(2) Im alpha).

    Physics: The Wigner function is a quasi-probability distribution in
    phase space (x, p coordinates). Key features:
    - Normalized: integral W(x,p) dx dp = 1
    - CAN be negative (non-classicality witness)
    - Positive Gaussian -> classical-like state (coherent)
    - Negative regions -> non-Gaussian nonclassicality (Fock, cat states)
    - Note: squeezed vacuum is nonclassical but Wigner-POSITIVE

    Parameters
    ----------
    state : qutip.Qobj
        Quantum state (ket or density matrix).
    xvec : numpy.ndarray, optional
        x-axis grid points. If omitted, uses np.linspace(-6, 6, 200).
    yvec : numpy.ndarray, optional
        p-axis grid points. If omitted, uses same values as xvec.
    title : str
        Plot title (supports LaTeX).
    ax : matplotlib.axes.Axes, optional
        Axes to plot on. If None, creates new figure.
    colormap : str
        Colormap name. If omitted, uses 'RdBu_r' (diverging with neutral midpoint).
    show_colorbar : bool
        Whether to add a colorbar.

    Returns
    -------
    ax : matplotlib.axes.Axes
        The axes object.
    W : numpy.ndarray
        The computed Wigner function values (for further analysis).
    """
    if xvec is None:
        xvec = np.linspace(-6, 6, 200)
    if yvec is None:
        yvec = xvec.copy()

    # Compute Wigner function using QuTiP
    W = qutip.wigner(state, xvec, yvec)

    if ax is None:
        fig, ax = plt.subplots(figsize=(6, 5))

    # Use explicit symmetric levels so zero is visually centered.
    wmax = np.max(np.abs(W))
    if wmax == 0:
        wmax = 1.0
    levels = np.linspace(-wmax, wmax, 101)

    im = ax.contourf(xvec, yvec, W, levels=levels, cmap=colormap, extend='both')

    if show_colorbar:
        plt.colorbar(im, ax=ax, label=r'$W(x, p)$')

    ax.set_xlabel(r'$x$ (position quadrature)', fontsize=12)
    ax.set_ylabel(r'$p$ (momentum quadrature)', fontsize=12)
    ax.set_title(title, fontsize=14)
    ax.set_aspect('equal')

    return ax, W


# =========================================================================
# Numerical Safety Helpers
# =========================================================================

def coherent_tail(N, alpha):
    """Analytic Poisson tail P(n >= N) for an ideal coherent state."""
    return float(poisson.sf(N - 1, abs(alpha) ** 2))


