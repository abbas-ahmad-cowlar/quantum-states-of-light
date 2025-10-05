"""Generate Phase 5 deliverables: Wigner gallery + master stats table."""
from pathlib import Path
import sys

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import qutip

PROJECT_ROOT = Path(__file__).resolve().parent
SRC_DIR = PROJECT_ROOT / "src"
FIG_DIR = PROJECT_ROOT / "figures"
FIG_DIR.mkdir(parents=True, exist_ok=True)

if str(SRC_DIR) not in sys.path:
    sys.path.insert(0, str(SRC_DIR))

from qo_utils import (
    photon_distribution, mean_photon_number, photon_variance,
    compute_g2_zero, mandel_Q, coherent_tail, thermal_tail,
    squeezed_vacuum_tail, wigner_normalization, wigner_sign_summary,
)

plt.rcParams.update({
    'font.size': 12, 'axes.labelsize': 14, 'axes.titlesize': 16,
    'figure.dpi': 150, 'savefig.dpi': 300,
    'text.usetex': False, 'mathtext.fontset': 'stix', 'font.family': 'STIXGeneral',
})

# ── DELIVERABLE A: 6-Panel Wigner Gallery ──────────────────────────

print("=== Generating 6-Panel Wigner Gallery ===\n")

from matplotlib.gridspec import GridSpec

xvec = np.linspace(-8, 8, 200)
N = 60

# Row 1: Fock states
states_row1 = [
    (qutip.basis(N, 0), r'Vacuum $|0\rangle$'),
    (qutip.basis(N, 1), r'Fock $|1\rangle$'),
    (qutip.basis(N, 3), r'Fock $|3\rangle$'),
]

# Row 2: Other state types
alpha = 3.0
r_sq = 1.0
alpha_sc = 2.0
squeezed_coherent = (
    qutip.displace(N, alpha_sc)
    * qutip.squeeze(N, r_sq)
    * qutip.basis(N, 0)
)
states_row2 = [
    (qutip.coherent(N, alpha), rf'Coherent $|\alpha={alpha}\rangle$'),
    (qutip.squeeze(N, r_sq) * qutip.basis(N, 0), rf'Squeezed vacuum $r={r_sq}$'),
    (squeezed_coherent, rf'Squeezed coherent $D({alpha_sc})S({r_sq})|0\rangle$'),
]

all_states = states_row1 + states_row2

# First pass: compute Wigner data and find global color scale
W_data = []
for state, title in all_states:
    W = qutip.wigner(state, xvec, xvec)
    norm = wigner_normalization(W, xvec, xvec)
    assert np.isclose(norm, 1.0, atol=5e-2), f"Wigner integral for {title} = {norm}"
    W_data.append((W, title))
    print(f"  {title}: W integral = {norm:.4f}")

wmax_global = max(np.max(np.abs(W)) for W, _ in W_data)

# Use GridSpec with a dedicated colorbar column to avoid layout overlap
fig = plt.figure(figsize=(20, 11))
gs = GridSpec(2, 4, figure=fig, width_ratios=[1, 1, 1, 0.05],
              wspace=0.3, hspace=0.35)

axes = [[fig.add_subplot(gs[r, c]) for c in range(3)] for r in range(2)]
cbar_ax = fig.add_subplot(gs[:, 3])

for idx, (W, title) in enumerate(W_data):
    row, col = idx // 3, idx % 3
    ax = axes[row][col]

    levels = np.linspace(-wmax_global, wmax_global, 101)
    im = ax.contourf(xvec, xvec, W, levels=levels, cmap='RdBu_r', extend='both')
    ax.set_title(title, fontsize=14, fontweight='bold')
    ax.set_xlabel(r'$x$', fontsize=12)
    ax.set_ylabel(r'$p$', fontsize=12)
    ax.set_aspect('equal')

fig.colorbar(im, cax=cbar_ax, label=r'$W(x,p)$')

fig.suptitle('Wigner Functions of Fundamental Quantum States of Light',
             fontsize=20, fontweight='bold')
plt.savefig(FIG_DIR / 'wigner_gallery_6panel.png', bbox_inches='tight', dpi=300)
plt.savefig(FIG_DIR / 'wigner_gallery_6panel.svg', bbox_inches='tight')
plt.show()
print("\n[PASS] Wigner gallery saved.\n")

# ── DELIVERABLE B: Extended Master Statistics Table ─────────────────

print("=== Computing Extended Master Statistics Table ===\n")

nbar_max = 5.0
thermal_N = int(np.ceil(np.log(1e-6) / np.log(nbar_max / (nbar_max + 1)))) + 1
N = max(60, thermal_N)
a = qutip.destroy(N)
xvec_sign = np.linspace(-8, 8, 151)

master_states = [
    ('Vacuum |0>',      qutip.basis(N, 0)),
    ('Fock |1>',        qutip.basis(N, 1)),
    ('Fock |3>',        qutip.basis(N, 3)),
    ('Fock |5>',        qutip.basis(N, 5)),
    ('Coherent a=1',    qutip.coherent(N, 1.0)),
    ('Coherent a=3',    qutip.coherent(N, 3.0)),
    ('Coherent a=5',    qutip.coherent(N, 5.0)),
    ('Squeezed r=0.5',  qutip.squeeze(N, 0.5) * qutip.basis(N, 0)),
    ('Squeezed r=1.0',  qutip.squeeze(N, 1.0) * qutip.basis(N, 0)),
    ('Thermal nb=1',    qutip.thermal_dm(N, 1.0)),
    ('Thermal nb=3',    qutip.thermal_dm(N, 3.0)),
    ('Thermal nb=5',    qutip.thermal_dm(N, 5.0)),
]

# Truncation check
for name, state in master_states:
    P = photon_distribution(state, n_max=N)
    analytic_tail = None
    if 'Coherent a=' in name:
        analytic_tail = coherent_tail(N, float(name.split('=')[1]))
    elif 'Thermal nb=' in name:
        analytic_tail = thermal_tail(N, float(name.split('=')[1]))
    elif 'Squeezed r=' in name:
        analytic_tail = squeezed_vacuum_tail(N, float(name.split('=')[1]))
    tail_display = P[-1] if analytic_tail is None else analytic_tail
    if tail_display > 1e-4:
        print(f"  WARNING: {name} has non-negligible tail at N={N}")
    else:
        print(f"  OK: {name} (tail={tail_display:.2e})")

header = (
    f"{'State':<16} {'<n>':>6} {'Dn':>6} {'g2(0)':>8} {'Q':>7} "
    f"{'Wigner':>9} {'Statistics':<14}"
)
print(f"\n{'='*72}")
print("  MASTER STATISTICS TABLE: QUANTUM STATES OF LIGHT")
print(f"{'='*72}")
print(header)
print("-" * 72)

for name, state in master_states:
    mn = mean_photon_number(state, a)
    var = photon_variance(state, a)
    dn = np.sqrt(max(var, 0))

    if mn < 1e-10:
        g2_str = "  ---  "
        Q_str = "  --- "
        stat_type = "---"
    else:
        g2 = compute_g2_zero(state, a)
        Q = mandel_Q(state, a)
        g2_str = f"{g2:>7.3f}"
        Q_str = f"{Q:>6.2f}"
        if Q < -0.01:
            stat_type = "Sub-Poisson"
        elif abs(Q) < 0.01:
            stat_type = "Poissonian"
        else:
            stat_type = "Super-Poisson"

    W_sign = qutip.wigner(state, xvec_sign, xvec_sign)
    _, w_min, _ = wigner_sign_summary(W_sign)
    wigner_sign = "negative" if w_min < -1e-4 else "nonneg"

    print(f"{name:<16} {mn:>6.2f} {dn:>6.2f} {g2_str} {Q_str} "
          f"{wigner_sign:>9} {stat_type:<14}")

print(f"{'='*72}")
print("\n[PASS] Master statistics table computed (all values from code).\n")

print("Phase 5 Step 5.1 COMPLETE.")
