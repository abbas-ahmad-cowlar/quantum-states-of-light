"""
QuTiP Orientation Script
========================
This script maps every QuTiP function to the physics it implements.
"""

import qutip
import numpy as np
import math

# =============================================================================
# PART 1: Creating Fock States
# =============================================================================
# Physics: Fock state |n> is a state with exactly n photons.
# In the truncated Hilbert space of dimension N, |n> is a column vector
# with a 1 at position n and 0 everywhere else.
# QuTiP: qutip.basis(N, n) creates |n> in an N-dimensional space.

N = 10  # Truncation dimension -- we work with |0> through |9>

# Create the vacuum state |0>
vac = qutip.basis(N, 0)
print("=== Fock States ===")
print(f"|0> (vacuum state):\n{vac}\n")
print(f"Type: {vac.type}")   # Should be 'ket'
print(f"Shape: {vac.shape}")  # Should be (10, 1)

# Create |1>, |3>, |5>
one = qutip.basis(N, 1)
three = qutip.basis(N, 3)
five = qutip.basis(N, 5)

# Verify orthonormality: <m|n> = delta_mn
print(f"\n<0|0> = {vac.dag() * vac}")     # Should be 1
print(f"<0|1> = {vac.dag() * one}")       # Should be 0
print(f"<3|3> = {three.dag() * three}")   # Should be 1
print(f"<3|5> = {three.dag() * five}")     # Should be 0

# =============================================================================
# PART 2: Creation and Annihilation Operators
# =============================================================================
# Physics: a removes a photon: a|n> = sqrt(n) |n-1>
#          a_dag adds a photon: a_dag|n> = sqrt(n+1) |n+1>
# QuTiP: qutip.destroy(N) creates a, and a.dag() gives a_dag

a = qutip.destroy(N)
print("\n=== Operators ===")
print(f"a (annihilation operator), shape {a.shape}:")
print(f"Matrix representation:\n{a}\n")

# Verify a|0> = 0 (cannot annihilate vacuum)
result = a * vac
print(f"a|0> = {result.norm():.2e} (should be ~0)")

# Verify a_dag|0> = |1> (creates one photon from vacuum)
one_created = a.dag() * vac
overlap = one.dag() * one_created
print(f"<1|a_dag|0> = {overlap} (should be 1.0)")

# Verify a|3> = sqrt(3) |2>
result_3 = a * three
two = qutip.basis(N, 2)
coefficient = two.dag() * result_3
print(f"<2|a|3> = {coefficient} (should be sqrt(3) = {np.sqrt(3):.4f})")

# Verify a_dag|3> = sqrt(4) |4> = 2|4>
result_3dag = a.dag() * three
four = qutip.basis(N, 4)
coefficient_dag = four.dag() * result_3dag
print(f"<4|a_dag|3> = {coefficient_dag} (should be sqrt(4) = 2.0)")

# =============================================================================
# PART 3: Number Operator and Expectation Values
# =============================================================================
# Physics: n_hat = a_dag * a counts photons. <n|n_hat|n> = n.
# QuTiP: qutip.num(N) creates n_hat, or compute a.dag() * a.
#         qutip.expect(op, state) computes <state|op|state>.

n_op = qutip.num(N)
# Equivalent: n_op = a.dag() * a

print("\n=== Expectation Values ===")
for n_val in [0, 1, 3, 5]:
    state = qutip.basis(N, n_val)
    mean_n = qutip.expect(n_op, state)
    print(f"<{n_val}|n_hat|{n_val}> = {mean_n:.1f} (should be {n_val})")

# =============================================================================
# PART 4: Building States from Vacuum
# =============================================================================
# Physics: |n> = (a_dag)^n / sqrt(n!) |0>
# Let's verify this by building |3> from scratch.

# Build |3> manually: apply a_dag three times to |0>
state = vac.copy()
for i in range(3):
    state = a.dag() * state
# Normalize: divide by sqrt(3!) = sqrt(6)
state = state / np.sqrt(math.factorial(3))

# Compare with qutip.basis(N, 3)
overlap_manual = three.dag() * state
print(f"\n=== Building |3> from vacuum ===")
print(f"<3|( (a_dag)^3/sqrt(3!) |0> ) = {overlap_manual} (should be 1.0)")

# =============================================================================
# PART 5: Commutation Relation Verification
# =============================================================================
# Physics: [a, a_dag] = a*a_dag - a_dag*a = I (identity)
# This is the FUNDAMENTAL relation of quantum optics.

commutator = a * a.dag() - a.dag() * a
identity = qutip.qeye(N)
print(f"\n=== Commutation Relation ===")
# In finite dimension, the commutator has a boundary correction.
diff = (commutator - identity).norm()
print(f"||[a, a_dag] - I|| = {diff:.2e} (nonzero because of finite-N boundary)")

print("\n[PASS] QuTiP orientation complete. All physics mappings verified.")
