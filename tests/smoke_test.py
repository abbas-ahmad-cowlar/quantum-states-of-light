"""
Smoke Test for QuTiP Installation
===================================
Verifies that QuTiP is properly installed and can perform basic quantum
optics operations: creating states, operators, and computing Wigner functions.
"""

import qutip
import numpy as np
import matplotlib
matplotlib.use('Agg')  # Non-interactive backend for testing
import matplotlib.pyplot as plt

# Test 1: Create vacuum state
vac = qutip.basis(10, 0)
assert vac.shape == (10, 1), "Vacuum state has wrong shape"
print("[PASS] Test 1: Vacuum state created successfully")

# Test 2: Create annihilation operator
a = qutip.destroy(10)
assert a.shape == (10, 10), "Annihilation operator has wrong shape"
print("[PASS] Test 2: Annihilation operator created successfully")

# Test 3: Verify a|0> = 0
result = a * vac
assert result.norm() < 1e-10, "a|0> should be zero"
print("[PASS] Test 3: a|0> = 0 (cannot annihilate vacuum)")

# Test 4: Verify a_dag|0> = |1>
one_photon = a.dag() * vac
expected = qutip.basis(10, 1)
assert (one_photon - expected).norm() < 1e-10, "a_dag|0> should be |1>"
print("[PASS] Test 4: a_dag|0> = |1> (creates one photon)")

# Test 5: Wigner function computes without error
xvec = np.linspace(-5, 5, 100)
W = qutip.wigner(vac, xvec, xvec)
assert W.shape == (100, 100), "Wigner function has wrong shape"
print("[PASS] Test 5: Wigner function computed successfully")

print("\n[PASS] All smoke tests passed. QuTiP is working correctly.")
