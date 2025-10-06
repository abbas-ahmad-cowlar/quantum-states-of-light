# Quantum States of Light

[![Python](https://img.shields.io/badge/Python-3.14%20tested-blue.svg)](https://python.org)
[![QuTiP](https://img.shields.io/badge/QuTiP-5.x-orange.svg)](https://qutip.org)
[![License](https://img.shields.io/badge/License-MIT-green.svg)](LICENSE)
[![Jupyter](https://img.shields.io/badge/Jupyter-Notebooks-F37626.svg)](https://jupyter.org)

Computational quantum optics examples using Python and
[QuTiP](https://qutip.org). The repository constructs common quantum states of
light, computes photon statistics and Wigner functions, and simulates standard
two-mode experiments such as Hong-Ou-Mandel interference and CHSH Bell
violation.

## Scope

The project covers:

- Quantum state construction: Fock, coherent, squeezed vacuum, squeezed
  coherent, and thermal states.
- Photon statistics: photon-number distributions, zero-delay coherence
  g^2(0), and Mandel Q.
- Phase-space analysis: Wigner functions, normalization checks, and Wigner
  negativity.
- Two-mode optics: beam splitter transformations and Hong-Ou-Mandel
  interference.
- Entanglement diagnostics: Bell states, concurrence, and CHSH violation.
- Reproducibility: fixed dependency records, validation tests, and notebook
  verification.

## Key Results

### Wigner Function Gallery

![Wigner functions of six fundamental quantum states of light](figures/wigner_gallery_6panel.png)

Six phase-space portraits on a shared X,P convention, including Wigner
negativity for Fock states and Gaussian positivity for coherent,
squeezed-vacuum, and squeezed-coherent examples.

### Photon Statistics

![g^2(0) bar chart showing zero-delay photon-statistics regimes](figures/g2_bar_chart_comparison.png)

Zero-delay photon-correlation regimes: antibunching, Poissonian light, and
bunching. The value of g^2(0) is a useful diagnostic, not a complete state
classifier.

### Hong-Ou-Mandel Interference

![HOM dip simulation showing P(1,1) = 0 for identical photons](figures/hom_effect_and_dip.png)

Identical photons at a 50:50 beam splitter show destructive two-photon
interference and a zero-coincidence dip at zero delay.

