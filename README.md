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

### Bell/CHSH Violation

![CHSH S-parameter exceeding the classical bound S = 2](figures/bell_chsh_violation.png)

The singlet-state CHSH calculation reaches S = 2 sqrt(2), exceeding the
local-realist bound while respecting no-signaling.

## Physics Background

The electromagnetic field in a cavity can be decomposed into normal modes, each
mode behaving as a quantum harmonic oscillator with Hamiltonian

```text
H = hbar omega (a_dag a + 1/2).
```

The creation and annihilation operators satisfy

```text
[a, a_dag] = 1
```

and generate the Fock-space basis used throughout the notebooks.

| State | Definition | g^2(0) | Key Property |
|-------|------------|--------|--------------|
| Fock \|n> | Number eigenstate | 1 - 1/n | Sub-Poissonian, non-classical |
| Coherent \|alpha> | D(alpha)\|0> | 1.0 | Poissonian photon statistics |
| Squeezed S(xi)\|0> | Squeezed vacuum | 3 + 1/<n> | Reduced uncertainty in one quadrature |
| Thermal rho_th | Mixed state | 2.0 | Super-Poissonian bunching |

<details>
<summary><strong>Extended Statistics Table</strong> (12 states, computed from code)</summary>

| State | <n> | Delta n | g^2(0) | Q | Wigner | Statistics |
|-------|----:|--------:|-------:|--:|--------|------------|
| Vacuum \|0> | 0.00 | 0.00 | - | - | nonnegative | - |
| Fock \|1> | 1.00 | 0.00 | 0.000 | -1.00 | negative | Sub-Poisson |
| Fock \|3> | 3.00 | 0.00 | 0.667 | -1.00 | negative | Sub-Poisson |
| Fock \|5> | 5.00 | 0.00 | 0.800 | -1.00 | negative | Sub-Poisson |
| Coherent alpha=1 | 1.00 | 1.00 | 1.000 | 0.00 | nonnegative | Poissonian |
| Coherent alpha=3 | 9.00 | 3.00 | 1.000 | 0.00 | nonnegative | Poissonian |
| Coherent alpha=5 | 25.00 | 5.00 | 1.000 | 0.00 | nonnegative | Poissonian |
| Squeezed r=0.5 | 0.27 | 0.83 | 6.683 | 1.54 | nonnegative | Super-Poisson |
| Squeezed r=1.0 | 1.38 | 2.56 | 3.724 | 3.76 | nonnegative | Super-Poisson |
| Thermal n_bar=1 | 1.00 | 1.41 | 2.000 | 1.00 | nonnegative | Super-Poisson |
| Thermal n_bar=3 | 3.00 | 3.46 | 2.000 | 3.00 | nonnegative | Super-Poisson |
| Thermal n_bar=5 | 5.00 | 5.48 | 2.000 | 5.00 | nonnegative | Super-Poisson |

</details>

## Notebooks

| # | Notebook | Topic | Main Checks |
|---|----------|-------|-------------|
| 1 | `01_fock_states.ipynb` | Fock states | Photon-number states, g^2(0), Wigner negativity |
| 2 | `02_coherent_states.ipynb` | Coherent states | Poisson statistics and displaced Gaussian Wigner functions |
| 3 | `03_squeezed_states.ipynb` | Squeezed states | Quadrature squeezing, uncertainty ellipses, even-photon structure |
| 4 | `04_photon_statistics.ipynb` | Photon statistics | g^2(0), Mandel Q, and comparison tables |
| 5 | `05_beam_splitter_hom.ipynb` | Beam splitters and HOM | Two-photon interference and P(1,1) = 0 |
| 6 | `06_entanglement_bell.ipynb` | Entanglement and Bell tests | CHSH violation, correlations, and concurrence |

## Reproduce

### Prerequisites

- Python 3.14 tested; see `environment_report.txt` for exact versions.
- pip.

### Installation

```bash
git clone https://github.com/abbas-ahmad-cowlar/quantum-states-of-light.git
cd quantum-states-of-light
python -m venv venv
venv\Scripts\activate
pip install -r requirements.txt
```

### Run notebooks

```bash
jupyter notebook notebooks/
```

### Run tests

```bash
python -m pytest tests/ -q
python verify_notebooks.py
```

## Project Structure

```text
quantum-states-of-light/
├── README.md
├── requirements.txt
├── requirements-lock.txt
├── environment_report.txt
├── LICENSE
├── verify_notebooks.py
├── src/
│   └── qo_utils.py
├── notebooks/
│   ├── 01_fock_states.ipynb
│   ├── 02_coherent_states.ipynb
│   ├── 03_squeezed_states.ipynb
│   ├── 04_photon_statistics.ipynb
│   ├── 05_beam_splitter_hom.ipynb
│   └── 06_entanglement_bell.ipynb
├── figures/
│   └── generated plots in PNG and SVG formats
└── tests/
    ├── test_qo_utils.py
    ├── test_core_functions.py
    ├── qutip_orientation.py
    ├── smoke_test.py
    └── truncation_check.py
```

## Classifying Quantum Light

```text
                    Measure g^2(0)
                          |
            +-------------+-------------+
            |             |             |
        g^2(0) < 1    g^2(0) = 1    g^2(0) > 1
            |             |             |
      Antibunching    Poissonian      Bunching
     Sub-Poissonian      Q = 0     Super-Poissonian
            |                           |
       Fock states                 Thermal light,
                                  squeezed vacuum
```

Use g^2(0) together with photon-number distributions, Wigner functions, mean
photon number, and quadrature statistics. A single scalar diagnostic does not
uniquely determine the state.

## Author

Syed Abbas Ahmad

## References

- Gerry and Knight, *Introductory Quantum Optics* (Cambridge, 2005)
- Fox, *Quantum Optics: An Introduction* (Oxford, 2006)
- [QuTiP Documentation](https://qutip.org)
