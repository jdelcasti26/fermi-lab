# `mtqite` — Multi-term QITE driver (MT-QITE)

`mtqite` is a high-level driver built on top of `qite_kernel` to run **multi-term / partitioned QITE** across:
- a scan of candidate imaginary time steps `Dt` (the `timeSeries` grid), and
- multiple **Trotter steps** (outer loop).

At each Trotter step, the class:
1. **measures and solves** the QITE linear systems (once per partition term, for the smallest `Dt`),
2. **reuses** those measurements to generate QITE update circuits for all candidate `Dt`,
3. **evaluates** all combinations of per-term circuits (according to a chosen permutation) and
4. **selects** the circuit combination that minimizes the energy, then evolves the state.

---

## What problem does it solve?

Given a Hamiltonian partitioned into `M` terms (as defined in `qite_kernel`),
MT-QITE constructs term-wise unitary approximations to imaginary-time evolution and then
**chooses**, at each Trotter step, the best combination of term-wise updates by **energy minimization**.

This implements a practical strategy:
- QITE measurements are expensive → measure once (for a reference `Dt`) and solve for multiple `Dt`.
- Multi-term ordering matters → search over per-term `Dt` combinations and pick the lowest-energy circuit.

---

## Key dependencies

`mtqite` assumes `qite_kernel` provides:
- `set_statevector(stateVector)`
- `measure_state(mterm)`
- `measureCount(mterm)`
- `solve_linear_system(mterm, Dt)` → returns `[thetas, qiteOps]`
- access to:
  - `QK.M` (number of partition terms)
  - `QK.PDH.twoN` (Hilbert space dimension, i.e. `2^n`)
  - `QK.PDH.flushMeasurements()`, `QK.PDH.measuresCO` (measurement counter)
  - `QK.HT.eValCO(state)` (energy estimator)
  - `QK.exact_eigenvalue()`, `QK.exact_eigenvectors()` (for accuracy/fidelity diagnostics)

`qiteOps` entries are expected to support `.matrixRep()` and have `.len` and `.pStr` fields (Pauli string container).

---

## Class overview

### Properties

#### Runtime configuration
- `QK`  
  Handle to the `qite_kernel` instance.
- `initialState`  
  Initial normalized statevector (size `2^n × 1`).
- `timeSeries`  
  Row vector of candidate **positive** time steps `Dt` (sorted and made unique).
- `noTIMES`  
  Number of time steps in `timeSeries`.
- `noTS`  
  Number of Trotter steps.

#### Simulation results
- `cellQiteOps(noTIMES, noTS, M)`  
  Cell array storing, for each `(Dt index, trotterStep, mterm)`, the list/array of Pauli operators defining the QITE update.
- `cellThetas(noTIMES, noTS, M)`  
  Same shape as above; stores the corresponding operator coefficients `theta`.
- `measureCountLS(noTIMES, noTS, M)`  
  Measurement cost recorded for building the **linear system** per term.  
  *Note:* In current code, only the first `Dt` (smallest in `timeSeries`) triggers measurements, but the count is stored at all `tn` indices only where measured.
- `measureCountEE(noTIMES, noTS)`  
  Placeholder for energy-estimation measurement counts (allocated; energy measurement cost is tracked per circuit grid in `estimate_energy_step`).

---

## Constructor

### `mtqite(QK, initialState, timeSeries, noTS)`

Validates and initializes the driver:
- checks `initialState` has dimension `QK.PDH.twoN × 1`
- normalizes `initialState`
- checks `timeSeries` is a row vector with positive entries
- sorts and uniques `timeSeries`
- allocates internal storage cells and counters

---

## Main workflow

### `estimate_energy(perm, symmetryType)`

Runs MT-QITE for `noTS` Trotter steps starting from `initialState`.

**Inputs**
- `perm`  
  Permutation vector of length `M` describing the order of term application in the circuit.
  Example for `M=3`: `[3 1 2]` means apply term 3, then 1, then 2.
- `symmetryType`  
  Supported values:
  - `"none"`
  - `"3rd_eq_1st"`: assumes term 3 is a mirror-symmetric copy of term 1.

**Algorithm (per Trotter step `ts`)**
1. `scan_step(...)`: for each independent term and each `Dt`, store `qiteOps` and `thetas`.
2. `estimate_energy_step(...)`: build all combinations of term-wise circuits across `Dt` grid, compute energy for each, select minimum.
3. Evolve state with the best circuit and proceed to the next step.

**Outputs**
- `stuSummary` (struct, scalar)  
  Contains arrays over Trotter steps:
  - `energy` (min energy per step)
  - `accuracy` (relative error vs exact eigenvalue)
  - `fidelity` (overlap with exact eigenvector set as implemented)
  - `minPosition` (chosen `tn` indices per term)
  - `opCount`, `PSCount` (cumulative operator/string depth)
  - `costLS`, `costEE` (cumulative measurement cost)
- `endState`  
  Final evolved statevector after `noTS` steps.

> The method prints diagnostics per step: energy, accuracy, depth, fidelity, and measurement totals.

---

## Precomputation and symmetry

### `scan_step(symmetryType, stateVector, trotterStep)`

Precomputes and stores QITE circuits for a fixed `stateVector` and Trotter step:
- sets state in `qite_kernel`
- flushes measurement counters
- for each term:
  - measures linear system **only once** (for smallest `Dt`)
  - solves linear system for each `Dt` and stores resulting `(qiteOps, thetas)`

**Symmetry**
- `"none"`: loops over all terms `mterm = 1..M`
- `"3rd_eq_1st"`: treats terms 1 and 2 as independent, and fills term 3 by mirroring term 1 via `mirror(...)`

### `mirror(arrPS)`

Returns a mirror-symmetric version of an array of Pauli strings by reversing each `.pStr{s}` string.
Used to impose a reflection symmetry where term 3 is derived from term 1.

---

## Energy-based selection

### `estimate_energy_step(perm, symmetryType, trotterStep, stateVector)`

Enumerates all combinations of time-step indices `(tn1, tn2, ..., tnM)` and selects the combination
that minimizes the energy.

**Grid construction**
- Builds a Cartesian product over `tn ∈ {1..noTIMES}` for each term position.
- Optional symmetry filter for `"3rd_eq_1st"`:
  - currently implemented only for `perm == [1 2 3]`, enforcing `tn1 == tn3`.

**Circuit evaluation**
For each grid point:
1. Build the circuit as a product of term circuits:
   \[
     U = U_{\text{term}(p)}(Dt_{p}) \cdots U_{\text{term}(1)}(Dt_{1})
   \]
2. Evolve: `evolvedState = U * stateVector`
3. Compute energy: `QK.HT.eValCO(evolvedState)`
4. Track measurement cost from `PDH.measuresCO` deltas.

**Selection and outputs**
- Picks the grid row with minimum energy.
- Rebuilds the chosen circuit to compute:
  - `endState`
  - `fidelity` w.r.t. exact eigenvectors (as implemented)
- Returns `stuSummary` for this step and prints diagnostics.

---

## Circuit construction details

### `calculate_circuit(timeDatapoint, trotterStep, mterm)`

Builds the unitary for a single term `mterm` at a given time index `timeDatapoint`.

If `qiteOps = {A_i}` and `thetas = {θ_i}`, then:
\[
U_{\text{term}}(Dt) = \prod_i \exp\left(-Dt \; \theta_i \; A_i\right)
\]
where `A_i` is obtained via `matrixRep()`.

Returns:
- `Circuit` (matrix)
- `noOps` (# exponentials/operators in product)
- `noPSs` (sum of Pauli-string lengths across operators)

---

## Utility

### `cartesianProduct(sets)`

Computes a Cartesian product grid from a cell-array of sets, using `ndgrid`.
Used to enumerate all `tn` combinations across terms.

---

## Example usage

```matlab
% Assume qite_kernel has been constructed elsewhere as QK
% and initial state psi0 is normalized (or will be normalized by mtqite)

timeSeries = [0.01 0.02 0.05 0.1];
noTS = 20;

MT = mtqite(QK, psi0, timeSeries, noTS);

perm = [3 1 2];              % term application order
symmetryType = "none";       % or "3rd_eq_1st"

[MT, summary, psiEnd] = MT.estimate_energy(perm, symmetryType);

disp(summary.energy);
disp(summary.minPosition);
