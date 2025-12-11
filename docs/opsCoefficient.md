# `opsCoefficient` Class Documentation

## Overview

The `opsCoefficient` class constructs all pairwise operator combinations needed for building the QITE linear-system coefficients. Given an ansatz

\[
\{\sigma_1, \sigma_2, \ldots, \sigma_N\},
\]

it produces, for every pair \(i \le j\),

\[
-\{\sigma_i, \sigma_j\} = -(\sigma_i \sigma_j + \sigma_j \sigma_i),
\]

stores the result as a `pauliString`, and precomputes its Pauli map.  
All results are stored in a 2D upper-triangular cell array `cellPS`.

This class is used by `qiteKernel` when assembling the coefficient matrix of the linear system.

---

## Constructor

```matlab
OC = opsCoefficient(PDH, ansatz)
```

### Arguments

| Argument | Type                            | Description                                              |
|----------|---------------------------------|----------------------------------------------------------|
| `PDH`    | `pauliDictionary` handle        | Shared Pauli dictionary                                  |
| `ansatz` | `1 x N` array of `pauliString`  | Ansatz operators relevant to the current Hamiltonian term |

The ansatz **must** be a row vector. A size mismatch causes an error.

Example:

```matlab
OC = opsCoefficient(PDH, ansatz);
```

---

## Properties

| Property       | Description                                                       |
|----------------|-------------------------------------------------------------------|
| `PDH`          | Handle to the shared `pauliDictionary`                             |
| `size_ansatz`  | Ansatz length                                                      |
| `cellPS`       | 2D cell array storing `- {ansatz(i), ansatz(j)}` for all `j >= i` |

`cellPS{i,j}` contains a single `pauliString` object.  
The lower triangle (`i > j`) is intentionally left empty.

---

## Internal Construction Logic

The constructor:

1. Checks that `ansatz` is a row vector.
2. Determines `size_ansatz`.
3. Allocates `cellPS`.
4. For every pair `i, j` with `j >= i`:

```matlab
tI = ansatz(i);
tJ = ansatz(j);

ps = -1 * anticommutator(tI, tJ);
ps = realPS(ps);
ps = ps.calcMap();

cellPS{i, j} = ps;
```

### Notes

- `anticommutator(a, b)` produces `a*b + b*a` as a `pauliString`.
- Only the **real** part is kept (`realPS(...)`) since imaginary parts correspond purely to phase and do not contribute to Hermitian operators.
- `calcMap()` precomputes the Pauli map and caches it inside `PDH`, allowing fast expectation-value calculations later.

---

## Relationship to Other Classes

### `pauliDictionary`

- Used to build maps for each resulting Pauli string.
- Ensures measurement caching and indexing works across the QITE pipeline.

### `pauliString`

- Used to represent all operator products.
- Provides:
  - `anticommutator`
  - `realPS`
  - `calcMap`

### `qiteKernel`

- Uses `opsCoefficient` to retrieve `- {σ_i, σ_j}` when assembling the linear system matrix for QITE.

---

## Example Usage

```matlab
% ansatz is 1 x N row vector of pauliString
OC = opsCoefficient(PDH, ansatz);

% Retrieve the operator -{sigma_2, sigma_5}
ps_25 = OC.cellPS{2, 5};

% Evaluate expectation
e = ps_25.eValLS(statevector);
```

The `opsCoefficient` object is therefore a compact precomputed operator pool used during QITE coefficient generation.
