# `pauliString` Class Documentation

## Overview

The `pauliString` class represents a finite linear combination of tensor-product Pauli operators acting on an \(n\)-qubit Hilbert space. A typical object encodes expressions such as

\(0.5\, XXYYZZ - 1.0\, XIIIIX + \tfrac{1}{3}\, ZZZZZZ\)

where each term is a Pauli string with a complex coefficient.

Internally:

- Wavefunctions are complex MATLAB vectors of size \(2^n\).
- Pauli operators are stored in symbolic form as strings, e.g. `"XXYYZZ"`.
- A shared `pauliDictionary` instance provides:
  - Qubit count and dimensions
  - Explicit one-qubit Pauli matrices
  - Precomputed action maps for each Pauli string
  - Cached measurement results and measurement counters

The `pauliString` class focuses on:

- Pauli algebra (sum, product, commutator, etc.)
- Efficient expectation-value evaluation using cached maps
- Optional full matrix representations (for debugging)

---

## Constructor

```matlab
ps = pauliString(PH, coeffs, strings)
```

**Arguments:**

| Argument  | Type                       | Description                                                                 |
|----------|----------------------------|-----------------------------------------------------------------------------|
| `PH`     | `pauliDictionary` handle   | Provides Pauli matrices, maps, measurement cache, and qubit count.         |
| `coeffs` | `len x 1` complex vector   | Coefficients of each Pauli term.                                           |
| `strings`| `len x 1` cell array       | Pauli strings such as `"XXI"`, `"IZZ"`, etc.                               |

**Behavior:**

- Truncates strings to match `PH.nqbs`.
- Validates characters (`I`, `X`, `Y`, `Z`).
- Sorts terms and merges duplicates.
- Removes zero-coefficient terms.

Example:

```matlab
PH = pauliDictionary(4);
coeffs  = [0.5; -1.0; 1/3];
strings = {"XXYY"; "XIII"; "ZZZZ"}';
ps = pauliString(PH, coeffs, strings);
```

---

## Properties

| Property        | Type                         | Description                                                        |
|----------------|------------------------------|--------------------------------------------------------------------|
| `PH`           | `pauliDictionary` handle      | Provides maps, Pauli matrices, cache, qubit count.                 |
| `updateNeeded` | logical                       | Whether maps must be generated via `calcMap`.                      |
| `len`          | integer                       | Number of terms in the Pauli sum.                                  |
| `coef`         | `len x 1` complex vector      | Coefficients.                                                       |
| `pStr`         | `len x 1` cell array          | Pauli strings (length = `PH.nqbs`).                                |

---

## Internal Utilities

### Sorting and merging terms

```matlab
o = o.intnSort()
```

- Sorts `pStr` lexicographically.
- Sums coefficients of duplicate strings.
- Removes zero-coefficient entries.

### Pauli string multiplication (with phases)

```matlab
[ph, chs] = o.intnMult(chs1, chs2)
```

- Computes the Pauli product qubit-by-qubit.
- Returns a phase `ph ∈ {1, -1, +i, -i}` and resulting string `chs`.

### Commutation test

```matlab
flag = o.intnComm(chs1, chs2)
```

Two strings commute iff the number of conflicting non-identity operators is even.

### Commutation structure

```matlab
[S, SMat] = o.commStructure()
```

- `S`: list of term pairs and commute flag.
- `SMat`: symmetric matrix with 1 (commuting) or 0 (anticommuting).

---

## Operator Overloading

### Addition and subtraction

```matlab
s = ps1 + ps2;
m = ps1 - ps2;
```

- Concatenates terms and re-simplifies via the constructor.

### Scalar multiplication

```matlab
m = alpha * ps;
```

### Product of two Pauli strings

```matlab
p = ps1 * ps2;
```

Implements:

\((\sum_i a_i A_i)(\sum_j b_j B_j) = \sum_{i,j} a_i b_j (A_i B_j)\)

### Commutator and anticommutator

```matlab
C = ps1.commutator(ps2);      % ps1*ps2 - ps2*ps1
A = ps1.anticommutator(ps2);  % ps1*ps2 + ps2*ps1
```

### Hermitian conjugate

```matlab
D = ps.dagger();
```

### Squaring

```matlab
sq = ps ^ 2;
```

Only `^2` is implemented.

---

## Expectation Value Evaluation

Expectation values use precomputed maps in `PH.DIC`. Each map describes how a Pauli string acts on a state vector without full matrix multiplication.

State `sv` may be:

- A single state vector (`2^n x 1`)
- A matrix of weighted states (each column has a probability weight)

### `eValLS` — linear system assembly

```matlab
e = ps.eValLS(sv);
```

- Reuses cached results when possible (`PH.DICOM`).
- Otherwise computes new measurements and caches them.
- Increments `PH.measuresLS`.

### `eValCO` — operator-coefficient containers

```matlab
e = ps.eValCO(sv);
```

Same as `eValLS`, but increments `PH.measuresCO`.

### `eValAZ` — analysis (no cache reuse)

```matlab
e = ps.eValAZ(sv);
```

- Never reuses cache.
- Always increments `PH.measuresAZ`.

---

## Map Calculation

```matlab
ps = ps.calcMap();
```

For each Pauli string:

1. If not in `PH.DIC`:
   - Allocate arrays `Idx`, `TR`, `TI`.
   - For each basis vector:
     - Compute tensor indices.
     - Determine output column and coefficient.
     - Store in the map arrays.
   - Add entry to `PH.DIC`.
   - Add measurement cache entry to `PH.DICOM`.

2. Set `updateNeeded = false`.

This avoids building full dense matrices.

---

## Matrix Representations

### Single-string matrix

```matlab
M = ps.strToMatrix(inputStr);
```

Builds the Kronecker product of one-qubit matrices.

### Full operator matrix

```matlab
M = ps.matrixRep();
```

Computes:

\(\displaystyle M = \sum_i coef(i) * strToMatrix(pStr{i})\)

---

## Display and String Utilities

### Console display

```matlab
disp(ps)
```

Shows each term `coef * pStr`.

### String representation

```matlab
s = ps.string();
```

Returns the printable representation.

### Empty/zero check

```matlab
flag = ps.isemptyorzero();
```

---

## Relationship to Other Classes

- **`pauliDictionary`**
  - Provides Pauli matrices, symbolic maps, measurement caches.
  - `pauliString` uses dictionary maps for all fast operations.

- **`qiteKernel`**
  - Builds operator pools using `pauliString`.
  - Uses `eValLS` / `eValCO` for assembling QITE linear systems.

- **`qite`**
  - Controls state updates and triggers measurement-cache resets.

Together, `pauliString` and `pauliDictionary` form the symbolic and caching backbone enabling efficient classical QITE emulation without full matrix operations.