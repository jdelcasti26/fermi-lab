# `hamiltonian` Class Documentation

## Overview

The `hamiltonian` class constructs a partitioned Hamiltonian and its corresponding QITE ansatz from model-specific master data files.

It:

- Loads a Hamiltonian model (fermionic or spin) defined in master files.
- Groups Hamiltonian terms into partitions chosen by the user.
- Builds one local Hamiltonian per partition element.
- Associates one ansatz per local Hamiltonian, based on domain rules.
- Provides exact diagonalization utilities for reference.

All Hamiltonian terms and ansatz operators are represented as `pauliString` objects and rely on a shared `pauliDictionary` instance.

---

## Data Sources and Configuration

The constructor reads:

- `params/cfg.json` (configuration file)
- master data files in the `master/` directory

The JSON config includes per-model information:

- `model` (string)
- `type` (`"fermions"` or `"spins"`)
- `hamiltonian` (file defining HH)
- `ansatz` (file defining AZ)
- `domain` (optional)
- `parameter` (optional)

Master files define:

- `HH` — base Hamiltonian terms (`pauliString` array)
- `AZ` — ansatz operators (`pauliString` array)
- `PP` — parameter scalars (optional)
- `DD` — domain definitions (optional)

The number of Hamiltonian terms must match across all definitions.

---

## Constructor

```matlab
H = hamiltonian(PDH, model, p, varargin)
```

### Arguments

| Argument    | Type                     | Description                                       |
|-------------|--------------------------|---------------------------------------------------|
| `PDH`       | `pauliDictionary` handle | Shared Pauli dictionary                           |
| `model`     | string                   | Model name (must exist in cfg.json)               |
| `p`         | scalar                   | Model parameter; use `1` if not needed            |
| `varargin`  | vectors                  | Partition term indices (one vector per partition) |

Example:

```matlab
partition = { [1:3], [4:15] };
H = hamiltonian(PDH, "H2", 1.0, partition{:});
```

This creates 2 grouped Hamiltonian terms:  
`H.H(1) = HH(1) + HH(2) + HH(3)`  
`H.H(2) = HH(4) + ... + HH(15)`

---

## Properties

| Property     | Description                                                                    |
|--------------|--------------------------------------------------------------------------------|
| `EPSEIG`     | Degeneracy tolerance for eigenvalues                                           |
| `PDH`        | Handle to the shared `pauliDictionary`                                         |
| `type`       | `"fermions"` or `"spins"`                                                      |
| `model`      | Model name                                                                     |
| `H`          | Array of partitioned Hamiltonian terms (each a `pauliString`)                  |
| `M`          | Number of partitions                                                           |
| `HT`         | Total Hamiltonian = sum of all `H(m)`                                          |
| `cellAZ`     | Cell array of ansatz lists, one per partition                                  |

---

## Hamiltonian Construction

### 1. Model lookup and validation

- Reads `cfg.json`
- Finds entry where `model` matches
- Ensures:
  - valid `type`
  - Hamiltonian file defined
  - ansatz file defined

### 2. Loading master definitions

```matlab
run(masterPath + stuConfig(idr).hamiltonian);   % defines HH
run(masterPath + stuConfig(idr).ansatz);        % defines AZ
```

Hermiticity checks:

- Each `HH(t)` must satisfy: `HH(t) - dagger(HH(t)) == 0`
- Each `AZ(t)` must satisfy: `AZ(t) + dagger(AZ(t)) == 0` (anti-Hermitian)

### 3. Parameter scaling (optional)

If `parameter` file is present:

```matlab
run(masterPath + stuConfig(idr).parameter);  % defines PP
```

Each Hamiltonian term may be scaled:

```
HH(i) = (p * PP(i)) * HH(i)   if PP(i) ≠ 0
```

### 4. Grouping into partition Hamiltonians

For each partition index list:

```matlab
arrH = HH(varargin{m});     % base terms
H(m) = arrH(1);
for s = 2:size(arrH,2)
    H(m) = H(m) + arrH(s);
end
```

The total Hamiltonian is:

```matlab
HT = H(1) + H(2) + ... + H(M);
```

---

## Ansatz Construction

The ansatz for each local Hamiltonian is stored in `cellAZ{m}`.

Two main cases:

---

### Case 1 — No domains defined

The global ansatz `AZ` is replicated for each partition:

```matlab
for m = 1:M
    cellAZ{m} = AZ;
end
```

---

### Case 2 — Domains defined

```matlab
run(masterPath + stuConfig(idr).domain);   % defines DD
```

The size of `DD` must equal the number of base Hamiltonian terms.

### Fermionic models (`type == "fermions"`)

- `DD` is a cell array.
- `DD(i)` contains ansatz indices associated with Hamiltonian term `HH(i)`.

For each partition:

1. Collect all domains for its terms
2. Merge and deduplicate the ansatz indices
3. Assign:

```matlab
cellAZ{m} = AZ(idd);
```

### Spin models (`type == "spins"`)

- `DD` is an array of one-term `pauliString` masks (e.g., `"IIXXII"`).
- `AZ` must contain only one-term `pauliString`s.

For each partition, a mask matches an ansatz string via:

```
template 'I' positions must also be 'I' in the target string
```

Using:

```matlab
flag = matchesI(template, target);
```

All matching ansatz terms are included.

---

## Exact Diagonalization Utilities

### All eigenvalues

```matlab
vals = H.exact_eigenvalues();
```

Computes eigenvalues of `HT.matrixRep()`.

### Ground-state energy

```matlab
e0 = H.exact_eigenvalue();
```

Minimum eigenvalue.

### Ground-state eigenvectors (possibly degenerate)

```matlab
V = H.exact_eigenvectors();
```

- Sorts eigenvalues
- Selects all eigenvectors with eigenvalue within tolerance `EPSEIG` of the minimum

---

## Relationship to Other Classes

- **`pauliDictionary`**  
  Provides Pauli matrices, maps, and caching.

- **`pauliString`**  
  Used to store all Hamiltonian terms and ansatz operators.

- **`qiteKernel`**  
  Consumes `H` and `cellAZ` to build operator pools and assemble QITE linear systems.

- **`qite`**  
  Uses the `hamiltonian` output to run QITE iterations and benchmarking.

---

## Example

```matlab
PDH = pauliDictionary(4);
partition = { [1:3], [4:15] };

H = hamiltonian(PDH, "H2", 1.0, partition{:});

local1  = H.H(1);
ansatz1 = H.cellAZ{1};

E0 = H.exact_eigenvalue();
V0 = H.exact_eigenvectors();
```

The `hamiltonian` class provides all model-specific logic needed to construct partitioned Hamiltonians and ansätze for QITE simulations.
