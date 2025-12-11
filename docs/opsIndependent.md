# `opsIndependent` Class Documentation

## Overview

The `opsIndependent` class constructs a **vector of operator combinations** between a given operator and each element of an ansatz. Given:

- an ansatz \(\{\sigma_1, \sigma_2, \ldots, \sigma_N\}\), and  
- a single operator \(O\) (typically the Hamiltonian \(H\) or its square \(H^2\)),

this class computes, for each ansatz element \(i\),

\[
[\;O, \sigma_i\;] = O \sigma_i - \sigma_i O,
\]

and stores the **real part** of the resulting `pauliString` in a cell array.  
Each element is also preprocessed via `.calcMap()` for fast expectation-value evaluation.

This operator vector is used in the QITE kernel when building terms involving commutators between the Hamiltonian (or related operators) and the ansatz generators.

---

## Constructor

```matlab
OI = opsIndependent(PDH, ansatz, op)
```

### Arguments

| Argument | Type                            | Description                                                |
|----------|---------------------------------|------------------------------------------------------------|
| `PDH`    | `pauliDictionary` handle        | Shared Pauli dictionary                                    |
| `ansatz` | `1 x N` array of `pauliString`  | Ansatz operators (row vector)                             |
| `op`     | `pauliString`                   | Operator to commute with each ansatz term (e.g. `H` or `H^2`) |

The `ansatz` **must** be a row array. If not, the constructor throws an error.

Example:

```matlab
% ansatz is a 1 x N array of pauliString
% op is a single pauliString (e.g., the Hamiltonian term or H^2)
OI = opsIndependent(PDH, ansatz, op);
```

---

## Properties

```matlab
classdef opsIndependent
    properties
        PDH           % Handle to Pauli dictionary
        size_ansatz   % Ansatz length (number of operators)
        cellPS        % Cell array of resulting Pauli strings
    end
end
```

### Description

| Property      | Description                                                        |
|---------------|--------------------------------------------------------------------|
| `PDH`         | Handle to shared `pauliDictionary`                                 |
| `size_ansatz` | Number of ansatz operators (`N`)                                   |
| `cellPS`      | `1 x N` cell array; `cellPS{i}` stores the i-th commutator result |

---

## Internal Construction Logic

The constructor performs the following steps:

1. **Shape check for ansatz**

   ```matlab
   [m, o.size_ansatz] = size(ansatz);
   if (m ~= 1)
       error("Expected row array for relevant pauliStrings (operators)");
   end
   ```

   The ansatz must be a row vector of `pauliString` objects.

2. **Allocation and computation**

   For each ansatz element \( \sigma_i \):

   ```matlab
   for i = 1:o.size_ansatz
       tI = ansatz(i);

       % Compute commutator [op, tI] = op*tI - tI*op
       ps = commutator(op, tI);

       % Keep only real contributions (imaginary parts do not require measurement)
       ps = realPS(ps);

       % Precompute map if the result is non-empty
       if ps.len > 0
           ps = ps.calcMap();
       end

       o.cellPS{i} = ps;
   end
   ```

### Notes

- `commutator(op, tI)` is a method of `pauliString` that returns `op * tI - tI * op` as a new `pauliString`.
- `realPS(...)` discards purely imaginary contributions (only real parts need to be measured as observables).
- `calcMap()` builds and caches the Pauli map in `PDH` so that later expectation values can be computed efficiently.
- If the resulting `pauliString` is empty (`len == 0`), no map is generated and the term is effectively skipped for measurement.

---

## Relationship to Other Classes

### `pauliDictionary`

- Used indirectly via `calcMap()`, to:
  - Register new Pauli maps (`DIC`)
  - Initialize measurement entries (`DICOM`)

### `pauliString`

- Provides:
  - `commutator(op, tI)`
  - `realPS(...)`
  - `calcMap()`
- Represents both the ansatz operators and the resulting commutators.

### `qiteKernel`

- Uses `opsIndependent` to access the commutator vector \([O, \sigma_i]\).
- Typically, `op` is a Hamiltonian term or a related operator (e.g. \(H^2\)),
  and these commutators enter the right-hand side or coefficients of the QITE linear system.

---

## Example Usage

```matlab
% Given:
%   PDH     : pauliDictionary handle
%   ansatz  : 1 x N array of pauliString
%   H_local : pauliString representing a local Hamiltonian term

OI = opsIndependent(PDH, ansatz, H_local);

% Access the commutator [H_local, sigma_k]
ps_k = OI.cellPS{k};

% Evaluate an expectation value with a state vector sv
% (for instance, inside a QITE kernel routine)
e = ps_k.eValLS(sv);     % or eValCO/eValAZ as appropriate
```

The `opsIndependent` class is thus a compact container for precomputed commutators between a given operator and each ansatz generator, ready to be used in QITE’s linear system assembly and related calculations.
