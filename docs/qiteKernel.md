# `qiteKernel` Class Documentation

## Overview

The `qiteKernel` class implements the core numerical engine of the QITE algorithm. It:

- Inherits from `hamiltonian`, gaining:
  - partitioned Hamiltonian terms `H(m)`
  - total Hamiltonian `HT`
  - an ansatz per partition `cellAZ{m}`
- Precomputes and stores all operator pools required to assemble QITE linear systems:
  - `[H, t_I]` via `opsIndependent`
  - `[H^2, t_I]` via `opsIndependent` (second-order formulas)
  - `-{t_I, t_J}` via `opsCoefficient`
- For a given reference state:
  - measures all necessary expectation values (`measure_state`)
  - builds the linear-system vectors (`independent_term`)
  - solves the linear system for a timestep `Dt` (`solve_linear_system`)
  - produces QITE parameters (`thetas`) and corresponding ansatz operators (`qiteOps`)
  - can evolve the state using these parameters (`evolve_state`)

The class supports `"1st"`, `"2nd"` and `"1st0"` QITE formula orders and allows choosing between MATLAB’s `pinv` or a custom sparse solver.

---

## Constructor

```matlab
QK = qiteKernel(PDH, options, model, p, varargin)
```

### Arguments

| Argument   | Type                     | Description |
|------------|--------------------------|-------------|
| `PDH`      | `pauliDictionary` handle | Shared Pauli dictionary |
| `options`  | struct                   | Must contain `.order` and `.Pinv` |
| `model`    | string                   | Hamiltonian model name (must exist in cfg.json) |
| `p`        | scalar                   | Model parameter (use 1 if not needed) |
| `varargin` | vectors                  | Partition groups of Hamiltonian terms |

### `options` fields

- `options.order` = `"1st"`, `"2nd"`, or `"1st0"`
- `options.Pinv` = `true` (use MATLAB `pinv`) or `false` (use sparse custom solver)

---

## Inherited and Added Properties

Inherited from `hamiltonian`:  
`EPSEIG`, `PDH`, `type`, `model`, `H`, `M`, `HT`, `cellAZ`.

Additional:

- `EPSLS` – tolerance for low-rank linear systems  
- `stateVector` – reference state used for all measurements  
- `H2` – squared Hamiltonian terms (used only in second-order QITE)  
- `measured(m)` – struct storing measurement results for partition term `m`:  
  - `c1`, `c2` – normalization coefficients  
  - `Mb1`, `Mb2` – independent-term vectors  
  - `SST` – symmetric LS coefficient matrix  
- `measureCount(m)` – number of Pauli-string measurements performed for term `m`  
- `tIH(m)` – `[H(m), t_I]` commutators as `opsIndependent`  
- `tIH2(m)` – `[H(m)^2, t_I]` commutators (`"2nd"` order only)  
- `tItJ(m)` – `-{t_I, t_J}` operators as `opsCoefficient`

---

## Precomputation in the Constructor

After calling the parent `hamiltonian` constructor:

- Validate `options.order` and `options.Pinv`
- Precompute:
  - `HT = HT.calcMap()`
  - `H(m) = H(m).calcMap()`
  - If second-order: `H2(m) = (H(m)^2).calcMap()`
- Build operator pools:
  - `tIH(m)  = opsIndependent(PDH, cellAZ{m}, H(m))`
  - `tIH2(m) = opsIndependent(PDH, cellAZ{m}, H2(m))`  (only for `"2nd"`)
  - `tItJ(m) = opsCoefficient(PDH, cellAZ{m})`

All operators are immediately mapped for fast expectation evaluation.

---

## Setting the Reference State

```matlab
QK = QK.set_statevector(newState);
```

- `newState` must be `2^n × 1`.
- The state is normalized for safety.

---

## Measuring a Partition Term

```matlab
QK = QK.measure_state(mterm);
```

For the partition term `mterm = 1..M`, this method:

1. Initializes measurement storage (`c1`, `c2`, `Mb1`, `Mb2`, `SST`)
2. Computes energy-related coefficients:
   - `"1st0"`: `c1 = 0`, `c2 = 0`
   - `"1st"`:  `c1 = Re(<H>)`, `c2 = 0`
   - `"2nd"`:  `c1 = Re(<H>)`, `c2 = Re(<H^2>)`
3. Computes:
   - `Mb1(i) = Re(<[H, t_I]>)`
   - `Mb2(i) = Re(<[H^2, t_I]>)` (second order only)
4. Builds the LS matrix:
   - `SST(i,j) = Re(< -{t_I, t_J} >)`
   - Lower triangle filled by symmetry
5. Updates measurement count:
   - `measureCount(mterm) = PDH.measuresLS - start_count`

This method is called whenever the reference state changes.

---

## Building the Linear-System Independent Term

```matlab
b = QK.independent_term(mterm, Dt);
```

With timestep `Dt`, define:

```
c = 1 - 2*Dt*c1 + 2*Dt^2*c2
```

Then:

```
b = c^(-1/2) * ( Mb1 - 0.5*Dt*Mb2 )
```

For `"1st"` and `"1st0"`, the expression reduces automatically since `Mb2 = 0`.

---

## Solving the Linear System

```matlab
[thetas, qiteOps] = QK.solve_linear_system(mterm, Dt);
```

Steps:

1. Build `b = independent_term(mterm, Dt)`
2. Solve:

   - With Moore–Penrose:

     ```
     X = pinv(SST) * b
     ```

   - With custom sparse solver:

     ```
     X = custom_pinv(SST, b, EPSLS)
     ```

3. Keep only significant components:

```
id = abs(X) > EPSLS
thetas  = X(id)'                % row vector
qiteOps = cellAZ{mterm}(id')    % row array of pauliString
```

Results define the QITE unitary update for that term.

---

## Evolving the State

```matlab
psi_new = QK.evolve_state(Dt, thetas, qiteOps);
```

For each `(theta_i, op_i)`:

1. Convert to matrix: `M = op_i.matrixRep()`
2. Build unitary: `U = expm(-Dt * theta_i * M)`
3. Update:

```
state = U * state
```

The method returns the evolved state but does not modify `stateVector`; the caller (typically `qite`) is responsible for updating it.

---

## Custom Sparse Linear Solver

```matlab
X = QK.custom_pinv(A, B, tol);
```

Solves `A*X = B` by:

- QR decomposition with pivoting
- Detecting effective rank via tolerance on transformed RHS
- Solving the reduced system
- Expanding back with permutation matrix `P`

Behavior:
- Zeros as many components as allowed by tolerance → shallower QITE circuits.

---

## Interaction with Other Classes

- **`pauliDictionary`**  
  Provides Pauli maps and measurement caches; all expectation values use `eValLS`.

- **`hamiltonian`**  
  Supplies local Hamiltonian terms and ansatz partitions.

- **`opsCoefficient`**  
  Provides `-{t_I, t_J}` operator pool.

- **`opsIndependent`**  
  Provides `[H, t_I]` and `[H^2, t_I]` operator pools.

- **`pauliString`**  
  Used for all operators; provides commutators, anticommutators, mapping, matrix representation.

- **`qite`**  
  High-level driver handling:
  - Trotter steps
  - Time steps
  - Measurement reuse
  - Energy tracking
  - State updates

---

## Typical Workflow

```matlab
PDH = pauliDictionary(noQubits);
partition = { [1:3], [4:15] };

options.order = "2nd";
options.Pinv  = false;

QK = qiteKernel(PDH, options, "H2", 1, partition{:});

psi0 = PDH.prepare_initial_state("H2");
QK   = QK.set_statevector(psi0);

QK = QK.measure_state(1);

Dt = 0.05;
[thetas, qiteOps] = QK.solve_linear_system(1, Dt);

psi1 = QK.evolve_state(Dt, thetas, qiteOps);
```

This is the computation performed repeatedly by the `qite` class for full QITE simulation.

