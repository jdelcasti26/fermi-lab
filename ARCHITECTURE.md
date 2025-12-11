# QITE Project Architecture Diagram

Below is a complete high-level architecture diagram showing how all classes, data files, and processes interact in this implementation of QITE.

```
                                    ┌──────────────────────────────┐
                                    │          User Script         │
                                    │        (run_qite.m)          │
                                    └───────────────┬──────────────┘
                                                    │
                                                    ▼
                                    ┌──────────────────────────────┐
                                    │      pauliDictionary.m       │
                                    │  - Defines I,X,Y,Z matrices  │
                                    │  - Pauli-string maps (Idx/TR/TI)
                                    │  - Measurement cache (DICOM) │
                                    │  - Measurement counters      │
                                    └───────────────┬──────────────┘
                                                    │
                                                    ▼
                    ┌────────────────────────────────────────────────────────┐
                    │                     pauliString.m                      │
                    │ - Symbolic Pauli algebra (+, -, *, comm, anticomm)    │
                    │ - Efficient map computation using pauliDictionary      │
                    │ - Converts strings to matrices when needed             │
                    │ - Expected-value evaluation via cached maps            │
                    └─────────────────────────────────┬──────────────────────┘
                                                      │
                                                      ▼
        ┌──────────────────────────────────────────────────────────────────────────────┐
        │                               hamiltonian.m                                  │
        │ - Loads model specification from params/cfg.json                              │
        │ - Reads master data files:                                                    │
        │      HH : Hamiltonian Pauli strings                                           │
        │      AZ : Ansatz Pauli strings                                                │
        │      DD : Domain masks / ansatz selectors                                     │
        │      PP : Optional parameter scaling                                          │
        │ - Partitions Hamiltonian via varargin arguments                                │
        │ - Produces:                                                                   │
        │      H(m)   : PauliString for each partition term                             │
        │      cellAZ : Ansatz per partition                                            │
        │      HT     : Total Hamiltonian                                               │
        └───────────────────────────────┬───────────────────────────────────────────────┘
                                        │
                                        ▼
              ┌──────────────────────────────────────────────────────────────────────┐
              │                             qiteKernel.m                             │
              │ Inherits from hamiltonian.m                                          │
              │                                                                      │
              │ Precomputation phase:                                                │
              │   - Compute H(m).map and H²(m).map                                   │
              │   - Build operator pools:                                            │
              │        opsIndependent: [H , t_i], [H² , t_i]                         │
              │        opsCoefficient:  -{t_i , t_j}                                 │
              │                                                                      │
              │ Measurement phase:                                                   │
              │   - measure_state():                                                 │
              │        c1, c2 normalizers                                            │
              │        Mb1, Mb2 independent vectors                                  │
              │        SST symmetric matrix                                           │
              │                                                                      │
              │ Linear system:                                                       │
              │   - independent_term(Dt)                                             │
              │   - solve_linear_system():                                           │
              │        pinv   or   custom_pinv                                       │
              │   → Returns theta parameters + operators for exponentials            │
              │                                                                      │
              │ Evolution: evolve_state(Dt, thetas, ops)                             │
              │                                                                      │
              └───────────────────────┬──────────────────────────────────────────────┘
                                      │
                                      ▼
                    ┌─────────────────────────────────────────────────────────┐
                    │                          qite.m                         │
                    │                                                         │
                    │ Controls the *full algorithm*:                          │
                    │   - Loop over timeSeries                                │
                    │   - Loop over Trotter steps                             │
                    │   - For each term: measure → solve → evolve             │
                    │   - Collect circuits (ops + thetas)                     │
                    │   - Measurement cost accounting                         │
                    │                                                         │
                    │ estimate_energy():                                      │
                    │   - Executes full QITE scan                             │
                    │   - Reconstructs circuits                               │
                    │   - Computes energies, fidelities                       │
                    │   - Produces summary and final state                    │
                    │   - Compares with exact ITE and exact eigenvalues       │
                    └─────────────────────────────────────────────────────────┘
```

---

# Data Flow Summary

```
cfg.json ──► hamiltonian.m
HH.m    ───► H(m) terms
AZ.m    ───► Ansatz operators
DD.m    ───► Domain → selects subset of AZ per term
PP.m    ───► Optional scaling of HH terms
```

```
hamiltonian.m → qiteKernel.m → qite.m → energy curves, circuit depth, fidelity, final state
```

---

# Measurement Reuse Logic

```
Trotter step ts = 1, partition term m = 1:
    First time Dt appears:
        measure_state()
        reuse = measured
Later:
    measured = reuse  (no new measurements)
```

After evolution:

```
PDH.flushMeasurements()   → marks all cached values non-fresh
```

---

# Operator Pools (used by qiteKernel)

```
opsIndependent:
    [H , t_i]
    [H² , t_i]   (if 2nd order)

opsCoefficient:
    - {t_i , t_j}
```

Each stored as a pauliString with map precomputed.

---

# QITE Linear System Structure

For a given partition term:

```
SST * theta = b
```

where:

- `SST(i,j)` = measurement of `- {t_i, t_j}`
- `b` built from:
  - `Mb1(i)` = measurement of [H, t_i]
  - `Mb2(i)` = measurement of [H², t_i] (2nd order)
  - time-dependent normalization c(Dt)

Solution:

```
theta = pinv(SST) * b
```

or sparse solution via custom_pinv.

---

# Evolution Pipeline

```
Initial state ψ₀
→ solve linear system
→ exponentiate operators exp(-Dt θ_i t_i)
→ apply sequentially
→ new state ψ₁
→ repeat for next Trotter step
```

---

If you'd like a **graphical version** (ASCII or Mermaid flowchart), I can generate that too.

Ready for item **(2) dependency graph** when you are!
