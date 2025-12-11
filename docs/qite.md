# `qite` Class Documentation

## Overview

The `qite` class is the **high-level driver** that runs the QITE algorithm using a `qiteKernel` instance.

It:

- Holds:
  - A handle to the QITE kernel (`QK`)
  - The **initial state**
  - A set of **time step sizes** to scan (`timeSeries`)
  - A **number of Trotter steps** (`noTS`)
- For each time step size and Trotter depth:
  - Runs the full QITE procedure (`scan`)
  - Stores all generated **circuits** (operators + parameters)
  - Computes **energy estimates** and **costs** (`estimate_energy`)
  - Tracks **operator depth** and **Pauli-string** depth
  - Computes **fidelity** with the exact ground state
- Produces:
  - A `stuSummary` structure with performance metrics
  - The **final optimal state**
  - A plot of QITE vs exact imaginary-time evolution

---

## Constructor

```matlab
Q = qite(QK, initialState, timeSeries, noTS)
```

### Arguments

| Argument      | Type                     | Description                                                  |
|---------------|--------------------------|--------------------------------------------------------------|
| `QK`          | `qiteKernel`             | QITE kernel (constructed with a given model, partition, etc.)|
| `initialState`| column vector            | Initial state (size `2^n × 1`)                              |
| `timeSeries`  | row vector               | Positive time step sizes \(\Delta\tau\) to scan             |
| `noTS`        | integer                  | Number of Trotter steps                                      |

### Behavior

- Verifies `initialState` has shape `(2^n, 1)` where `n` is the number of qubits in `QK.PDH`.
- Normalizes `initialState` for safety.
- Verifies `timeSeries` is:
  - a row vector
  - strictly positive (`timeSeries > 0`)
- Sorts and deduplicates `timeSeries`:

```matlab
o.timeSeries = sort(unique(timeSeries));
o.noTIMES    = size_timeSeries(2);   % number of time step sizes
o.noTS       = noTS;
```

- Allocates containers:

  ```matlab
  o.cellQiteOps    = cell(o.noTIMES, o.noTS, o.QK.M);
  o.cellThetas     = cell(o.noTIMES, o.noTS, o.QK.M);
  o.measureCountLS = zeros(o.noTIMES, o.noTS, o.QK.M);
  o.measureCountEE = zeros(o.noTIMES, o.noTS);
  ```

---

## Properties

```matlab
classdef qite
    properties
        % Configuration
        QK            % qiteKernel handle
        initialState  % normalized initial state
        timeSeries    % row vector of time steps
        noTIMES       % number of distinct time steps
        noTS          % number of Trotter steps

        % Simulation results
        cellQiteOps      % {tn, ts, mterm}: ansatz operators (pauliString row array)
        cellThetas       % {tn, ts, mterm}: corresponding parameters (row vector)
        measureCountLS   % linear-system measurement counts (tn, ts, mterm)
        measureCountEE   % energy measurement counts (tn, ts)
        arrayE           % estimated energies (tn, ts)
        arrOpCount       % circuit operator counts (tn, ts)
        arrPSCount       % circuit Pauli-string counts (tn, ts)
    end
end
```

---

## Running the QITE Scans

### `scan` method

```matlab
Q = Q.scan();
```

This method runs the **QITE procedure** for:

- each time step size in `timeSeries`, and
- each Trotter step up to `noTS`, and
- each partition term `mterm = 1..QK.M`.

#### Workflow

For each time step `Dt` (indexed by `tn`):

1. Reset the kernel state:

   ```matlab
   o.QK = o.QK.set_statevector(o.initialState);
   o.QK.PDH = o.QK.PDH.flushMeasurements();
   ```

2. For each Trotter step `ts`:
3. For each partition term `mterm`:

   - Special treatment for `ts == 1` and `mterm == 1` at the **first time**:

     ```matlab
     if ts == 1 && mterm == 1
         if Dt == o.timeSeries(1)
             % Measure LS once and reuse for all Dt
             o.QK = o.QK.measure_state(mterm);
             o.QK.reuse = o.QK.measured;
             o.measureCountLS(tn, ts, mterm) = o.QK.measureCount(mterm);
         end

         % Reuse previously measured data
         o.QK.measured = o.QK.reuse;
     else
         % For other (ts, mterm), measure freshly
         o.QK = o.QK.measure_state(mterm);
         o.measureCountLS(tn, ts, mterm) = o.QK.measureCount(mterm);
     end
     ```

   - Solve the linear system:

     ```matlab
     [thetas, qiteOps] = o.QK.solve_linear_system(mterm, Dt);
     ```

   - Store circuit parameters:

     ```matlab
     o.cellQiteOps{tn, ts, mterm} = qiteOps;
     o.cellThetas{tn, ts, mterm}  = thetas;
     ```

   - Evolve the current state:

     ```matlab
     evolvedState = o.QK.evolve_state(Dt, thetas, qiteOps);
     o.QK = o.QK.set_statevector(evolvedState);
     o.QK.PDH = o.QK.PDH.flushMeasurements();
     ```

This fills `cellQiteOps`, `cellThetas`, and `measureCountLS` for all combinations of `(time step, Trotter step, partition term)`.

---

## Estimating Energy and Summarizing Results

### Signature

```matlab
[Q, stuSummary, finalState] = Q.estimate_energy();
```

This is the **main user-facing method**.

It:

1. Allocates result arrays:

   ```matlab
   o.arrayE     = zeros(o.noTIMES, o.noTS);
   o.arrOpCount = zeros(o.noTIMES, o.noTS);
   o.arrPSCount = zeros(o.noTIMES, o.noTS);
   ```

2. Runs the full scan:

   ```matlab
   o = o.scan();
   ```

3. For each time index `tn` and Trotter step `ts`:
   - Builds the circuit up to step `ts` (using `calculate_circuit`)
   - Applies it to `initialState`
   - Measures the energy with `QK.HT.eValCO(state)`
   - Counts energy measurements and circuit depth

#### Energy evaluation

For each `(tn, ts)`:

```matlab
state = o.calculate_circuit(tn, ts) * o.initialState;

o.QK.PDH = o.QK.PDH.flushMeasurements();
measures_ref = o.QK.PDH.measuresCO;

o.arrayE(tn, ts) = o.QK.HT.eValCO(state);
o.measureCountEE(tn, ts) = o.QK.PDH.measuresCO - measures_ref;
```

#### Circuit depth counting

For each `(tn, ts)`:

```matlab
for mterm = 1:o.QK.M
    size_Ops = size(o.cellQiteOps{tn, ts, mterm}, 2);
    for i = 1:size_Ops
        o.arrOpCount(tn, ts) = o.arrOpCount(tn, ts) + 1;
        o.arrPSCount(tn, ts) = o.arrPSCount(tn, ts) + ...
            o.cellQiteOps{tn, ts, mterm}(i).len;
    end
end
```

`arrOpCount` gives total number of exponentials;  
`arrPSCount` counts the Pauli strings involved.

---

## Summary Structure

The method builds a summary struct:

```matlab
stuSummary.energy     % best energy per Trotter step
stuSummary.accuracy   % relative error vs exact ground energy
stuSummary.fidelity   % fidelity vs ground-state subspace
stuSummary.opCount    % total operator count (up to that Trotter step)
stuSummary.PSCount    % total Pauli-string count (up to that Trotter step)
stuSummary.minPosition% time index at which minimum energy was found
stuSummary.costLS     % cumulative LS measurement cost
stuSummary.costEE     % cumulative energy measurement cost
```

### Reference exact values

```matlab
exactEigenvalue   = o.QK.exact_eigenvalue();
exactEigenvectors = o.QK.exact_eigenvectors();
```

For each Trotter step `ts`:

1. **Best energy** and its time index:

   ```matlab
   [stuSummary.energy(ts), atrow] = min(o.arrayE(:, ts));
   stuSummary.minPosition(ts) = atrow;
   ```

2. **Accumulated circuit depth** up to that Trotter step:

   ```matlab
   stuSummary.opCount(ts) = sum(o.arrOpCount(atrow, 1:ts));
   stuSummary.PSCount(ts) = sum(o.arrPSCount(atrow, 1:ts));
   ```

3. **Accuracy**:

   ```matlab
   stuSummary.accuracy(ts) = (stuSummary.energy(ts) - exactEigenvalue) / ...
                             abs(exactEigenvalue);
   ```

4. **Fidelity** with exact ground-space**:

   ```matlab
   endState = o.calculate_circuit(atrow, ts) * o.initialState;
   stuSummary.fidelity(ts) = (endState' * exactEigenvectors) * ...
                             (endState' * exactEigenvectors)';
   ```

5. **Measurement costs** for that Trotter step:

   ```matlab
   stuSummary.costLS(ts) = sum(sum(o.measureCountLS(:, ts, :)));
   stuSummary.costEE(ts) = sum(o.measureCountEE(:, ts));
   ```

Then costs are **cumulatively summed** over Trotter steps:

```matlab
stuSummary.costLS = cumsum(stuSummary.costLS);
stuSummary.costEE = cumsum(stuSummary.costEE);
```

---

## Final State and Logging

After building `stuSummary`:

- The method finds overall best energy:

  ```matlab
  [minE, atstep] = min(stuSummary.energy);
  ```

- Computes final state:

  ```matlab
  finalState = o.calculate_circuit( ...
      stuSummary.minPosition(atstep), atstep) * o.initialState;
  ```

- Prints summary info to the console, including:
  - QITE ground energy
  - accuracy and fidelity
  - operator/string depth
  - effective Trotter steps
  - total measurements

---

## Exact Imaginary-Time Evolution (Reference)

For comparison, `estimate_energy` computes the exact imaginary-time expectation:

```matlab
HM = o.QK.HT.matrixRep();
ExactITE = zeros(o.noTIMES, 1);

for tm = 1:o.noTIMES
    Dtm = o.timeSeries(tm);
    num = real(o.initialState' * expm(-Dtm * HM) * HM * expm(-Dtm * HM) * o.initialState);
    den = real(o.initialState' * expm(-Dtm * HM) * expm(-Dtm * HM) * o.initialState);
    ExactITE(tm) = num / den;
end
```

Then it generates a plot of:

- `o.arrayE(:, ts)` for all Trotter steps `ts`
- `ExactITE` as reference

and labels the axes and legend accordingly.

---

## Circuit Reconstruction

### `calculate_circuit` method

```matlab
Circuit = Q.calculate_circuit(timeDatapoint, trotterStep);
```

Builds the **full QITE circuit** (as a matrix) for:

- a given time step index `timeDatapoint`, and
- up to a Trotter step `trotterStep`.

#### Algorithm

1. Start with identity:

   ```matlab
   Circuit = eye(o.QK.PDH.twoN);
   ```

2. For each Trotter step `ts = 1..trotterStep`:
3. For each partition term `mterm = 1..QK.M`:
4. For each operator in `cellQiteOps{timeDatapoint, ts, mterm}`:

   ```matlab
   opMatrix = o.cellQiteOps{timeDatapoint, ts, mterm}(i).matrixRep();
   theta    = o.cellThetas{timeDatapoint, ts, mterm}(i);

   Circuit = expm(- o.timeSeries(timeDatapoint) * theta * opMatrix) * Circuit;
   ```

This yields a **full matrix representation** of the QITE circuit, which can be applied to any state.

---

## Typical Usage Example

```matlab
% 1. Prepare kernel and qite instance
noQubits = 4;
model    = "H2";

PDH = pauliDictionary(noQubits);
partition = { [1:3], [4:15] };

options.order = "2nd";
options.Pinv  = false;

QK = qiteKernel(PDH, options, model, 1, partition{:});

psi0 = PDH.prepare_initial_state(model);
timeSeries = 0.05:0.05:0.50;
noTS = 10;

Q = qite(QK, psi0, timeSeries, noTS);

% 2. Run QITE and obtain summary and final state
[Q, summary, finalState] = Q.estimate_energy();

% 3. Inspect results
summary.energy
summary.accuracy
summary.fidelity
summary.opCount
summary.PSCount
summary.costLS
summary.costEE
```

The `qite` class thus provides the orchestration and analysis layer on top of `qiteKernel`, handling:

- time-step scanning,
- Trotter-step stacking,
- circuit construction and storage,
- measurement and depth accounting,
- and comparison against exact imaginary-time evolution.
