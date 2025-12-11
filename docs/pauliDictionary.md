*pauliDictionary* Class Documentation

## Overview

The *pauliDictionary* *handle* class provides all Pauli-string infrastructure required by the QITE-based algorithms and related simulation routines. Its responsibilities include:

- defining one-qubit Pauli matrices 𝐼, 𝑋, 𝑌, 𝑍
and additional one-qubit operators,

- converting between integer indices and Pauli strings,

- computing tensor products of operator factors,

- storing and reusing efficient “Pauli maps” for expectation-value evaluation,

- caching measurement results of Pauli strings for a given state,

- keeping counters of measurement usage during QITE,

- providing convenient predefined initial states for supported model families.

Wavefunctions are represented as complex MATLAB vectors of size $2^𝑛$.

Pauli strings (e.g. "XXYY") are stored in compact symbolic form, and their action on a wavefunction is encoded via maps with entries in {+1,−1,+1𝑖, −1𝑖}, exploiting the structure of tensor-product Pauli operators and avoiding explicit $2^n x 2^n$ matrices.

This class is used by pauliString, qiteKernel, and ultimately the high-level QITE-based algorithms.

## Constructor
PDH = pauliDictionary(nqubits)

Creates a dictionary object for an 𝑛-qubit system.

### Properties initialized:

- nqbs — number of qubits
- twoN = $2^n$ — wavefunction dimension
- fourN = $4^n$ — size of full Pauli basis
- Pauli matrices 𝐼, 𝑋, 𝑌, 𝑍

Additional 1-qubit operators:

- *A* — Hadamard gate
- *C* — creation operator
- *D* — annihilation opeator
- *N* — number operator

Dictionaries (accessible througout the whole project):
- DIC — dictionary of Pauli maps (used to compute expectation values efficiently)
- DICOM — dictionary of measurement results and “fresh” flags

Measurement counters
- Internal counters measuresLS, measuresCO, measuresAZ

## Pauli-String Encoding (index to string)
s = numToString(i)

Converts an integer index i (1-based) into the corresponding Pauli string of length nqbs.
The index is interpreted in base-4, where:

| Digit | Pauli |
|-------|--------|
|   0   |   I    |
|   1   |   X    |
|   2   |   Y    |
|   3   |   Z    |

<u>Examples (for 3 qubits):</u>

PDH.numToString(1)   % returns 'III'

PDH.numToString(2)   % returns 'IIX'

PDH.numToString(16)  % returns 'ZZZ'

## Pauli-String Encoding (string to index)
n = stringToNum(s)

Inverse mapping: converts a Pauli string s into its integer index.

<u>Example:</u>

PDH.stringToNum('IIX')   % returns 2

Indexing is again 1-based.

## Tensor Products
K = kp(MultsCell)

Computes the Kronecker (tensor) product of a row cell array of square matrices:

K = MultsCell{1} ⊗ MultsCell{2} ⊗ ... ⊗ MultsCell{end}


Restrictions:

MultsCell must be 1×N.

Used internally to build explicit operators when needed (rarely used in the QITE flow thanks to Pauli maps).

## Measurement Cache and Counters

Expectation values of Pauli strings can appear many times during QITE. To avoid redundant work, *pauliDictionary* caches:

- the computed expectation value (calculated by class *pauliString*),
- whether the value is still “fresh” for the current state.

These are stored in DICOM.

Counters track how many distinct measurements were actually needed:

measuresLS — measurements for linear-system assembly

measuresCO — for coefficient-container generation

measuresAZ — for multi-sweep analysis operations

### Reset and increment methods:

resetMeasuresLS, resetMeasuresCO, resetMeasuresAZ

incMeasuresLS, incMeasuresCO, incMeasuresAZ

#### flushMeasurements()

Marks all cached measurements as not fresh.
This should be called after the wavefunction is updated so that expectation values are recomputed only when necessary.

## Model-Specific Initial States for a selction of models
state = prepare_initial_state(modelFam)

Returns a normalized initial wavefunction for a number of predefined model families.

This is useful for running reproducible QITE benchmarks.

Supported combinations include (examples):

| Model family                        | Qubits | Description                               |
|-------------------------------------|--------|-------------------------------------------|
| "H2"                                |   4    | Molecular hydrogen HF initial state   |
| "ISIOB04"                           |   4    | Ising family all-ones state           |
| "ISIOB06", "ISYOB06", "HEYOB06"     |   6    | Various 6-qubit benchmark states           |
| "HEYOB08"                           |   8    | Strongly structured 8-qubit benchmark     |
| "H4CHAIN..."                        |   8    | 4-electron H chain (Hartree–Fock)           |

If the pair (modelFam, nqbs) is unknown, an error is thrown.

<u>Example:</u>

PDH = pauliDictionary(4);

psi0 = PDH.prepare_initial_state("H2");

## Internal Properties
| Property            | Description                                                 |
|---------------------|-------------------------------------------------------------|
| **PAULI**           | Cell array `{I, X, Y, Z}`                                   |
| **DIC**             | Dictionary mapping Pauli strings to their action maps       |
| **DICOM**           | Dictionary mapping Pauli strings to cached measurement results |
| **measuresLS/CO/AZ** | Measurement counters (for linear system, energy, variational sweeps) |

## Relationship to Other Classes

pauliString uses pauliDictionary to store Pauli maps and measurements.


## Summary

*pauliDictionary* is the foundational layer for all Pauli-operator manipulation in this QITE and QITE-related implementations:

- It provides Pauli-string encoding/decoding.

- It enables fast expectation-value computation.

- It caches measurements to dramatically reduce repeated work.

- It supports reproducible initial states for several quantum-simulation models.

Together with *pauliString*, it forms an efficient symbolic layer on top of MATLAB vectors, avoiding full matrix representations and enabling scalable classical emulation of QITE-based algorithms.