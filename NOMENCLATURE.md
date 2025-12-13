# Nomenclature and Naming Conventions

This document describes the nomenclature policy used throughout the codebase. It has been observed whenever possible except for pauliDictionary.m and pauliString.m which are legacy classes from an older project. 
The goal is to ensure **legibility, consistency, and ease of navigation**, especially for algorithmic and research-oriented code.

---

## General Naming Style

- **Variables, vectors, and objects** use **Java-style camelCase** notation  
  Example: `initialState`, `timeSeries`, `measureCountLS`

- **Functions and methods** use **snake_case notation**  
  Example: `estimate_energy()`, `scan_step()`, `calculate_circuit()`

- **Dimensions of arrays, structs, or tables** are explicitly named using  
  `size_variableName`  
  Example: `size_vector`, `size_perm`, `size_grid`

This explicit naming avoids ambiguity and improves readability when handling multidimensional data.

---

## Prefixes for Disambiguation

To clarify the nature of composite variables, the following prefixes are used:

- `cellXXX` : cell arrays  
- `arrXXX`  : numeric arrays  
- `stuXXX`  : structures  
- `tabXXX`  : tables  

Example:
- `cellQiteOps`
- `arrStuSummary`
- `stuSummary`

---

## Boolean Variables

- Boolean (logical) variables are always prefixed with:

  `flag`

Example:
- `flagConverged`
- `flagSymmetryEnabled`

---

## Indices and Counters

- Indices and loop counters are named with the prefix:

  `id`

This makes indexing intent immediately clear.

Example:
- `id`
- `idGrid`
- `idTerm`

---

## Standard Abbreviations

The following abbreviations are used consistently in variable and function names:

- `LS`     : Linear System  
- `AZ`     : Ansatz  
- `Op`, `Ops` : Operator(s), operation(s)  
- `PS`     : Pauli String  
- `params` : Parameters  
- `cfg`    : Configuration  
- `ini`    : Initial  
- `perm`   : Permutation  
- `acc`    : Accumulated  

Example:
- `measureCountLS`
- `cellQiteOps`
- `perm`
- `accEnergy`

---

## Design Philosophy

This nomenclature policy is designed to:
- Make mathematical intent explicit
- Reduce cognitive load when navigating complex algorithms
- Facilitate collaboration and long-term maintainability
- Align code structure with the conceptual structure of the underlying physics

Consistency is preferred over brevity: names are chosen to be readable and unambiguous rather than minimal.
