# MT-QITE: Classical Emulation of Multiple-Time Quantum Imaginary-Time Evolution

For a detailed theoretical introduction and examples see:

 J. Del Castillo, M. Granath, E. van Nieuwenburg<br>
 Multiple-Time Quantum Imaginary Time Evolution arXiv:2512.10875 [quant-ph]<br>
 https://arxiv.org/abs/2512.10875<br>

This repository implements a **high-performance, noiseless classical simulation** of the **Multiple-Time Quantum Imaginary-Time Evolution (MT-QITE)** algorithm, based on the QITE algorithm, and several related variants.  
It is built entirely in MATLAB and designed for **research-grade experimentation**, focusing on:

- Efficient symbolic Pauli-string algebra  
- Cached expectation-value computation  
- Partitioned Hamiltonians and flexible ansätze  
- First-order, second-order, and hybrid QITE update rules  
- Detailed measurement accounting  
- Exact eigenvalue references and ITE comparison curves  

The codebase is modular and extensible, allowing new Hamiltonians, ansätze, and models to be added easily using simple master data files and a JSON configuration system.

---

# ✨ Features

- Fast expectation-value evaluation without \(2^n \times 2^n\) matrices  
- Automatic caching of Pauli measurement results  
- Arbitrary Hamiltonian partitioning  
- Fermionic and spin model support  
- Custom or pseudoinverse-based QITE linear-system solver  
- Reusable measurement structures for efficient Trotter scans  
- Exact diagonalization for reference energies and eigenstates  
- Full reconstruction of QITE circuits (operator sequences + parameters)  
- Measurement-cost profiling and Pauli-string usage statistics  

---

# 📂 Project Structure

```
master/           % Hardcoded model definitions
params/           % JSON configuration + domain definitions
pauliDictionary.m % Core Pauli infrastructure
pauliString.m     % Symbolic Pauli-string algebra
hamiltonian.m     % Builds partitioned Hamiltonians + ansätze
opsCoefficient.m  % Builds {-t_i , t_j} operator pool
opsIndependent.m  % Builds [H , t_i] and [H² , t_i] pools
qiteKernel.m      % Linear system assembly + evolution steps
qite.m            % Full QITE algorithm (Trotter scans, energy estimation)
run_qite.m        % Example driver script
```

---

# 🧠 Core Concepts

## Pauli Dictionary  
`pauliDictionary.m` efficiently represents Pauli maps so expectation values can be computed as element-wise multiplications instead of full matrix operations.

It also handles:
- Pauli-string encoding (`I/X/Y/Z`)
- Measurement caching with freshness flags
- Tracking how many *distinct* Pauli measurements were required

---

## Pauli String Algebra  
`pauliString.m` implements:
- Addition, subtraction, multiplication  
- Real/imag decomposition  
- Commutator and anticommutator  
- Efficient map computation  
- Expected value evaluation via dictionary lookups  
- Conversion to full matrices (only when necessary)

---

## Hamiltonian Construction  
`hamiltonian.m` loads model specifications from:

- **Hamiltonian** master file (`HH`)  
- **Ansatz** master file (`AZ`)  
- **Domain** file (`DD`)  
- **Parameter** file (`PP`) when applicable  
- **cfg.json** to map model names to files

It supports:
- Fermionic and spin models  
- Arbitrary partitioning (`{[1,3], [2,4,5], ...}`)  
- Different ansatz assignments per partition term  

---

## Operator Pools for QITE  
QITE requires evaluating expressions like:

- \([H, t_i]\)
- \([H^2, t_i]\)
- \(-\{t_i, t_j\}\)

These are precomputed using two support classes:

- `opsIndependent.m` → builds all `[H, t_i]` and `[H², t_i]` terms  
- `opsCoefficient.m` → builds all `- {t_i , t_j}` terms  

Each result is a `pauliString` with a cached map.

---

# ⚙️ QITE Kernel

`qiteKernel.m` performs the heavy lifting:

- Measures all operators needed for the linear system  
- Reuses measurements across time steps when allowed  
- Builds the independent term \(b(Dt)\)  
- Solves the linear system using:
  - Moore–Penrose pseudoinverse, or
  - `custom_pinv` (sparse solution for reduced circuit depth)
- Evolves the state with:
  \[
  |\psi'\rangle = \prod_i e^{-Dt\,\theta_i t_i}\,|\psi\rangle
  \]

It maintains internal structures to track:

- Coefficients (`c1`, `c2`)  
- Independent vectors (`Mb1`, `Mb2`)  
- Symmetric system matrix (`SST`)  
- Measurement reuse buffers  
- Measurement counts  

---

# 🚀 QITE Algorithm (High-Level)

`qite.m` manages:

1. **Scanning over time step sizes**
2. **Performing Trotter steps**
3. **Running QITE for each Hamiltonian term**
4. **Accumulating circuit definitions**  
5. **Computing final energies**
6. **Computing exact ITE reference**
7. **Fidelity evaluation vs exact eigenstates**
8. **Cost accounting (measurements, depth, Pauli usage)**

It returns:

- `summary` struct (energy, fidelity, accuracy, cost, optimal time point)
- End state  
- Full operator-level QITE circuit  

---

# ▶️ Example Usage

```matlab
% run_qite.m
noQubits = 4;
modelFam = "H2";

PDH = pauliDictionary(noQubits);
partition = {[1:3], [4:15]};

options.order = "2nd";
options.Pinv  = false;

initialState = PDH.prepare_initial_state(modelFam);

QK = qiteKernel(PDH, options, modelFam, 1, partition{:});
Q  = qite(QK, initialState, 0.05:0.05:0.50, 10);

[Q, summary, endState] = Q.estimate_energy();
```

---

# 📊 Outputs

The algorithm reports:

- **Ground-state energy estimate**
- **Deviation from exact eigenvalue**
- **Fidelity w.r.t. minimum eigenvector(s)**
- **Circuit depth (ops + Pauli strings)**
- **Total measurement cost**
- **Optimal timestep and Trotter step**

A final plot compares:

- QITE energy vs. imaginary-time β
- Exact ITE energy curve

---

# 🗂 Adding New Models

Add to `params/cfg.json`:

```json
{
  "model": "MyModel",
  "type": "fermions",
  "hamiltonian": "/_HHF04.m",
  "ansatz": "/_AZF04.m",
  "domain": "/_DDF04.m",
  "parameter": "/_PPF04.m"
}
```

Then create the master data files:

- `master/_HHF04.m`
- `master/_AZF04.m`
- `master/_DDF04.m`
- `master/_PPF04.m` (optional)

Each should define the appropriate arrays of `pauliString` objects.

---

# 🔧 Requirements

- MATLAB R2023a or newer  
- Statistics/Optimization toolbox recommended for pseudoinverse  
- No external dependencies  

---

# 🧪 Testing

You can validate correctness by checking:

- Hermiticity of Hamiltonian terms (done automatically)
- Anti-Hermiticity of ansatz terms (checked automatically)
- Comparison of QITE vs exact ITE curves
- Fidelity vs exact eigenvectors

---

# 📜 License

MIT License

---

## Citation

If you use this code in academic work, please cite the associated paper:

> Julio Del Castillo, Mats Granath, Evert van Nieuwenburg,  
> **Multiple-Time Quantum Imaginary Time Evolution**,  
> *arXiv preprint* (2025).  
> https://arxiv.org/abs/2512.10875

---

# 🙌 Acknowledgements

This project draws inspiration from:

- The original QITE algorithm (Motta et al., 2020)  
- Classical imaginary-time simulations  
- Quantum chemistry ansatz structures  

And was developed with careful focus on **performance**, **correctness**, and **readability**.

### BibTeX
```bibtex
@article{Del_Castillo_Multiple-Time_Quantum_Imaginary_2025,
author = {Del Castillo, Julio and Granath, Mats and van Nieuwenburg, Evert},
journal = {arXiv},
title = {{Multiple-Time Quantum Imaginary Time Evolution}},
url = {https://arxiv.org/abs/2512.10875},
year = {2025}
}

