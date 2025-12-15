# Reproducing the Published Results

This repository accompanies the results reported in the preprint  
**arXiv:2512.10875**.

The numerical data and figures presented in the manuscript can be reproduced using the MATLAB source files provided here, following the instructions below.

---

## Directory Structure

To reproduce the results, place all available source files in a common root directory with the following structure:

/src  
/examples  
/params  
/master  
/initials  
/results  

The role of each directory is as follows:

- **/src**  
  Core MATLAB source files implementing the QITE and MT-QITE algorithms.

- **/examples**  
  Driver scripts (`reported_XXXX.m`) corresponding to specific figures or datasets reported in the manuscript.

- **/params**  
  Configuration files (JSON and auxiliary `.m` files) defining model parameters and runtime options.

- **/master**  
  Hard-coded Hamiltonian, ansatz, and domain definitions used to construct the physical models.

- **/initials**  
  Batches of initial state vectors used when multiple initializations are required.

- **/results**  
  Precomputed `.mat` files containing the data used to generate the figures, provided for convenience.

---

## Running the Simulations

From the root directory, execute the desired `reported_XXXX.m` script corresponding to the figure or result you wish to reproduce. These scripts will:

1. Initialize the Hamiltonian and ansatz definitions.
2. Run the QITE or MT-QITE simulation.
3. Generate intermediate `.mat` data files when required.
4. Produce the published figures.

All plotting is handled by dedicated MATLAB scripts named `qiteGraphs*.m`.

---

## Notes on Computational Cost

- Some simulations, in particular those involving **8-qubit models**, may require substantial computational time.
- Depending on your hardware, additional partitioning or batching strategies may be necessary.
- In the specific case of the **H₄ chain**, data generation and figure rendering are intentionally separated into two distinct scripts.

---

## Reproducibility

The provided `.mat` files in `/results` correspond exactly to the data used in the manuscript figures. Users may either regenerate the data from scratch or directly load these files to reproduce the plots.

For questions or clarifications regarding reproducibility, please refer to the accompanying arXiv manuscript.
