## Example: Running QITE on a 4-qubit H₂ Hamiltonian

```matlab
% Sample script to run QITE
clear;
tic;

noQubits = 4;
modelFam = "H2";

% Create Pauli dictionary (manages Pauli string maps & cached measurements)
PDH = pauliDictionary(noQubits);

% Partition of Hamiltonian into two terms
partition = {[1:3], [4:15]};

% QITE kernel options
options.order = "2nd";   % use 2nd–order formulas
options.Pinv  = false;   % do not use Moore-Penrose algorithm

% Prepare initial state using predefined model
initialState = PDH.prepare_initial_state(modelFam);

% Build QITE kernel (one ansatz per Hamiltonian term acording to local domains)
QK = qiteKernel(PDH, options, modelFam, 1, partition{:});

% Construct QITE object with time steps and Trotter steps
Q = qite(QK, initialState, 0.05:0.05:0.50, 10);

% Run QITE and estimate the energy
[Q, summary, endState] = Q.estimate_energy();

totalTime = toc;
