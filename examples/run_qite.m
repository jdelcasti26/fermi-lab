% Sample script to run QITE
clear;
tic;

noQubits = 4;
modelFam = "H2";

% Define pauliDictionary handle
PDH = pauliDictionary(noQubits);
partition = {[1:3], [4:15]};

% QITE kernel options
options.order = "2nd";
options.Pinv = false;

initialState = PDH.prepare_initial_state(modelFam);

% Define qiteKernel handle and ansätze according to partition
QK = qiteKernel(PDH, options, modelFam, 1, partition{:});
% Define qite object with initialState and scan parameters
Q = qite(QK, initialState, [0.05 : 0.05 : 0.50], 10);
% Run qite
[Q, summary, endState] = Q.estimate_energy();

totalTime = toc;
