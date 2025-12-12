% Sample script to run MT-QITE
clear;
tic;

noQubits = 4;
modelFam = "H2";
% Define pauliDictionary handle
PDH = pauliDictionary(noQubits);
part = {[1:3], [4:15]};

% QITE kernel options
options.order = "2nd";
options.Pinv = false;

% MT-QITE options
symmetryType = "none";

initialState = PDH.prepare_initial_state(modelFam);

% Define qiteKernel handle and ansätze according to partition
QK = qiteKernel(PDH, options, modelFam, 1, part{:});
% Define mtqite object with initialState and scan parameters
MTQ = mtqite(QK, initialState, [0.05 : 0.05 : 0.50], 2);
% Run MT-QITE
[MTQ, stuSummary, endState] = MTQ.estimate_energy([1, 2], symmetryType);

totalTime = toc;
