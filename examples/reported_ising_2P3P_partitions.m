% rep*.m Compare MT-QITE vs QITE with 2 and 3-term partitions. Reported data
clear;
tic;

% Model parameters
modelFam = "ISIOB06";
modelParam = 4;
noQubits = 6;
outputFile = "ISIOB06.mat";

% Time scan parameters
Dt = 5 / 120;
noTIMES = 12;
noTS = 10;
timeSeries = [Dt : Dt : noTIMES * Dt];

% QITE kernel options
options.order = "2nd";
options.Pinv = true;

% Define dictionary handle and initial state
PDH = pauliDictionary(noQubits);

initialState = PDH.prepare_initial_state(modelFam);

% Partition and partition permutaion
partition2P = {[1,2,7,8, 5,6,10,11], [3,4,9]}; % ISI06 2P
permutation2P = [1, 2];

partition3P = {[1,2,7,8], [3,4,9], [5,6,10,11]}; % ISI06 3P
permutation3P = [1, 2, 3];

%%%%%%% 2P CALCULATIONS
% MT-QITE options
symmetryType = "none";

% Define kernel handle (common for every initialState given the partition)
QK = qiteKernel(PDH, options, modelFam, modelParam, partition2P{:});

% Define mtqite object
MTQ = mtqite(QK, initialState, timeSeries, noTS);
% Define qite object
Q = qite(QK, initialState, timeSeries, noTS);
% Run mtqite
[MTQ, stuSummaryMTQ_2P, endState] = MTQ.estimate_energy(permutation2P, symmetryType);
% Run qite
[Q, stuSummaryQ_2P, endState] = Q.estimate_energy();


%%%%%%% 3P CALCULATIONS
% MT-QITE options
symmetryType = "3rd_eq_1st";

% Define kernel handle (common for every initialState given the partition)
QK = qiteKernel(PDH, options, modelFam, modelParam, partition3P{:});

% Define mtqite object
MTQ = mtqite(QK, initialState, timeSeries, noTS);
% Define qite object
Q = qite(QK, initialState, timeSeries, noTS);
% Run mtqite
[MTQ, stuSummaryMTQ_3P, endState] = MTQ.estimate_energy(permutation3P, symmetryType);
% Run qite
[Q, stuSummaryQ_3P, endState] = Q.estimate_energy();

% Plot graphs for arrays of summary structs 
qiteGraphs3P2P(stuSummaryMTQ_3P, stuSummaryQ_3P, stuSummaryMTQ_2P, stuSummaryQ_2P,"PSCount", "fidelity");
qiteGraphs3P2P(stuSummaryMTQ_3P, stuSummaryQ_3P, stuSummaryMTQ_2P, stuSummaryQ_2P,"costLS", "fidelity");

totalTime = toc;

save(outputFile, 'stuSummaryMTQ_2P', 'stuSummaryQ_2P', 'stuSummaryMTQ_3P', 'stuSummaryQ_3P', 'totalTime', 'noTS');

