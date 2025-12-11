% rep*.m Compare MT-QITE vs QITE Heisenberg model on 8 qubits. Reported data
clear;
tic;

% Model parameters
modelFam = "HEYOB08";
batchFile = "initials/iniStatesHEYOB08.mat";
outputFile = "HEY08OB.mat";
modelParam = 1;
noQubits = 8;

% Time scan parameters
Dt = 5 / 120;
noTIMES = 12;
noTS = 10;
timeSeries = [Dt : Dt : noTIMES * Dt];

% QITE kernel options
options.order = "2nd";
options.Pinv = false;

% MT-QITE options
symmetryType = "none";

% Partition and partition permutaion
partition = {[1,2,3,4,9,10,11,12,17,18,19,20], [5,6,7,8,13,14,15,16,21,22,23,24]};
permutation = [1, 2];

% Define dictionary handle and initial state
PDH = pauliDictionary(noQubits);

load(batchFile);
noSamples = size(arrInitialStates, 2);


% Define kernel handle (common for every initialState given the partition)
QK = qiteKernel(PDH, options, modelFam, modelParam, partition{:});


for s = 1 : noSamples
    initialState = arrInitialStates(:, s);

    % Define mtqite object
    MTQ = mtqite(QK, initialState, timeSeries, noTS);
    % Define qite object
    Q = qite(QK, initialState, timeSeries, noTS);
    % Run mtqite
    [MTQ, arrStuSummaryMTQ(s), endState] = MTQ.estimate_energy(permutation, symmetryType);
    % Run qite
    [Q, arrStuSummaryQ(s), endState] = Q.estimate_energy();
end

% Plot graphs for arrays of summary structs 
% (if more than one sample plot error bars)
qiteGraphs(arrStuSummaryMTQ, arrStuSummaryQ, "PSCount", "fidelity");
qiteGraphs(arrStuSummaryMTQ, arrStuSummaryQ, "costTotal", "fidelity");

totalTime = toc;

save(outputFile, 'QK', 'MTQ', 'Q', 'arrStuSummaryMTQ', 'arrStuSummaryQ', 'totalTime', 'noSamples', 'noTS');

