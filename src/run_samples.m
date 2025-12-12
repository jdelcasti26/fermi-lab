% Compare MT-QITE vs QITE. Scan batch of initial state samples
clear;
tic;

% Model parameters
modelFam = "ISIOB04";
batchFile = "initials/iniStatesISIOB04.mat";
outputFile = "ISIOB04.mat";
modelParam = 1;
noQubits = 4;

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
partition = {[1,2,5,6], [3,4,7]};
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
    [Q, arrStuSummaryQ(s)] = Q.estimate_energy();
end

% Plot graphs for arrays of summary structs 
% (if more than one sample plot error bars)
qiteGraphs(arrStuSummaryMTQ, arrStuSummaryQ, "PSCount", "fidelity");
qiteGraphs(arrStuSummaryMTQ(1), arrStuSummaryQ(1), "costTotal", "fidelity");

totalTime = toc;

save(outputFile, 'QK', 'MTQ', 'Q', 'arrStuSummaryMTQ', 'arrStuSummaryQ', 'totalTime', 'noSamples', 'noTS');

