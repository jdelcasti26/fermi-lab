% Sample script to get H-H-H-H binding curve. Reported data
clear;

noQubits = 8;
modelFam = "H4CHAIN";
partition1P = {[1:110]};
permutation1P = [1];
% Partitions are selected with bespoke code whose results are hardcoded
% here. The procedure is the following: calculate the commutator structure
% of ansatz operator terms agains hamiltonian terms in matrix form 
% (1 commute / 0 do not commute). Reorder rows and colums to make the 
% matrix as block diagonal as possible (e.g. using k-means).
% Partition according to lists of reordered ansatz terms. 
partition2P = {[3,4,7,8,10,12,13,14,17,26,29,35,38,41,45,47,50,54,58,59,60,61,66,68,69,71,74,...
               75,76,77,79,85,87,90,92,95,97,98,101,102,103,106,109],...
              [1,2,5,6,9,11,15,16,18,19,20,21,22,23,24,25,27,28,30,31,32,33,34,36,37,39,40,42,43,44,46,48,...
               49,51,52,53,55,56,57,62,63,64,65,67,70,72,73,78,80,81,82,83,84,86,88,89,91,93,94,96,99,100,104,105,107,108,110]};
permutation2P = [1, 2];
partition3P = {[16,17,18,20,27,28,29,32,36,37,38,39,41,46,47,48,50,51,55,56,63,64,65,67,70,71,72,73,81,85,86,87,90,91,92,96,100,104,107],...
               [1,2,5,6,9,11,15,19,21,22,23,24,25,30,31,33,34,40,42,43,44,49,52,53,57,62,78,80,82,83,84,88,89,93,94,99,105,108,110],...
               [3,4,7,8,10,12,13,14,26,35,45,54,58,59,60,61,66,68,69,74,75,76,77,79,95,97,98,101,102,103,106,109]};
permutation3P = [1, 2, 3];

% Time scan parameters
Dt = 5/120;
noTIMES = 12;
noTS = 10;
timeSeries = [Dt : Dt : noTIMES * Dt];

% QITE kernel options
options.order = "1st";
options.Pinv = false;
% MT-QITE options
symmetryType = "none";


% Interatomic distance values
noPoints = 13;
% Model names and filenames modelFam + fileSuffix
fileSuffix = ["060", "065", "070", "075", "080", "085", "090", "095", "100", "105", "110", "115", "120"];

PDH = pauliDictionary(noQubits);
initialState = PDH.prepare_initial_state(modelFam + fileSuffix(1)); % All models start with the HF state

for i = 1 : noPoints
%for i = {{ALPHA}} : {{ALPHA}}
    tic;
    %%% TRIVIAL PARTITION QITE + MT-QITE
    % Define qiteKernel handle and ansätze for trivial partition
    QK = qiteKernel(PDH, options, modelFam + fileSuffix(i), 1, partition1P{:});
    % Define qite object with initialState and scan parameters
    Q = qite(QK, initialState, timeSeries, noTS);
    % Run qite
    [Q, stuSummaryQ1P, endState] = Q.estimate_energy();
    % Define mtqite object with initialState and scan parameters
    MTQ = mtqite(QK, initialState, timeSeries, noTS);
    % Run MT-QITE
    [MTQ, stuSummaryMTQ1P, endState] = MTQ.estimate_energy(permutation1P, symmetryType);
    clear QK;
    clear Q;
    clear MTQ;

    %XX TWO-TERM PARTITION MT-QITE
    % Define qiteKernel handle and ansätze for two-term partition
    QK = qiteKernel(PDH, options, modelFam + fileSuffix(i), 1, partition2P{:});
    % Define mtqite object with initialState and scan parameters
    MTQ = mtqite(QK, initialState, timeSeries, noTS);
    % Run MT-QITE
    [MTQ, stuSummaryMTQ2P, endState] = MTQ.estimate_energy(permutation2P, symmetryType);
    clear QK;
    clear MTQ;

    %XX THREE-TERM PARTITION MT-QITE
    % Define qiteKernel handle and ansätze for three-term partition
    QK = qiteKernel(PDH, options, modelFam + fileSuffix(i), 1, partition3P{:});
    % Define mtqite object with initialState and scan parameters
    MTQ = mtqite(QK, initialState, timeSeries, noTS);
    % Run MT-QITE
    [MTQ, stuSummaryMTQ3P, endState] = MTQ.estimate_energy(permutation3P, symmetryType);
    clear QK;
    clear MTQ;

    fileName = modelFam + fileSuffix(i) + ".mat";
    totalTime = toc;
    save(fileName, 'stuSummaryQ1P', 'stuSummaryMTQ1P', 'stuSummaryMTQ2P', 'stuSummaryMTQ3P', 'totalTime');
end



