% Sample script to run QITE
clear;
tic;

%noQubits = 4;
%modelFam = "ISIOB04";

modelFam = "HEYOB06";
noQubits = 6;

% Define pauliDictionary handle
PDH = pauliDictionary(noQubits);
%partition = {[1,4,5,7], [2,3,6]}; % ISI06
%partition = {[1,2,5,6], [3,4,7]}; % ISI06
partition2P = {[1,2,3,7,8,9,13,14,15], [4,5,6,10,11,12,16,17,18]}; % Heisenberg (symmetric)
partition1P = {[1:18]};

% QITE kernel options
options.order = "2nd";
options.Pinv = true;

% MT-QITE options
symmetryType = "none";

Dt = 5/120;
noTIMES = 12;
noTS = 10;
timeSeries = [Dt : Dt : noTIMES * Dt];

initialState = PDH.prepare_initial_state(modelFam);

% Define qiteKernel handle and ansätze according to partition
QK = qiteKernel(PDH, options, modelFam, 1, partition1P{:});

% Define qite object with initialState and scan parameters
Q = qite(QK, initialState, timeSeries, noTS);
% Run qite
[Q, stuSummaryQ1P, endState] = Q.estimate_energy();

% Define sm-qite object with initialState and scan parameters
SMQ = smqite(QK, initialState, timeSeries, noTS);
% Run sm-qite
[SMQ, stuSummarySMQ1P, endState] = SMQ.estimate_energy();

% Define ac-qite object with initialState and scan parameters
ACQ = acqite(QK, initialState, timeSeries, noTS);
% Run ac-qite
[ACQ, stuSummaryACQ1P, endState] = ACQ.estimate_energy();

% Define mtqite object with initialState and scan parameters
MTQ = mtqite(QK, initialState, timeSeries, noTS);
% Run MT-QITE
[MTQ, stuSummaryMTQ1P, endState] = MTQ.estimate_energy([1], symmetryType);

clear QK;
clear Q;
clear MTQ;
clear ACQ;
clear SMQ;

% Define qiteKernel handle and ansätze according to partition
QK = qiteKernel(PDH, options, modelFam, 1, partition2P{:});

% Define qite object with initialState and scan parameters
Q = qite(QK, initialState, timeSeries, noTS);
% Run qite
[Q, stuSummaryQ2P, endState] = Q.estimate_energy();

% Define sm-qite object with initialState and scan parameters
SMQ = smqite(QK, initialState, timeSeries, noTS);
% Run sm-qite
[SMQ, stuSummarySMQ2P, endState] = SMQ.estimate_energy();

% Define ac-qite object with initialState and scan parameters
ACQ = acqite(QK, initialState, timeSeries, noTS);
% Run ac-qite
[ACQ, stuSummaryACQ2P, endState] = ACQ.estimate_energy();

% Define mtqite object with initialState and scan parameters
MTQ = mtqite(QK, initialState, timeSeries, noTS);
% Run MT-QITE
[MTQ, stuSummaryMTQ2P, endState] = MTQ.estimate_energy([1, 2], symmetryType);

clear QK;
clear Q;
clear MTQ;
clear ACQ;
clear SMQ;

qiteGraphs3P2Pb(stuSummaryMTQ2P, stuSummaryQ2P, stuSummaryACQ2P, stuSummarySMQ2P, "trotterSteps", "fidelity");
qiteGraphs3P2Pb(stuSummaryMTQ2P, stuSummaryQ2P, stuSummaryACQ2P, stuSummarySMQ2P, "PSCount", "fidelity");
qiteGraphs3P2Pb(stuSummaryMTQ2P, stuSummaryQ2P, stuSummaryACQ2P, stuSummarySMQ2P, "costLS", "fidelity");

qiteGraphs3P2Pb(stuSummaryMTQ1P, stuSummaryQ1P, stuSummaryACQ1P, stuSummarySMQ1P, "trotterSteps", "fidelity");
qiteGraphs3P2Pb(stuSummaryMTQ1P, stuSummaryQ1P, stuSummaryACQ1P, stuSummarySMQ1P, "PSCount", "fidelity");
qiteGraphs3P2Pb(stuSummaryMTQ1P, stuSummaryQ1P, stuSummaryACQ1P, stuSummarySMQ1P, "costLS", "fidelity");

totalTime = toc;
save("similarWorkHei.mat");
