% Heatmap Hubbard model 6 qubits. Reported data
clear;
tic;

% Graph parameters
Field = 'Electrostatic repulsion U';
Depth = 'Trotter step';
Heat = 'Infidelity per measurement';

outputFile = "heatmapHubbard.mat";

% Model parameters
% Hubbard model
modelParams = [0.25, 0.50, 0.75, 1.0, 1.25, 1.50, 1.75, 2.00];

noParams = size(modelParams, 2);

% Model parameters
modelFam = "HUBOB06";
noQubits = 6;

% Time scan parameters
Dt = 5 / 120;
noTIMES = 12;
noTS = 10;
timeSeries = [Dt : Dt : noTIMES * Dt];

% QITE options
options.order = "2nd";
options.Pinv = true;

% MT-QITE options
symmetryType = "none";

% Partition and partition permutaion
partition = {[2:10], [11:18]}; % Hubbard
permutation = [1, 2];

% Define dictionary handle and initial state
PDH = pauliDictionary(noQubits);

initialState = PDH.prepare_initial_state(modelFam);

% Scan the different parameter values and return stuSummary
% arrays for both MT-QITE and QITE
for p = 1 : noParams
    modelParam = modelParams(p);

    % Define kernel handle
    QK = qiteKernel(PDH, options, modelFam, modelParam, partition{:});
    % Define mtqite object
    MTQ = mtqite(QK, initialState, timeSeries, noTS);
    % Define qite object
    Q = qite(QK, initialState, timeSeries, noTS);
    % Run mtqite
    [MTQ, arrStuSummaryMTQ(p), endState] = MTQ.estimate_energy(permutation, symmetryType);
    % Run qite
    [Q, arrStuSummaryQ(p), endState] = Q.estimate_energy();
end

% Prepare data for heatmap - MT-QITE
arrHeatmap = zeros(p * noTS, 3);
tabHeatmapMTQ = array2table(arrHeatmap, 'VariableNames', {Field, Depth, Heat});
for p = 1 : noParams
    % Inform parameter value on first column
    tabHeatmapMTQ.(Field)((p-1) * noTS + 1: (p-1) * noTS + noTS) = ones(noTS, 1) * modelParams(p);
    % Inform Trotter step on second column
    tabHeatmapMTQ.(Depth)((p-1) * noTS + 1: (p-1) * noTS + noTS) = (1 : noTS)';
    % Inform logarithmic infidelity per measurement on third column
    % Prepare measurements per step
    measuresPerStep = arrStuSummaryMTQ(p).costLS + arrStuSummaryMTQ(p).costEE;
    measuresPerStep = [measuresPerStep(1); diff(measuresPerStep)];
    tabHeatmapMTQ.(Heat)((p-1) * noTS + 1: (p-1) * noTS + noTS) = - log10(1 - arrStuSummaryMTQ(p).fidelity) ./ measuresPerStep;
end
% Prepare data for heatmap - QITE
tabHeatmapQ = array2table(arrHeatmap, 'VariableNames', {Field, Depth, Heat});
for p = 1 : noParams
    % Inform parameter value on first column
    tabHeatmapQ.(Field)((p-1) * noTS + 1: (p-1) * noTS + noTS) = ones(noTS, 1) * modelParams(p);
    % Inform Trotter step on second column
    tabHeatmapQ.(Depth)((p-1) * noTS + 1: (p-1) * noTS + noTS) = (1 : noTS)';
    % Inform logarithmic infidelity per measurement on third column
    % Prepare measurements per step
    measuresPerStep = arrStuSummaryQ(p).costLS + arrStuSummaryQ(p).costEE;
    measuresPerStep = [measuresPerStep(1); diff(measuresPerStep)];
    tabHeatmapQ.(Heat)((p-1) * noTS + 1: (p-1) * noTS + noTS) = - log10(1 - arrStuSummaryQ(p).fidelity) ./ measuresPerStep;
end

% Heatmap elaboration
% Colormap
ncs = 256;               % 256 colors
gamma = 0.3;             % < 1 expands the low end; > 1 expands the high end

cm  = parula(ncs);
x   = linspace(0, 1, ncs).'; % original positions along the map
xw  = x .^ gamma;          % warped positions (gamma<1 -> richer lows)
cmw = [interp1(x, cm(:,1), xw, 'linear') ...
       interp1(x, cm(:,2), xw, 'linear') ...
       interp1(x, cm(:,3), xw, 'linear')];

clear h;
fig = figure();
tcl = tiledlayout(fig, 1, 2);
n = 2;  % number of heatmaps
h = gobjects(n, 1); 
ax = nexttile(tcl); 
h(1) = heatmap(tabHeatmapQ, Field, Depth, 'ColorVariable', Heat, ...
               'ColorbarVisible', 'off', 'GridVisible', 'off', 'CellLabelColor', 'none');
h(1).YDisplayData = flipud(h(1).YDisplayData);

colormap(h(1), cmw);
h(1).Title = 'QITE';
h(1).XLabel = Field;
h(1).YLabel = Depth;
h(1).FontName = 'Times';
h(1).FontSize = 20;  

ax = nexttile(tcl); 
h(2) = heatmap(tabHeatmapMTQ, Field, Depth, 'ColorVariable', Heat, ...
               'ColorbarVisible', 'off', 'GridVisible', 'off', 'CellLabelColor', 'none');
h(2).YDisplayData = flipud(h(2).YDisplayData);

colormap(h(2), cmw);
h(2).Title = 'MT-QITE';
h(2).XLabel = Field;
h(2).YLabel = Depth;
h(2).FontName = 'Times';
h(2).FontSize = 20;  


% Equate color limits in all heatmaps
colorLims = vertcat(h.ColorLimits);
globalColorLim = [min(colorLims(:,1)), max(colorLims(:,2))];
set(h, 'ColorLimits', globalColorLim)
% Create global colorbar that uses the global color limits
ax = axes(tcl, 'visible', 'off', 'Colormap', h(1).Colormap, 'CLim', globalColorLim);
cb = colorbar(ax);
cb.Layout.Tile = 'East';
cb.FontSize = 16;

totalTime = toc;
save(outputFile);