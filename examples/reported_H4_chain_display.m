% Sample script to get H-H-H-H binding curve. Reported data
clear;
tic;

flagcore = 1; % 1 to incluide nuclear repulsion, 0 otherwise
trs = 6; % Number of effective Trotter steps included in displayed data

noQubits = 8;
modelFam = "H4CHAIN";

noTS = 10; % Total number of Trotter steps available

partition1P = {[1:110]};


% Interatomic distance values
noPoints = 13;
distance = [0.60 : 0.05 : 1.20]; % Available interatomic distances
% Model names and filenames modelFam + fileSuffix
fileSuffix = ["060", "065", "070", "075", "080", "085", "090", "095", "100", "105", "110", "115", "120"];
% Nuclear repulsion from PySCF
arrEcore = [3.821835412200, 3.527848072800, 3.275858924743, 3.057468329760, 2.866376559150, 2.697766173318, 2.547890274800, 2.413790786653, 2.293101247320, 2.183905949829, 2.084637497564, 1.994001084626, 1.910917706100];

PDH = pauliDictionary(noQubits);
initialState = PDH.prepare_initial_state(modelFam + distance(1)); % All models start with the HF state

summary = [];

for i = 1 : noPoints
    d = distance(i);
    % Hartree-Fock reference energy check
    H = hamiltonian(PDH, modelFam + fileSuffix(i), 1, partition1P{:});
    Ehf = initialState' * H.HT.matrixRep() * initialState;
    Efci = H.exact_eigenvalue();
    Ecore = arrEcore(i);
    Ehf = Ehf + flagcore * Ecore;
    Efci = Efci + flagcore * Ecore;
    load("data\" + modelFam + fileSuffix(i) + ".mat");
    Eq1P = stuSummaryQ1P.energy(trs) + flagcore * Ecore;
    Emtq1P = stuSummaryMTQ1P.energy(trs) + flagcore * Ecore;
    Emtq2P = stuSummaryMTQ2P.energy(trs) + flagcore * Ecore;
    Emtq3P = stuSummaryMTQ3P.energy(trs) + flagcore * Ecore;
    summary = [summary; d, Ehf, Efci, Emtq1P, Emtq2P, Emtq3P, Eq1P];

    fprintf("d = %f, HF = %f, FCI = %f\n", d, Ehf + flagcore * Ecore, Efci + flagcore * Ecore);

    clear stuSummaryQ1P;
    clear stuSummaryMTQ1P;
    clear stuSummaryMTQ2P;
    clear stuSummaryMTQ3P;
    clear totalTime;
end

% Graphs

palette = [ ...
    0.000 0.000 0.000;  % Black
    0.902 0.624 0.000;  % Orange        (230,159,0)
    0.337 0.706 0.914;  % Sky blue      (86,180,233)
    0.000 0.620 0.451;  % Bluish green  (0,158,115)
    0.941 0.894 0.259;  % Yellow        (240,228,66)
    0.000 0.447 0.698;  % Blue          (0,114,178)
    0.835 0.369 0.000;  % Vermilion     (213,94,0)
    0.800 0.475 0.655];

noPoints = size(summary, 1);

figure
hold on
hh1 = plot(summary(:,1), ones(noPoints, 1) * 0.00159, '--', 'Color', palette(6,:)); % Definition of chemical accuracy 1 Kcal/mol = 1.59e-3 Hartree/molecule
hh4 = plot(summary(:,1), summary(:,4) - summary(:,3), '*', 'Color', palette(2,:));
hh5 = plot(summary(:,1), summary(:,5) - summary(:,3), 'o', 'Color', palette(4,:));
hh6 = plot(summary(:,1), summary(:,6) - summary(:,3), '*', 'Color', palette(4,:));
hh7 = plot(summary(:,1), summary(:,7) - summary(:,3), 'o', 'Color', palette(2,:));
hold off

xl = xlabel(strcat("H-H distance (", char(197), ')'));
yl = ylabel('Error \DeltaE (hartree / mol)');
lgd = legend([hh1, hh7, hh4, hh5, hh6], {'Chemical acuracy', 'QITE', 'MT-QITE 1P', 'MT-QITE 2P', 'MT-QITE 3P'});
xl.FontSize = 16;
yl.FontSize = 16;
lgd.FontSize = 16;
ax.FontSize = 18;

figure
hold on
h1 = plot(summary(:,1), summary(:,2), '--', 'Color', palette(6,:));
h2 = plot(summary(:,1), summary(:,3), '-', 'Color', palette(6,:));
h3 = plot(summary(:,1), summary(:,4), '*', 'Color', palette(2,:));
h4 = plot(summary(:,1), summary(:,5), 'o', 'Color', palette(4,:));
h5 = plot(summary(:,1), summary(:,6), '*', 'Color', palette(4,:));
h6 = plot(summary(:,1), summary(:,7), 'o', 'Color', palette(2,:));
hold off

xl = xlabel(strcat("H-H distance (", char(197), ')'));
yl = ylabel('Energy (hartree)');
lgd = legend([h1, h2, h6, h3, h4, h5], {'Hartree-Fock', 'Exact', 'QITE', 'MT-QITE 1P', 'MT-QITE 2P', 'MT-QITE 3P'});
xl.FontSize = 16;
yl.FontSize = 16;
lgd.FontSize = 16;
ax.FontSize = 18;
