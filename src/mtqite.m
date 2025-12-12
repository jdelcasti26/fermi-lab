% mtqite.m
classdef mtqite
    properties
        % Runtime configuration
        QK                 % Handle to qite_kernel 
        initialState       % Initial state for algorithm
        timeSeries         % Vector containing (positive) time step sizes to be scanned
        noTIMES            % Number of different time step sizes
        noTS               % Number of Trotter steps
        % Simulation results
        cellQiteOps        % Circuit definition for all times, along all Trotter steps, all terms in partition
        cellThetas         % Circuit parameters for the above
        measureCountLS     % Array with performed measurements to obtain linear system per time, per Trotter step per term
        measureCountEE     %% Array with performed measurements to obtain estimated energy per time, per Trotter step
    end % properties
    methods
        function o = mtqite(QK, initialState, timeSeries, noTS)
            % Inform handle to qite_kernel
            o.QK = QK;
            % Dimensonality check, statevector mast be 2^no_qubits x 1
            size_vector = size(initialState);
            if(size_vector(2) ~= 1 || size_vector(1) ~= o.QK.PDH.twoN)
                error("Inconsistent dimensions for statevector");
            else
                % Renormalization for safety
                o.initialState = initialState / norm(initialState);
            end       
            % Dimensonality check for timeSeries
            size_timeSeries = size(timeSeries);
            if(size_timeSeries(1) ~= 1)
                error("Expected a row vector for timeSeries");
            end
            if(sum(timeSeries <= 0))
                error("Expected positive values in timeSeries");
            end
            % Sort time step sizes for safety
            o.timeSeries = sort(unique(timeSeries));
            o.noTIMES = size_timeSeries(2);
            o.noTS = noTS;
            % Allocate container structures
            o.cellQiteOps = cell(o.noTIMES, o.noTS, o.QK.M);
            o.cellThetas = cell(o.noTIMES, o.noTS, o.QK.M);
            o.measureCountLS = zeros(o.noTIMES, o.noTS, o.QK.M);
            o.measureCountEE = zeros(o.noTIMES, o.noTS);
        end
        function o = scan_step(o, symmetryType, stateVector, trotterStep)
            % Initialize state in qite_kernel
            o.QK = o.QK.set_statevector(stateVector);
            o.QK.PDH = o.QK.PDH.flushMeasurements();
            if(symmetryType == "none")
                for mterm = 1 : o.QK.M % For each term in the partition
                    tn = 0; % Different time counter
                    for Dt = o.timeSeries % For every different time step size
                        tn = tn + 1;
                        % Measure only for first time step size regardless of
                        % trotterStep, term in partition mterm
                        if(Dt == o.timeSeries(1))
                            % Measure qite's linear system
                            o.QK = o.QK.measure_state(mterm);
                            % Inform performed number of distinct measurements
                            o.measureCountLS(tn, trotterStep, mterm) = o.QK.measureCount(mterm);
                            %fprintf("Paso %d\n", mterm);
                        end
                        % Obtain RTE configuration
                        [thetas, qiteOps] = o.QK.solve_linear_system(mterm, Dt);
                        % Inform circuit parameters in container
                        o.cellQiteOps{tn, trotterStep, mterm} = qiteOps;
                        o.cellThetas{tn, trotterStep, mterm} = thetas;
                    end % for Dt
                end % for mterm
            end % if(symmetryType == "none")
            if(symmetryType == "3rd_eq_1st")
                for mterm = 1 : 2 % For each independent term (3 will be copied from 1)
                    tn = 0; % Different time counter
                    for Dt = o.timeSeries % For every different time step size
                        tn = tn + 1;
                        % Measure only for first time step size regardless of
                        % trotterStep, term in partition mterm
                        if(Dt == o.timeSeries(1))
                            % Measure qite's linear system
                            o.QK = o.QK.measure_state(mterm);
                            % Inform performed number of distinct measurements
                            o.measureCountLS(tn, trotterStep, mterm) = o.QK.measureCount(mterm);
                            %fprintf("Paso %d\n", mterm);
                        end
                        % Obtain RTE configuration
                        [thetas, qiteOps] = o.QK.solve_linear_system(mterm, Dt);
                        % Inform circuit parameters in container
                        o.cellQiteOps{tn, trotterStep, mterm} = qiteOps;
                        o.cellThetas{tn, trotterStep, mterm} = thetas;
                    end % for Dt
                end % for mterm
                % For mterm = 3, copy from mterm = 1 with mirror symmetry
                % Inform circuit parameters in container
                tn = 0; % Different time counter
                for Dt = o.timeSeries % For every different time step size
                    tn = tn + 1;                
                    o.cellQiteOps{tn, trotterStep, 3} = o.mirror(o.cellQiteOps{tn, trotterStep, 1});
                    o.cellThetas{tn, trotterStep, 3} = o.cellThetas{tn, trotterStep, 1};
                end
            end % if(symmetryType == "3rd_eq_1st")            
        end
        % Mirror-symmetric array of Pauli strings
        function mirrorArrPS = mirror(o, arrPS)
            size_arrPS = size(arrPS, 2);
            mirrorArrPS = arrPS;
            for i = 1 : size_arrPS
                psLength = arrPS(i).len;
                for s = 1 : psLength
                    mirrorArrPS(i).pStr{s} = reverse(arrPS(i).pStr{s});
                end
            end
        end
        % Function estimate_energy_step
        % Calculates all the possible combinations of circuits U(3) * U(1) * U(2)
        % given a permutation like [3, 1, 2] of terms in the hamiltonian
        % Energy is estimated and minimized in order to find the optimum
        % configuration to evolve to the next statevector
        function [o, stuSummary, endState] = estimate_energy_step(o, perm, symmetryType, trotterStep, stateVector)
            % Check size of perm vector
            size_perm = size(perm);
            if(size_perm(1) ~= 1 || size_perm(2) ~= o.QK.M)
                error("Check the permutation parameter perm. Expected [p2, p1, ...]");
            end
            % Generate grid structure with combinations of time step sizes
            % gridTimesn : master grid, row vectors of length size_perm(2) with time step number tn 
            % gridTimesDt : row vectors of length size_perm(2) with time step values Dt
            % gridOpCount : operator count of the current circuit
            % gridPSCount : pauli string count of the current circuit
            % gridEnergy : energy estimate for current circuit
            % gridMeasurementsEE : measurement cost
            
            % Calculate grid of tn1, tn2, ...  coordinates
            % according to the specified permutation
            sets = cell(1, 0);
            for p = 1 : size_perm(2)
                sets{end + 1} = [1 : o.noTIMES];
            end
            gridTimesn = o.cartesianProduct(sets);
            % Filter redundant combinations by symmetry
            if(symmetryType == "3rd_eq_1st")
               if(sum(perm - [1, 2, 3]) == 0)
                   gridTimesn = gridTimesn(gridTimesn(:,1) == gridTimesn(:,3), :);
               else
                   error("Specified symmetry incompatible with specified partition");
               end
            end
            % Read size of grid (number of combinations)
            size_grid = size(gridTimesn, 1);
            % Allocate the rest of the grids
            gridTimesDt = zeros(size_grid, size_perm(2));
            gridOpCount = zeros(size_grid, 1);
            gridPSCount = zeros(size_grid, 1);
            gridEnergy = zeros(size_grid, 1);
            gridMeasurementsEE = zeros(size_grid, 1);
            % Calculate grid entries
            for id = 1 : size_grid
                circuit = eye(o.QK.PDH.twoN); 
                for p = 1 : size_perm(2)
                    gridTimesDt(id, p) = o.timeSeries(gridTimesn(id, p));
                    [circuitTerm, noOps, noPSs] = o.calculate_circuit(gridTimesn(id, p), trotterStep, p);
                    circuit = circuitTerm * circuit;
                    gridOpCount(id) = gridOpCount(id) + noOps;
                    gridPSCount(id) = gridPSCount(id) + noPSs;
                end
                % Energy calculation
                evolvedState = circuit * stateVector;
                o.QK.PDH = o.QK.PDH.flushMeasurements();
                measuresRef = o.QK.PDH.measuresCO;
                gridEnergy(id) = o.QK.HT.eValCO(evolvedState);
                gridMeasurementsEE(id) = o.QK.PDH.measuresCO - measuresRef;
            end         

            % Reference values
            exactEigval = o.QK.exact_eigenvalue();
            exactEigvecs = o.QK.exact_eigenvectors();
            
            [stuSummary.energy, atrow] = min(gridEnergy(:));
            stuSummary.opCount = gridOpCount(atrow);
            stuSummary.PSCount = gridPSCount(atrow);
            % Accuracy
            stuSummary.accuracy = (stuSummary.energy - exactEigval) / abs(exactEigval);
            % Calculate endstate and fidelity
            circuit = eye(o.QK.PDH.twoN); 
            for p = 1 : size_perm(2)
                gridTimesDt(atrow, p) = o.timeSeries(gridTimesn(atrow, p));
                [circuitTerm, ~, ~] = o.calculate_circuit(gridTimesn(atrow, p), trotterStep, p);
                circuit = circuitTerm * circuit;
            end
            endState = circuit * stateVector;
            stuSummary.fidelity = (endState' * exactEigvecs) * (endState' * exactEigvecs)';
            % Inform in which datapoint minimum has occured
            stuSummary.minPosition = gridTimesn(atrow,:);
            % Inform cost
            stuSummary.costLS = sum(sum(o.measureCountLS(:, trotterStep, :)));
            stuSummary.costEE = sum(gridMeasurementsEE);

            fprintf("MT-QITE step %d, energy: %f \n", trotterStep, stuSummary.energy);
            fprintf("Accuracy: %14.5e \n", stuSummary.accuracy);
            fprintf("Depth in implemented operators %d and strings: %d \n", stuSummary.opCount, stuSummary.PSCount);
            fprintf("Fidelity: %14.5e \n", stuSummary.fidelity);
            fprintf("Linear system %d, energy %d and total %d mesurements.\n", stuSummary.costLS , stuSummary.costEE, stuSummary.costLS + stuSummary.costEE);
            fprintf("Minimum found at datapoint(s) no.: %s \n", strjoin(string(stuSummary.minPosition), ' '));

        end
        function [o, stuSummary, endState] = estimate_energy(o, perm, symmetryType)
            if(symmetryType ~= "none" && symmetryType ~= "3rd_eq_1st")
                error("Requested symmetry not implemented");
            end
            stateVector = o.initialState;
            for ts = 1 : o.noTS
                o = o.scan_step(symmetryType, stateVector, ts);
                [o, arrStuSummary(ts), stateVector] = o.estimate_energy_step(perm, symmetryType, ts, stateVector);
            end
            stuSummary = table2struct(struct2table(arrStuSummary), 'ToScalar', true);
            % Cumulate depths and costs
            stuSummary.opCount = cumsum(stuSummary.opCount);
            stuSummary.PSCount = cumsum(stuSummary.PSCount);
            stuSummary.costLS = cumsum(stuSummary.costLS);
            stuSummary.costEE = cumsum(stuSummary.costEE);            
            endState = stateVector;
        end
        % Function calculate_circuit()
        function [Circuit, noOps, noPSs] = calculate_circuit(o, timeDatapoint, trotterStep, mterm)
            % Initialize with identity matrix
            Circuit = eye(o.QK.PDH.twoN);
            size_Ops = size(o.cellQiteOps{timeDatapoint, trotterStep, mterm}, 2);
            noOps = 0;
            noPSs = 0;
            for i = 1 : size_Ops
                % Calculate Pauli string in matrix form
                opMatrix = o.cellQiteOps{timeDatapoint, trotterStep, mterm}(i).matrixRep();
                theta = o.cellThetas{timeDatapoint, trotterStep, mterm}(i);
                Circuit = expm(- o.timeSeries(timeDatapoint) * theta * opMatrix) * Circuit;
                noOps = noOps + 1;
                noPSs = noPSs + o.cellQiteOps{timeDatapoint, trotterStep, mterm}(i).len;
            end
        end % calculate_circuit
        % Cartesian product of a variable number of sets
        % Sets has format {[1:2], [1:3], [1,5]}
        function result = cartesianProduct(o, sets)
            c = cell(1, numel(sets));
            [c{:}] = ndgrid(sets{:});
            result = cell2mat(cellfun(@(v)v(:), c, 'UniformOutput', false));
        end        
    end % methods
end % class definition