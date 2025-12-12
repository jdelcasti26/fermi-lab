classdef smqite
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
        measureCountEE     % Array with performed measurements to obtain estimated energy per time, per Trotter step
        arrayE             % Energy array (per time datapoint, per Trotter step)
        arrOpCount         % Operator (in the ansatz) count, per time datapoint, per Trotter step 
        arrPSCount         % Pauli string count, per time datapoint, per Trotter step 
    end % properties
    methods
        function o = smqite(QK, initialState, timeSeries, noTS)
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
        function o = scan(o)
            tn = 0; % Different time counter
            for Dt = o.timeSeries % For every different time step size
                % Initialize state in qite_kernel
                o.QK = o.QK.set_statevector(o.initialState);
                o.QK.PDH = o.QK.PDH.flushMeasurements();
                tn = tn + 1;
                for ts = 1 : o.noTS % For each Trotter step
                    for mterm = 1 : o.QK.M % For each term in the partition
                        if(ts == 1 && mterm == 1)
                            if(Dt == o.timeSeries(1))
                                % Measure qite's linear system
                                o.QK = o.QK.measure_state(mterm);
                                % Copy results to be reused later
                                o.QK.reuse = o.QK.measured;
                                % Inform performed number of distinct measurements
                                o.measureCountLS(tn, ts, mterm) = o.QK.measureCount(mterm);
                            end
                            o.QK.measured = o.QK.reuse;
                            % Obtain RTE configuration
                            [thetas, qiteOps] = o.QK.solve_linear_system(mterm, Dt);
                            % Inform circuit parameters in container
                            % Container accumulates Qite's operations for
                            % all Trotter steps, all terms in the partition
                            if(ts == 1 && mterm == 1) % Always the case
                                o.cellQiteOps{tn, ts, mterm} = qiteOps;
                                o.cellThetas{tn, ts, mterm} = thetas;
                            end
                            % Evolve state with accumulated Qite ops from initialState
                            evolvedState = o.calculate_circuit(tn, ts, mterm) * o.initialState;
                            % Updatet state implies stored measurement values are no longer valid
                            o.QK = o.QK.set_statevector(evolvedState);
                            o.QK.PDH = o.QK.PDH.flushMeasurements();
                        else
                            % Measure qite's linear system
                            o.QK = o.QK.measure_state(mterm);
                            % Inform performed number of distinct measurements
                            o.measureCountLS(tn, ts, mterm) = o.QK.measureCount(mterm);
                            % Obtain RTE configuration
                            [thetas, qiteOps] = o.QK.solve_linear_system(mterm, Dt);
                            % Inform circuit parameters in container
                            % Container accumulates Qite's operations for
                            % all Trotter steps, all terms in the partition
                            if(ts == 1 && mterm == 1)
                                % Code never reached
                                o.cellQiteOps{tn, ts, mterm} = qiteOps;
                                o.cellThetas{tn, ts, mterm} = thetas;
                            elseif(mterm > 1)
                                [o.cellQiteOps{tn, ts, mterm}, o.cellThetas{tn, ts, mterm}] = o.accumulate_qite_ops(o.cellQiteOps{tn, ts, mterm - 1}, o.cellThetas{tn, ts, mterm - 1}, qiteOps, thetas);
                            elseif(mterm == 1 && ts > 1)
                                [o.cellQiteOps{tn, ts, mterm}, o.cellThetas{tn, ts, mterm}] = o.accumulate_qite_ops(o.cellQiteOps{tn, ts - 1, mterm}, o.cellThetas{tn, ts - 1, mterm}, qiteOps, thetas);
                            else
                                error("Internal indexing error");
                            end
                            % Evolve state with accumulated Qite ops from initialState
                            evolvedState = o.calculate_circuit(tn, ts, mterm) * o.initialState;
                            % Updatet state implies stored measurement values are no longer valid
                            o.QK = o.QK.set_statevector(evolvedState);
                            o.QK.PDH = o.QK.PDH.flushMeasurements();                            
                        end % if ts 1, mterm 1
                    end % for mterm
                end % for ts
            end % for t
        end
        function [o, stuSummary, finalState] = estimate_energy(o)
            o.arrayE = zeros(o.noTIMES, o.noTS);
            o.arrOpCount = zeros(o.noTIMES, o.noTS);
            o.arrPSCount = zeros(o.noTIMES, o.noTS);
            % Scan time, partition terms, Trotter steps
            o = o.scan();
            for tn = 1 : o.noTIMES
                for ts = 1 : o.noTS
                    % Energy calculation
                    state = o.calculate_circuit(tn, ts, o.QK.M) * o.initialState;
                    o.QK.PDH = o.QK.PDH.flushMeasurements();
                    measures_ref = o.QK.PDH.measuresCO;
                    o.arrayE(tn, ts) = o.QK.HT.eValCO(state);
                    o.measureCountEE(tn, ts) = o.QK.PDH.measuresCO - measures_ref;
                    % Circuit depth count    
                    size_Ops = size(o.cellQiteOps{tn, ts, o.QK.M}, 2);
                    o.arrOpCount(tn, ts) = size_Ops;
                    for i = 1 : size_Ops
                        o.arrPSCount(tn, ts) = o.arrPSCount(tn, ts) + o.cellQiteOps{tn, ts, o.QK.M}(i).len;
                    end
                end
            end
            % Summary data to be returned
            stuSummary.energy = zeros(o.noTS, 1);
            stuSummary.accuracy = zeros(o.noTS, 1);
            stuSummary.fidelity = zeros(o.noTS, 1);
            stuSummary.opCount = zeros(o.noTS, 1);
            stuSummary.PSCount = zeros(o.noTS, 1);
            stuSummary.minPosition = zeros(o.noTS, 1);
            stuSummary.costLS = zeros(o.noTS, 1);
            stuSummary.costEE = zeros(o.noTS, 1);
            % Reference values
            exactEigenvalue = o.QK.exact_eigenvalue();
            exactEigenvectors = o.QK.exact_eigenvectors();
            for ts = 1 : o.noTS
                [stuSummary.energy(ts), atrow] = min(o.arrayE(:, ts));
                stuSummary.opCount(ts) = o.arrOpCount(atrow, ts);
                stuSummary.PSCount(ts) = o.arrPSCount(atrow, ts);
                % Accuracy
                stuSummary.accuracy(ts) = (stuSummary.energy(ts) - exactEigenvalue) / abs(exactEigenvalue);
                % Calculate endstate and fidelity
                endState = o.calculate_circuit(atrow, ts, o.QK.M) * o.initialState;
                stuSummary.fidelity(ts) = (endState' * exactEigenvectors) * (endState' * exactEigenvectors)';
                % Inform in which datapoint minimum has occured
                stuSummary.minPosition(ts) = atrow;
                % Inform cost
                stuSummary.costLS(ts) = sum(sum(o.measureCountLS(:, ts, :)));
                stuSummary.costEE(ts) = sum(o.measureCountEE(:, ts));
            end
            % Cumulate cost
            stuSummary.costLS = cumsum(stuSummary.costLS);
            stuSummary.costEE = cumsum(stuSummary.costEE);

            [minE, atstep] = min(stuSummary.energy);
            % Final state calculation
            finalState = o.calculate_circuit(stuSummary.minPosition(atstep), atstep, o.QK.M) * o.initialState;
            fprintf("SM-QITE ground energy: %f \n", minE);
            fprintf("Accuracy: %14.5e \n", stuSummary.accuracy(atstep));
            fprintf("Fidelity: %14.5e \n", stuSummary.fidelity(atstep));
            fprintf("Depth in implemented operators %d and strings: %d \n", stuSummary.opCount(atstep), stuSummary.PSCount(atstep));
            fprintf("Effective number of Trotter steps: %d \n", atstep);
            fprintf("Linear system %d, energy %d and total %d mesurements.\n", stuSummary.costLS(atstep) , stuSummary.costEE(atstep), stuSummary.costLS(atstep) + stuSummary.costEE(atstep));
            fprintf("Minimum found at datapoint no.: %d \n", stuSummary.minPosition(atstep));

            % Exact ITE for reference in graph
            ExactITE = zeros(o.noTIMES, 1);
            % Hamiltonian matrix
            HM = o.QK.HT.matrixRep();
            for tm = 1 : o.noTIMES
                Dtm = o.timeSeries(tm);
                ExactITE(tm) = real(o.initialState' * expm(-Dtm * HM) * HM * expm(-Dtm * HM) * o.initialState) / real(o.initialState' * expm(-Dtm * HM) * expm(-Dtm * HM) * o.initialState);
            end
            % Graph
            Beta = o.timeSeries';
            %szz = size(o.Earray(:, 1 : o.effnts));
            plot(Beta, o.arrayE(:, :), Beta, ExactITE, '+');
            %plot(o.effnts * Beta, o.Earray(:, o.effnts), Beta, o.Eexactite, '+', o.t0 + o.Ddt * o.dtptbeta, o.Emin, 'ro');
            %leg = string([1 : o.effnts]);
            leg = sprintf("Trotter step %d", 1);
            for ts = 2 : o.noTS
                leg(end + 1) = sprintf("Trotter step %d", ts);
            end
            leg(end + 1) = sprintf("Exact ITE");
            
            tl = title('Classical and Quantum ITE');
            xl = xlabel('\Delta\tau');
            yl = ylabel('Energy');
            lgd = legend(leg);
            tl.FontSize = 14;
            xl.FontSize = 14;
            yl.FontSize = 14;
            lgd.FontSize = 14;
        end % estimateEnergy
        function Circuit = calculate_circuit(o, timeDatapoint, trotterStep, mterm)
            % Initialize with identity matrix
            Circuit = eye(o.QK.PDH.twoN);
            size_Ops = size(o.cellQiteOps{timeDatapoint, trotterStep, mterm}, 2);
            for i = 1 : size_Ops
                % Calculate Pauli string in matrix form
                opMatrix = o.cellQiteOps{timeDatapoint, trotterStep, mterm}(i).matrixRep();
                theta = o.cellThetas{timeDatapoint, trotterStep, mterm}(i);
                Circuit = expm(- o.timeSeries(timeDatapoint) * theta * opMatrix) * Circuit;
            end
        end % calculate_circuit
        % Function to accumulate Qite's operators and parameters for
        % step-merged type algorithms. Parameters 'thetas' are accumulated
        % according to pauliStrings in array 'qiteOps'. Arrays are expected
        % to be row arrays
        function [accQiteOps, accThetas] = accumulate_qite_ops(o, baseQiteOps, baseThetas, qiteOps, thetas)
            % Dimensions
            size_base = size(baseQiteOps, 2);
            size_new = size(qiteOps, 2);
            % Allocate results
            accQiteOps = baseQiteOps;
            accThetas = baseThetas;
            % For every pauliString in new qiteOps
            % loopaccQiteOps to check if new pauliString is already present, if so, accumulate parameters thetas 
            for new = 1 : size_new
                flagExists = false;
                for base = 1 : size_base
                    % Check if pauiStrings are equal
                    if(isemptyorzero(baseQiteOps(base) - qiteOps(new)))
                        accThetas(base) = accThetas(base) + thetas(new);
                        flagExists= true;
                        break;
                    end
                end
                if(~flagExists)
                    % Append operator
                    accQiteOps = [accQiteOps, qiteOps(new)];
                    accThetas = [accThetas, thetas(new)];
                end
            end
        end
    end % methods
end % class definition