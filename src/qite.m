classdef qite
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
        function o = qite(QK, initialState, timeSeries, noTS)
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
                            o.cellQiteOps{tn, ts, mterm} = qiteOps;
                            o.cellThetas{tn, ts, mterm} = thetas;
                            % Evolve state
                            evolvedState = o.QK.evolve_state(Dt, thetas, qiteOps);
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
                            o.cellQiteOps{tn, ts, mterm} = qiteOps;
                            o.cellThetas{tn, ts, mterm} = thetas;
                            % Evolve state
                            evolvedState = o.QK.evolve_state(Dt, thetas, qiteOps);
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
                    state = o.calculate_circuit(tn, ts) * o.initialState;
                    o.QK.PDH = o.QK.PDH.flushMeasurements();
                    measures_ref = o.QK.PDH.measuresCO;
                    o.arrayE(tn, ts) = o.QK.HT.eValCO(state);
                    o.measureCountEE(tn, ts) = o.QK.PDH.measuresCO - measures_ref;
                    % Circuit depth count    
                    for mterm = 1 : o.QK.M
                        size_Ops = size(o.cellQiteOps{tn, ts, mterm}, 2);
                        for i = 1 : size_Ops                    
                            o.arrOpCount(tn, ts) = o.arrOpCount(tn, ts) + 1;
                            o.arrPSCount(tn, ts) = o.arrPSCount(tn, ts) + o.cellQiteOps{tn, ts, mterm}(i).len;
                        end
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
                stuSummary.opCount(ts) = sum(o.arrOpCount(atrow, 1 : ts));
                stuSummary.PSCount(ts) = sum(o.arrPSCount(atrow, 1 : ts));
                % Accuracy
                stuSummary.accuracy(ts) = (stuSummary.energy(ts) - exactEigenvalue) / abs(exactEigenvalue);
                % Calculate endstate and fidelity
                endState = o.calculate_circuit(atrow, ts) * o.initialState;
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
            finalState = o.calculate_circuit(stuSummary.minPosition(atstep), atstep) * o.initialState;
            fprintf("QITE ground energy: %f \n", minE);
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
        function Circuit = calculate_circuit(o, timeDatapoint, trotterStep)
            % Initialize with identity matrix
            Circuit = eye(o.QK.PDH.twoN);
            for ts = 1 : trotterStep
                for mterm = 1 : o.QK.M
                    size_Ops = size(o.cellQiteOps{timeDatapoint, ts, mterm}, 2);
                    for i = 1 : size_Ops
                        % Calculate Pauli string in matrix form
                        opMatrix = o.cellQiteOps{timeDatapoint, ts, mterm}(i).matrixRep();
                        %o.cellQiteOps{timeDatapoint, ts, mterm}(i)
                        theta = o.cellThetas{timeDatapoint, ts, mterm}(i);
                    	Circuit = expm(- o.timeSeries(timeDatapoint) * theta * opMatrix) * Circuit;
                    end
                end
            end
        end % calculate_circuit
    end % methods
end % class definition