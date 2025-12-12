classdef acqite
    properties
        % Cosntants
        EPSTIME = 1.0E-12; % Time steps accuracy
        % Runtime configuration
        QK                 % Handle to qite_kernel 
        initialState       % Initial state for algorithm
        timeSeries         % Vector containing (positive) time step sizes to be scanned
        noTIMES            % Number of different time step sizes
        noTS               % Number of Trotter steps
        % Simulation results
        cellQiteOps        % Circuit definition for 1 time, along all Trotter steps, all terms in partition
        cellThetas         % Circuit parameters for the above
        cellAccQiteOps     % Accumulated Qite ops. per Trotter step, for optimum time interval
        cellAccThetas
        measureCountLS     % Array with performed measurements to obtain linear system per time, per Trotter step per term
        measureCountEE     % Array with performed measurements to obtain estimated energy per time, per Trotter step
        arrayE             % Energy array (per time datapoint, per Trotter step)
        arrOpCount         % Operator (in the ansatz) count, per time datapoint, per Trotter step 
        arrPSCount         % Pauli string count, per time datapoint, per Trotter step 
    end % properties
    methods
        function o = acqite(QK, initialState, timeSeries, noTS)
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
            % Check that timeSeries starts with a finite value
            if(timeSeries(1) == 0)
                error("Time series must start with a non-zero value");
            end
            % Check that timeSeries is uniformly distributed
            if (~all(abs(diff(timeSeries) - median(diff(timeSeries))) <= o.EPSTIME))
                error("Expected a uniformly distributed time series");
            end
            % Sort time step sizes for safety
            o.timeSeries = sort(unique(timeSeries));
            o.noTIMES = size_timeSeries(2);
            o.noTS = noTS;
            % Allocate container structures. Only first time datapoint is measured and stored
            o.cellQiteOps = cell(1, o.noTS, o.QK.M);
            o.cellThetas = cell(1, o.noTS, o.QK.M);
            % Accumulated cells detailed by Trotter step only
            o.cellAccQiteOps = cell(o.noTS);
            o.cellAccThetas = cell(o.noTS);            
            o.measureCountLS = zeros(o.noTIMES, o.noTS, o.QK.M);
            o.measureCountEE = zeros(o.noTIMES, o.noTS);
        end
        function o = scan_step(o, stateVector, trotterStep)
            o.QK = o.QK.set_statevector(stateVector);
            o.QK.PDH = o.QK.PDH.flushMeasurements();
            tn = 1; % Only the first datapoint is measured, all the rest are proportional
            for mterm = 1 : o.QK.M % For each term in the partition
                % Measure qite's linear system
                o.QK = o.QK.measure_state(mterm);
                % Inform performed number of distinct measurements
                o.measureCountLS(tn, trotterStep, mterm) = o.QK.measureCount(mterm);
                % Obtain RTE configuration
                [thetas, qiteOps] = o.QK.solve_linear_system(mterm, o.timeSeries(tn));
                % Inform circuit parameters in container
                o.cellQiteOps{tn, trotterStep, mterm} = qiteOps;
                o.cellThetas{tn, trotterStep, mterm} = thetas;
                % Evolve state
                evolvedState = o.QK.evolve_state(o.timeSeries(tn), thetas, qiteOps);
                % Updatet state implies stored measurement values are no longer valid
                o.QK = o.QK.set_statevector(evolvedState);
                o.QK.PDH = o.QK.PDH.flushMeasurements();
            end % for mterm
        end % function scan_step
        function [o, stuSummary, endState] = estimate_energy(o)
            o.arrayE = zeros(o.noTIMES, o.noTS);
            o.arrOpCount = zeros(o.noTS);
            o.arrPSCount = zeros(o.noTS);

            stateVector = o.initialState;
            for ts = 1 : o.noTS
                o = o.scan_step(stateVector, ts);
                [o, arrStuSummary(ts), stateVector] = o.estimate_energy_step(ts, stateVector);
            end
            stuSummary = table2struct(struct2table(arrStuSummary), 'ToScalar', true);
            % Cumulate depths and costs
            stuSummary.opCount = cumsum(stuSummary.opCount);
            stuSummary.PSCount = cumsum(stuSummary.PSCount);
            stuSummary.costLS = cumsum(stuSummary.costLS);
            stuSummary.costEE = cumsum(stuSummary.costEE);            
            endState = stateVector;
        end 
        function [o, stuSummary, stateVector] = estimate_energy_step(o, trotterStep, stateVector)
            % Accumulate by mterm available QiteOps
            accQiteOps = o.cellQiteOps{1, trotterStep, 1}; % First time datapoint, first mterm
            accThetas = o.cellThetas{1, trotterStep, 1};
            for mterm = 2 : o.QK.M
                [accQiteOps, accThetas] = o.accumulate_qite_ops(accQiteOps, accThetas, o.cellQiteOps{1, trotterStep, mterm}, o.cellThetas{1, trotterStep, mterm});
            end
            for tn = 1 : o.noTIMES
                evolvedState = o.calculate_circuit(accQiteOps, accThetas) * stateVector;
                o.QK.PDH = o.QK.PDH.flushMeasurements();
                measuresRef = o.QK.PDH.measuresCO;
                o.arrayE(tn, trotterStep) = o.QK.HT.eValCO(evolvedState);
                o.measureCountEE(tn, trotterStep) = o.QK.PDH.measuresCO - measuresRef;
                if(tn > 1) % Multiply the unitaries
                    if((o.arrayE(tn, trotterStep) > o.arrayE(tn - 1, trotterStep)) || tn == o.noTIMES)
                        o.cellAccQiteOps{trotterStep} = accQiteOps;
                        o.cellAccThetas{trotterStep} = accThetas;
                        stuSummary.energy = min(o.arrayE(tn, trotterStep), o.arrayE(tn - 1, trotterStep));
                        if(o.arrayE(tn, trotterStep) > o.arrayE(tn - 1, trotterStep))
                            stuSummary.minPosition = tn - 1;
                        else
                            stuSummary.minPosition = tn;
                        end
                        break;
                    end
                    [accQiteOps, accThetas] = o.accumulate_qite_ops(accQiteOps, accThetas, accQiteOps, accThetas);
                end
            end
            % Count accumulated Qite operations and store them in cell
            size_ops = size(accQiteOps, 2);
            o.cellAccQiteOps{trotterStep} = accQiteOps;
            o.cellAccThetas{trotterStep} = accThetas;
            stuSummary.opCount = size_ops;
            stuSummary.PSCount = 0;
            for op = 1 : size_ops
                stuSummary.PSCount = stuSummary.PSCount + accQiteOps(op).len;
            end
            % Accuracy and fidelity
            % Reference values
            exactEigval = o.QK.exact_eigenvalue();
            exactEigvecs = o.QK.exact_eigenvectors();
            stuSummary.accuracy = (stuSummary.energy - exactEigval) / abs(exactEigval);
            % end state
            stateVector = o.calculate_circuit(accQiteOps, accThetas) * stateVector;
            stuSummary.fidelity = (stateVector' * exactEigvecs) * (stateVector' * exactEigvecs)';
            % stuSummary.minPosition already informed
            % Inform cost
            stuSummary.costLS = sum(sum(o.measureCountLS(:, trotterStep, :)));
            stuSummary.costEE = sum(o.measureCountEE(:, trotterStep));

            fprintf("AC-QITE step %d, energy: %f \n", trotterStep, stuSummary.energy);
            fprintf("Accuracy: %14.5e \n", stuSummary.accuracy);
            fprintf("Depth in implemented operators %d and strings: %d \n", stuSummary.opCount, stuSummary.PSCount);
            fprintf("Fidelity: %14.5e \n", stuSummary.fidelity);
            fprintf("Linear system %d, energy %d and total %d mesurements.\n", stuSummary.costLS , stuSummary.costEE, stuSummary.costLS + stuSummary.costEE);
            fprintf("Minimum found at datapoint(s) no.: %s \n", strjoin(string(stuSummary.minPosition), ' '));
        end
        function Circuit = calculate_circuit(o, qiteOps, thetas)
            % Initialize with identity matrix
            Circuit = eye(o.QK.PDH.twoN);
            size_Ops = size(qiteOps, 2);
            for i = 1 : size_Ops
                % Calculate Pauli string in matrix form
                opMatrix = qiteOps(i).matrixRep();
                theta = thetas(i);
                Circuit = expm(- o.timeSeries(1) * theta * opMatrix) * Circuit;
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