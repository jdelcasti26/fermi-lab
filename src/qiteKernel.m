%==========================================================================
% qiteKernel.m
%
%  Author: J. Del Castillo
%  Project: MT-QITE / fermi-lab
%  For theoretical details, see: arXiv:2512.10875
%
%==========================================================================
% Project   : QITE — Quantum Imaginary Time Evolution
% Class     : qiteKernel
% Author    : [Your Name]
% Reference : [Optional arXiv link]
%
% Description
% -----------
% Core kernel class implementing the measurement-based Quantum Imaginary
% Time Evolution (QITE) algorithm.
%
% This class encapsulates all low-level numerical and algebraic operations
% required by QITE, independently of higher-level control logic (time scans,
% Trotterization strategies, symmetry handling, etc.), which are implemented
% in wrapper classes such as qite, mtqite, etc.
%
% It is the common QITE kernel calculations for a number of QITE-based
% algorithms, like pure QITE, MT-QITE, ACQITE and sm-QITE
%
% The kernel operates on a reference quantum state |ψ⟩ and computes the
% measurement data required to build and solve the QITE linear system
%
%     SST · θ = b
%
% for a given Hamiltonian term, ansatz, and imaginary-time step Δt.
%
% Design Principles
% -----------------
% • Handle class inheriting from `hamiltonian`, ensuring shared hamiltonian partition
%  ansätze and efficient reuse across algorithmic layers.
% • Measurement-driven formulation: all quantities entering the linear
%   system are obtained from expectation values on the current stateVector.
% • Explicit separation between:
%     - measurement generation,
%     - linear-system construction,
%     - linear-system solution,
%     - state evolution.
%
% State Management
% ----------------
% The kernel maintains an internal reference `stateVector`, which:
% • is set or updated externally via `set_statevector()`,
% • serves as the state on which all QITE measurements are performed,
% • must remain unchanged if measurement reuse is intended.
%
% Measurement Handling and Reuse
% ------------------------------
% When `measure_state(mterm)` is invoked, a structure `measured(mterm)` is
% populated with all coefficients required to build the linear system:
%
%   measured(mterm).c1, c2   : normalization coefficients
%   measured(mterm).Mb1, Mb2 : independent vectors
%   measured(mterm).SST      : symmetric coefficient matrix
%
% These measurements are independent of Δt and may therefore be reused to
% solve multiple linear systems for different imaginary-time step sizes,
% provided the reference stateVector has not been updated.
%
% Measurement reuse is explicitly supported via the `reuse` structure and is
% exploited by higher-level algorithms (e.g. qite, mtqite) to reduce the
% total measurement cost.
%
% Linear System Solution
% ----------------------
% The linear system is solved by `solve_linear_system()` using either:
% • the Moore–Penrose pseudoinverse (pinv), or
% • a custom QR-based solver (`custom_pinv`) designed to minimize the number
%   of non-zero coefficients, thereby reducing circuit depth.
%
% A numerical tolerance EPSLS is used to identify and discard negligible
% coefficients.
%
% State Evolution
% ---------------
% The method `evolve_state()` constructs the QITE unitary
%
%     U(Δt) = ∏_i exp(-Δt · θ_i · P_i)
%
% from the non-zero solution coefficients and applies it to the current
% reference state, producing the evolved state |ψ'⟩.
%
% Supported QITE Orders
% ---------------------
% The kernel supports:
% • First-order QITE ("1st")
% • Second-order QITE ("2nd")
% • First-order with zero-order normalization ("1st0")
%
% Measurement operators for each order are precomputed at construction time.
%
% Intended Usage
% --------------
% This class is not meant to be used directly by the end user. It is the
% computational backbone of higher-level QITE workflows, which manage:
% • time-step scans,
% • Trotter steps,
% • Hamiltonian term ordering,
% • symmetry exploitation,
% • energy estimation and benchmarking.
%
%==========================================================================
classdef qiteKernel < hamiltonian
    properties
        % Inherited
        % EPSEIG, PDH, type, model, H(), M, HT, cellAZ{}
        EPSLS = 1.0e-7;     % Tolerance for low rank linear system inversion 
        options             % options struct defining:
                            % .order = ("1st", "2nd", "1st0")
                            % .Pinv = (true, false) if the Moore-Penrose
                            % pseudoinverse method is to be applied to invert linear system
        stateVector         % Reference stateVector on which qite measurements are performed 
        % Internal storage to be reused. 
        H2                  % Squared haimtonian. Array (1,...,M) of pauliStrings
        measured            % struct array with coefficients, vectors and matrix resulting from
                            % measurements given a stateVector (one struct per term 1, ..., M)
                            % .c1 .c2 Order 1 and 2 normalizing coefficients
                            % .Mb1 .Mb2 Order 1 and 2 measurements entering the linear system's independent coefficients
                            % .SST coefficient matrix of the linear system (symmetric)
        reuse               % Identical structure (see above) to store measurements that
                            % are to be reused, especially in qite mode,
                            % for the first Trotter step and the first hamiltonian term
        measureCount        % Array (1,...,M) with the number of linear system measurements in latest round
        % Cells with one entry per term in
        % the Hamiltonian model according to partition definition in constructor's varargin
        tItJ                % Array of opsCoefficients with 2D cell of pauliStrings with definition of operators -{t_I, t_J}
        tIH                 % Array of opsIndependents with 1D cell of pauliStrings with definition of operators [H, t_I]
        tIH2                % Array of opsIndependents with 1D cell of pauliStrings with definition of operators [H^2, t_I]
    end % properties
    methods
        % Class constructor
        % First argument is handle to pauliDictionary
        function o = qiteKernel(PDH, options, model, p, varargin)
            
            % Hamiltonian construction
            o@hamiltonian(PDH, model, p, varargin{:});
            
            % Order of QITE formulas can be 1, 2 or 1 with 0-order norm c
            if(options.order == "1st" || options.order == "2nd" ||  options.order == "1st0")
                o.options.order = options.order;
            else
                error("Invalide qite order. Must be '1st', '2nd' or '1st0'");
            end
            if(islogical(options.Pinv)) 
                o.options.Pinv = options.Pinv;
            else
                error("Invalid option.Pinv. It must be boolean");
            end       
            % Prepare measuring maps in pauliStrings and measurement operators
            o.H2 = pauliString.empty;           
            o.tIH = opsIndependent.empty;      
            o.tIH2 = opsIndependent.empty;         
            o.tItJ = opsCoefficient.empty; 
            % Measuring maps
            o.HT = o.HT.calcMap();
            for m = 1 : o.M
                % Measuring maps
                o.H(m) = o.H(m).calcMap();
                if(o.options.order == "2nd")
                    o.H2(m) = o.H(m) ^ 2;
                    o.H2(m) = o.H2(m).calcMap();
                end
                % Define measuremnt operator strings (maps defined in constructors)
                fprintf("Building strings tIH ...\n");
                o.tIH(m) = opsIndependent(o.PDH, o.cellAZ{m}, o.H(m));
                if(o.options.order == "2nd")
                    fprintf("Building strings tIH2 ...\n");
                    o.tIH2(m) = opsIndependent(o.PDH, o.cellAZ{m}, o.H2(m));
                end
                fprintf("Building strings tItJ ...\n");
                o.tItJ(m) = opsCoefficient(o.PDH, o.cellAZ{m});
            end % for m
            o.measured = struct.empty;
            o.measureCount = [];
        end % end constructor
        % Set statevector
        function o = set_statevector(o, newState)
            % Dimensonality check, stateVector mast be 2^no_qubits x 1
            size_newState = size(newState);
            if(size_newState(2) ~= 1 || size_newState(1) ~= o.PDH.twoN)
                error("Inconsistent dimensions for statevector");
            else
                % Renormalization for safety
                o.stateVector = newState / norm(newState);
            end            
        end
        % measureState : core method with QITE's measurements calculations
        % mterm is one of the 1, .... M terms in the partitioned hamiltonian
        function o = measure_state(o, mterm)
            %fprintf("Meas. state %d\n", mterm);
            % Measurement tracking
            start_count = o.PDH.measuresLS;
            % Allocate space for linear system's parameters (each hterm may have a different size)
            size_LS = size(o.cellAZ{mterm}, 2);
            %fprintf("Meas. state ansatz size %d\n", size_LS);
            o.measured(mterm).c1 = 0;
            o.measured(mterm).c2 = 0;
            o.measured(mterm).Mb1 = zeros(size_LS, 1);
            o.measured(mterm).Mb2 = zeros(size_LS, 1);
            o.measured(mterm).SST = zeros(size_LS, size_LS);

            if(o.options.order == "1st0")
                o.measured(mterm).c1 = 0;
                o.measured(mterm).c2 = 0;
            elseif(o.options.order == "1st") 
                o.measured(mterm).c1 = real(o.H(mterm).eValLS(o.stateVector));
                o.measured(mterm).c2 = 0;
            elseif(o.options.order == "2nd")
                o.measured(mterm).c1 = real(o.H(mterm).eValLS(o.stateVector));
                o.measured(mterm).c2 = real(o.H2(mterm).eValLS(o.stateVector)); 
            end
            %fprintf("Meas. state c1 %f\n", o.measured(mterm).c1);
    
            % Sweep rows
            for i = 1 : size_LS
                % Mb calculation
                % First order QITE formulas (also used in second order)
                o.measured(mterm).Mb1(i) = real(o.tIH(mterm).cellPS{i}.eValLS(o.stateVector));
                %%% DEBUG
                %if(i == 1 && mterm == 1)
                %    fprintf("Mb1(1)\n");
                %    o.tIH(mterm).cellPS{i}.disp();
                %end
                %%%
                % Second order QITE formulas
                if(o.options.order == "2nd") % Otherwise zero
                    o.measured(mterm).Mb2(i) = real(o.tIH2(mterm).cellPS{i}.eValLS(o.stateVector));
                end
                % Double loop to calculate matrix SST
                for j = 1 : size_LS
                    if(j == i) % Warning : in type = "spins" only i = 1, j = 1 needs to be calculated
                        o.measured(mterm).SST(i, j) = real(o.tItJ(mterm).cellPS{i, j}.eValLS(o.stateVector));
                    end
                    if(j > i)
                        o.measured(mterm).SST(i, j) = real(o.tItJ(mterm).cellPS{i, j}.eValLS(o.stateVector));
                    elseif(j < i)
                        o.measured(mterm).SST(i, j) = o.measured(mterm).SST(j, i);
                    end
                end
            end
            o.measureCount(mterm) = o.PDH.measuresLS - start_count;
            %fprintf("Count %d: %f\n", mterm, o.measureCount(mterm));
        end % end measureState
        % Calculate linear system independent term
        function b = independent_term(o, mterm, Dt)
            % Normalizing constant (dependent on time for both 1st and 2nd order, not for 1st0)
            c = 1 - 2 * Dt * o.measured(mterm).c1 + 2 * Dt ^ 2 *  o.measured(mterm).c2;
            % Independent term (dependet on time for second order only)   
            b = (c ^ (-0.5)) * ( o.measured(mterm).Mb1 - 0.5 * Dt *  o.measured(mterm).Mb2); % With 0.5 2nd-order correction
        end % end calcLSparams   
        % Solve linear system SST * theta = Mb
        % Return time Dt, non-zero thethas and associated operators from
        % the ansatz in qite_ops, so that the qite unitary is defined
        % applying a product series of exp(Dt * theta * qite_op) (see evolve_state())
        function [thetas, qiteOps] = solve_linear_system(o, mterm, Dt)
            % Calculate independent coefficients (with measured vectors for Dt)
            b = o.independent_term(mterm, Dt);
            % Solve equations with the Moore-Penrose pseudoinverse
            % algorithm for numerical stability or custom_pinv to minimize
            % the number of non-zero solutions (minimize circuit depth)
            if(o.options.Pinv == true)
                X = pinv(o.measured(mterm).SST) * b; % Tried and tested
            else
                X = o.custom_pinv(o.measured(mterm).SST, b, o.EPSLS);
            end
            id_nonzero = (abs(X) > o.EPSLS);
            %fprintf("LS %d t = %f: %s \n", mterm, Dt, strjoin(string(id_nonzero), ' '));
            thetas = X(id_nonzero)'; % Row vector
            qiteOps = o.cellAZ{mterm}(id_nonzero');  % Row array   
            %% DEBUG
            %if(mterm == 1)
                %SST = o.measured(mterm).SST;
                %Mb = b;
                %A = X;
                %[sza, szb] = size(SST);
                %for ia = 1 : sza
                %    for ib = 1 : szb
                %        fprintf("%f ", SST(ia, ib));
                %    end
                %    fprintf("Mb %f\n", Mb(ia));
                %end
                %for ia = 1 : sza
                %    fprintf("A %f\n", A(ia));
                %end
            %end
            %% DEBUG

        end
        % Function evolve_state : applies the obtained qite unitaries and
        % updates o.stateVector
        function evolvedState = evolve_state(o, Dt, thetas, qiteOps)
            no_operations = size(thetas, 2);
            newState = o.stateVector;
            for i = 1 : no_operations
                % Calculate Pauli string in matrix form
                opMatrix = qiteOps(i).matrixRep();
                UnitaryOp = expm(- Dt * thetas(i) * opMatrix);
                newState = UnitaryOp * newState; % 1.000 can be 2 for Zassenhaus formula
            end % for      
            % Upadate final statevector
            evolvedState = newState;
        end
        % Custom linear system solution for a non full rank matrix of coefficients A, independent
        % terms in column B, maximizing the number of zeroes in result X with tolerance tol
        function X = custom_pinv(o, A, B, tol)
            size_A = size(A);
            size_B = size(B);
            if(size_A(1) ~= size_A(2))
                error("Matrix A must be square");
            end
            if(size_B(1) ~= size_A(1))
                error("Independent parameters dimension mismatch");
            end
            [Q, R, P] = qr(A);

            % Solve system R * Xtil = inv(Q) * B
            % Find linearly dependent rows
            Qm1B = inv(Q) * B;
            Xtil = zeros(size_B(1), 1);
            cut = size_A(1);
            for i = size_A(1) : -1 : 1
                if(abs(Qm1B(i)) < tol)
                    cut = i;
                else
                    break;
                end
            end
            if(cut < size_A(1))
                cut = cut - 1;
            end
            Xtil(1 : cut) = inv(R(1 : cut, 1 : cut)) * Qm1B(1 : cut); % Rest below cut zeros
            X = P * Xtil;
        end % end custom_pinv       
    end % methods
end % classdef