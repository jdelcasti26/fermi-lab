% pauliDictionary.m
% PAULIDICTIONARY  Pauli-string maps, measurement cache and initial states.
%
%   PDH = pauliDictionary(nqubits) creates a Pauli dictionary object for
%   an nqubit system. The class:
%
%     • Stores one-qubit Pauli matrices I, X, Y, Z and some derived
%       one-qubit operators (Hadamard, creation, annihilation, number).
%     • Provides conversions between integer indices and Pauli strings,
%       e.g. 1 ↔ 'III', 2 ↔ 'IIX', ...
%     • Offers a Kronecker-product helper for building tensor-product
%       operators from one-qubit factors.
%     • Maintains a dictionary (DIC) of Pauli-string “maps” used to compute
%       expectation values efficiently, without constructing full 2^n×2^n
%       matrices.
%     • Maintains a dictionary (DICOM) that caches measurement results
%       (expectation values) of Pauli strings for a given state, together
%       with flags indicating whether a given measurement is “fresh”.
%     • Tracks separate measurement counters for different parts of the
%       algorithm (linear system, container generation, analysis).
%     • Provides convenience routines to build problem-specific initial
%       states for predefined model families (prepare_initial_state).
%
%   Internally, nqubit Pauli strings (e.g. 'XXYY') are represented as
%   compact labels, and their action on wavefunctions is encoded via maps
%   with entries in {+1, -1, +1i, -1i}. This avoids explicit construction
%   of dense Pauli matrices and enables fast expectation-value evaluation.
%
%   Example:
%
%     PDH = pauliDictionary(4);
%     s   = PDH.numToString(2);      % returns 'IIX' for 4 qubits
%     n   = PDH.stringToNum('IIX');  % returns 2
%
%   See also: qiteKernel, pauliString, hamiltonian, qite
classdef pauliDictionary < handle
    properties
        nqbs        % Number of qubits
        twoN        % Operator size : 2 ^ nqbs
        fourN       % Linear system maximum size / Pauli basis size : 4 ^ nqbs
        I           % Pauli matrix I
        X           % Pauli matrix X
        Y           % Pauli matrix Y
        Z           % Pauli matrix Z
        % Additional 1-qubit operators
        A           % Hadamard gate A
        C           % Creation operator (a')
        D           % Anihilation operator (a)
        N           % Number operator (a'a)
        % Dictionary for Pauli string operators in index + coefficient form
        DIC
        % Dictionary for measurements
        DICOM
        % Class internal properties
        PAULI       % Full pauli basis cell: {I, X, Y, Z} in matrix form
        % Measurement count: 
        % - involved in linear system generation (measuresLS)
        % - container generation (measuresCO)
        % - Multiple-time sweeps in analyzeh (measuresAZ)
        measuresLS
        measuresCO
        measuresAZ
    end % end properties
    methods
        function o = pauliDictionary(nqubits)
            o.nqbs = nqubits;
            o.twoN = 2 ^ o.nqbs; % Operator size
            o.fourN = 4 ^ o.nqbs; % Linear system / Pauli basis size
            % Hardcoded Pauli matrices
            o.I = eye(2);
            o.X = [0 1; 1 0];
            o.Y = [0 -1j; 1j 0];
            o.Z = [1 0; 0 -1];
            % Internal structures
            o.PAULI = cell(1,4);
            o.PAULI{1} = o.I;
            o.PAULI{2} = o.X;
            o.PAULI{3} = o.Y;
            o.PAULI{4} = o.Z;
            % Hardcoded 1 qubit operators. Hadamard, creation, anihilation and number
            o.A = (1 / sqrt(2)) * [1 1; 1 -1];
            o.C = [0 0;1 0]; % (1/2)(X - iY)
            o.D = [0 1;0 0]; % (1/2)(X + iY)
            o.N = [0 0; 0 1]; % (1/2)(I - Z)
            % Pauli dictionaries
            % Dictionary for maps
            o.DIC = dictionary(string.empty, struct.empty);
            % Dictionary for measurements and control flags
            o.DICOM = dictionary(string.empty, struct.empty);
            % Measurement count
            o.measuresLS = 0;
            o.measuresCO = 0;
            o.measuresAZ = 0;
        end
        % NUMTOSTRING Convert integer index to Pauli string.
        %
        %   s = numToString(o, i) returns the Pauli string corresponding to
        %   the index i, using a base-4 encoding over {I,X,Y,Z}. For
        %   example, for nqbs = 3:
        %     i = 1  → 'III'
        %     i = 2  → 'IIX'
        %     i = 3  → 'IIY'
        %   and so on.
        function s = numToString(o, i)
            % Calculate sequence number in base 4 (string format, st)
            st = dec2base(i - 1, 4);
            % Pad zeroes on the left up to no. of qubits
            lst = length(st);
            for c = (lst + 1) : o.nqbs
                st = strcat('0', st); 
            end    
            % Initialize resulting string
            s = '';
            for q = 1 : o.nqbs
                PauliIndex = str2num(st(q)) + 1;
                if(PauliIndex == 1)
                    s = strcat(s, 'I');
                elseif(PauliIndex == 2)
                    s = strcat(s, 'X');
                elseif(PauliIndex == 3)
                    s = strcat(s, 'Y');    
                elseif(PauliIndex == 4)
                    s = strcat(s, 'Z');
                else
                    error("Internal error");
                end
            end
        end % end numToString
        % STRINGTONUM Convert Pauli string to integer index.
        %
        %   n = stringToNum(o, s) returns the integer index corresponding
        %   to the Pauli string s, using the same encoding as NUMTOSTRING.
        %   Indexing is 1-based to match MATLAB conventions.        
        function n = stringToNum(o, s)
            % Return value allocation
            n = 0;
            for c = 1 : o.nqbs
                % I = 1/0, X = 2/1, Y = 3/2, Z = 4/3 
                if (s(c) == 'X')
                    n = n + 1 * 4^(o.nqbs - c);
                elseif(s(c) == 'Y')
                    n = n + 2 * 4^(o.nqbs - c);
                elseif(s(c) == 'Z')
                    n = n + 3 * 4^(o.nqbs - c);  
                end
            end 
            % Array indexes in Matlab begin with 1!
            n = n + 1;  
        end
        % KP  Kronecker product of a row cell array of matrices.
        %
        %   K = kp(o, MultsCell) returns the tensor product
        %       MultsCell{1} ⊗ MultsCell{2} ⊗ ... ⊗ MultsCell{end}.
        %
        %   MultsCell must be a 1×N cell array of square matrices.
        function kmatrix = kp(o, MultsCell)
            szMc = size(MultsCell);
            if(szMc(1) > 1)
                error("Multiplicands array must be a row cell");
            else
                % number of multiplicands
                Nm = szMc(2);
            end
            % Initialize result with last multiplicans
            kmatrix = MultsCell{Nm};
            for i = Nm - 1 : -1 : 1
                kmatrix = kron(MultsCell{i}, kmatrix);
            end
        end % end kp   
        % RESETMEASURESLS Reset linear-system measurement counter to zero.
        % RESETMEASURESCO Reset enrgy measurement counter.
        % RESETMEASURESAZ Reset variational measurement counter.
        %
        % INCMEASURESLS / CO / AZ increment the respective counters.
        %
        % FLUSHMEASUREMENTS marks all cached measurements as "not fresh"
        % in the DICOM dictionary, so they will be recomputed for the next
        % updated state vector.        
        function o = resetMeasuresLS(o)
            o.measuresLS = 0;
        end
        function o = resetMeasuresCO(o)
            o.measuresCO = 0;
        end
        function o = resetMeasuresAZ(o)
            o.measuresAZ = 0;
        end        
        function o = incMeasuresLS(o)
            o.measuresLS = o.measuresLS + 1;
        end  
        function o = incMeasuresCO(o)
            o.measuresCO = o.measuresCO + 1;
        end           
        function o = incMeasuresAZ(o)
            o.measuresAZ = o.measuresAZ + 1;
        end     
        % Function to erase all measurements after the state vector has
        % been updated
        function o = flushMeasurements(o)
            K = keys(o.DICOM);
            sz = size(K);
            stuini.fresh = false;
            stuini.value = 0;
            for i = 1 : sz(1)
                o.DICOM(K(i)) = stuini;
            end
        end
        % PREPARE_INITIAL_STATE Build a model-specific initial state vector.
        %
        %   state = prepare_initial_state(o, modelFam) returns a normalized
        %   state vector of length 2^nqbs corresponding to a predefined
        %   initial state for the selected model family, e.g.:
        %
        %     • "H2"       with nqbs = 4
        %     • "ISIOB04"  with nqbs = 4
        %     • "ISIOB06", "ISYOB06", "HEYOB06" with nqbs = 6
        %     • "HEYOB08", "H4CHAIN..." with nqbs = 8
        %
        %   If the requested (modelFam, nqbs) pair is unknown, an error is
        %   raised. The returned state is always normalized.
        function state = prepare_initial_state(o, modelFam)
            state = zeros(o.twoN, 1);
            % Four-qubit models
            if(modelFam == "H2" && o.nqbs == 4)
                state(4, 1) = 1;
            elseif(modelFam == "ISIOB04" && o.nqbs == 4)
                state = ones(o.twoN, 1);
            % Six-qubit models
            elseif((modelFam == "ISIOB06" || modelFam == "ISYOB06") && o.nqbs == 6)
                state = ones(o.twoN, 1);
            elseif(modelFam == "HEYOB06" && o.nqbs == 6)
                state = [0; -176.7767; 0; 0; 51.7767; 0; 0; 88.3883; 114.2767; 0; 0; 132.5825; 0; 145.5267; 154.6796; 0; 161.1517; 0; 0; 165.7282; 0; 168.9642; 171.2524; 0; 0; 172.8704; 174.0146; 0; 174.8236; 0; 0; 424.8236; -424.8236; 0; 0; -174.8236; 0; -174.0146; -172.8704; 0; 0; -171.2524; -168.9642; 0; -165.7282; 0; 0; -161.1517; 0; -154.6796; -145.5267; 0; -132.5825; 0; 0; -114.2767; -88.3883; 0; 0; -51.7767; 0; 0; 176.7767; 0];
            elseif(modelFam == "HUBOB06" && o.nqbs == 6)
                state(4, 1) = 1;
                state(7, 1) = 1;
                state(10, 1) = 1;
                state(13, 1) = 1;
                state(19, 1) = 1;
                state(25, 1) = 1;
                state(34, 1) = 1;
                state(37, 1) = 1;
                state(49, 1) = 1;
            % Eight-qubit models
            elseif(modelFam == "HEYOB08" && o.nqbs == 8)
                state = [0.000000;0.000000;0.000000;0.125000;0.000000;0.000000;-0.125000;0.000000;0.000000;0.000000;0.125000;0.000000;0.000000;0.000000;0.000000;-0.125000;0.000000;0.000000;-0.125000;0.000000;0.000000;0.000000;0.000000;0.125000;0.000000;0.000000;0.000000;-0.125000;0.000000;0.000000;0.125000;0.000000;0.000000;0.000000;0.125000;0.000000;0.000000;0.000000;0.000000;-0.125000;0.000000;0.000000;0.000000;0.125000;0.000000;0.000000;-0.125000;0.000000;0.000000;0.000000;0.000000;-0.125000;0.000000;0.000000;0.125000;0.000000;0.000000;0.000000;-0.125000;0.000000;0.000000;0.000000;0.000000;0.125000;0.000000;0.000000;0.125000;0.000000;0.000000;0.000000;0.000000;-0.125000;0.000000;0.000000;0.000000;0.125000;0.000000;0.000000;-0.125000;0.000000;0.000000;0.000000;0.000000;-0.125000;0.000000;0.000000;0.125000;0.000000;0.000000;0.000000;-0.125000;0.000000;0.000000;0.000000;0.000000;0.125000;0.000000;0.000000;0.000000;0.125000;0.000000;0.000000;-0.125000;0.000000;0.000000;0.000000;0.125000;0.000000;0.000000;0.000000;0.000000;-0.125000;0.000000;0.000000;-0.125000;0.000000;0.000000;0.000000;0.000000;0.125000;0.000000;0.000000;0.000000;-0.125000;0.000000;0.000000;0.125000;0.000000;0.000000;0.125000;0.000000;0.000000;-0.125000;0.000000;0.000000;0.000000;0.125000;0.000000;0.000000;0.000000;0.000000;-0.125000;0.000000;0.000000;-0.125000;0.000000;0.000000;0.000000;0.000000;0.125000;0.000000;0.000000;0.000000;-0.125000;0.000000;0.000000;0.125000;0.000000;0.000000;0.000000;0.125000;0.000000;0.000000;0.000000;0.000000;-0.125000;0.000000;0.000000;0.000000;0.125000;0.000000;0.000000;-0.125000;0.000000;0.000000;0.000000;0.000000;-0.125000;0.000000;0.000000;0.125000;0.000000;0.000000;0.000000;-0.125000;0.000000;0.000000;0.000000;0.000000;0.125000;0.000000;0.000000;0.125000;0.000000;0.000000;0.000000;0.000000;-0.125000;0.000000;0.000000;0.000000;0.125000;0.000000;0.000000;-0.125000;0.000000;0.000000;0.000000;0.000000;-0.125000;0.000000;0.000000;0.125000;0.000000;0.000000;0.000000;-0.125000;0.000000;0.000000;0.000000;0.000000;0.125000;0.000000;0.000000;0.000000;0.125000;0.000000;0.000000;-0.125000;0.000000;0.000000;0.000000;0.125000;0.000000;0.000000;0.000000;0.000000;-0.125000;0.000000;0.000000;-0.125000;0.000000;0.000000;0.000000;0.000000;0.125000;0.000000;0.000000;0.000000;-0.125000;0.000000;0.000000;0.125000;0.000000;0.000000;0.000000];
            elseif(startsWith(modelFam, "H4CHAIN") && o.nqbs == 8)
                state(16, 1) = 1; % Hartree-Fock, four electrons 
            else
                error("Requested initial state unknown");
            end
            % Normalize
            state = state / norm(state);
        end
    end % end methods
end