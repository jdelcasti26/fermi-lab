% Class to define the Pauli definitions: Pauli basis, dictionary and their associated methods:
% - tensor product multiplication
% - string to index translation ???
% - index to string translation
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
        % Convert series index into string (e.g. 1 = 'III', 2 = 'IIX', 3 = 'IIY', ..., 5 = 'IXI' ...
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
        % function [kmatrix] = kp(MultiplicandsCell) returns the tensor
        % product of matrices arg1 x arg2 x arg3 x arg4 ... defined
        % in MultsCell
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
        % Convenience function with some sample initial states for a
        % selected set of model families
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