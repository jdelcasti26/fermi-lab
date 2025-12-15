%==========================================================================
% opsCoefficient.m
%
%  Author: J. Del Castillo
%  Project: MT-QITE / fermi-lab
%  For theoretical details, see: arXiv:2512.10875
%
%
% Purpose
% -------
% This class constructs and stores the operators entering the
% QITE linear system, obtained from pairwise combinations of ansatz operators.
% For a given ansatz term {t_I}, it evaluates the real Pauli-string contributions
% resulting from the (anti)commutator structure required by imaginary-time
% evolution:
%
%     C_{IJ} ~ Re{ - {t_I, t_J} }
%
% Only the upper-triangular part (J >= I) is computed and stored, exploiting
% the symmetry of the resulting coefficient matrix.
%
% Each matrix element corresponds to a Pauli string array (or empty string) that
% contributes to the linear system defining the QITE update directions.
%
% Design notes
% ------------
% • The resulting Pauli strings are stored in a 2D cell array `cellPS`,
%   indexed as (i, j) with j >= i.
% • Complex (purely imaginary) contributions are discarded by construction,
%   as only real-valued operators enter the QITE linear system.
% • Each resulting Pauli string is immediately mapped to its internal
%   dictionary representation via `calcMap()` for efficient reuse.
% • This class is typically instantiated once per ansatz update and reused
%   during linear-system assembly inside qite_kernel.
%
% Intended usage
% --------------
% This class is a low-level building block used internally by qite_kernel
% to assemble the operator overlap matrix appearing in the QITE linear system.
% It is not intended to be used directly by high-level drivers.
%
%==========================================================================
classdef opsCoefficient
    properties
        PDH
        size_ansatz        % Ansatz length (new size linear system)
        cellPS             % 2D cell of resulting Pauli strings
        %pStr              % 2D array of resulting Pauli strings
        % Dictionary for resulting sigma * sigma string, and associated
        % array ARRSS (size_ansatz * size_ansatz) to store resulting pauliString in
        % signed index
        %DICSS       % Key is an integer with the Pauli string index sequence times phase (only real +1 or -1)
        %ARRSS
        %PSEUR      % WIP Pauli string
    end
    methods
        % Constructor: 
        % First argument is handle to pauliDictionary
        function o = opsCoefficient(PDH, ansatz)
            [m, o.size_ansatz] = size(ansatz);
            if(m ~= 1)
                error("Expected row vector for relevant ansatz operators");
            end

            o.PDH = PDH;
            % WARNING int storage for 15 qubits maximum
            %o.DICSS = dictionary(int32.empty, pauliString.empty); 
            %o.ARRSS = zeros(o.size_ansatz, o.size_ansatz, 'int32');

            % Allocate array of pauliStrings
            %o.pStr = pauliString.empty;
            o.cellPS = pauliString.empty;
            for i = 1 : o.size_ansatz
                % Pauli operator i
                tI = ansatz(i);                
                for j = i : o.size_ansatz
                    % Pauli operator j
                    tJ = ansatz(j);
                    %o.pStr{i, j} = realPS(SIG_i * SIG_j); % If resulting product has an i phase ignore
                    o.cellPS{i, j} = realPS(-1 * anticommutator(tI, tJ)); % If resulting product has an i phase ignore
                    o.cellPS{i, j} = o.cellPS{i, j}.calcMap();
                    % Calculate resulting index times phase (resulting
                    % Pauli string only has one term)
                    % Empty check
                    %if(PSR.len == 0)
                        % If empty pauliString assign 'fake' index = 0
                    %    indPSR = int32(0);
                    %else
                        % Commented behavior for single Pauli string ansatz
                        %%indPSR = int32(PSR.coef(1) * o.PH.stringToNum(PSR.pStr{1}));
                        % End of commented behaviour
                        % NEW Behaviour for multiple Pauli string ansatz
                        %indPSR = int32(1000 * i + j);
                    %end
                    % Create entry in ARRSS
                    %o.ARRSS(i, j) = indPSR;
                    %if(~isKey(o.DICSS, indPSR))
                    %    % Calculate associated map and add to dictionary
                    %    if(PSR.len > 0)
                    %        PSR = PSR.calcMap();
                    %    end
                    %    o.DICSS(indPSR) = PSR;
                    %end
                end
            end
            %timesigma = toc;
            %fprintf("Elapsed %d seconds\n", timesigma);
        end % end constructor
    end
end