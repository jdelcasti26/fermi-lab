% =========================================================================
% opsIndependent.m
%
%  Author: J. Del Castillo
%  Project: MT-QITE / fermi-lab
%  For theoretical details, see: arXiv:2512.10875
%
% Description:
%   This class constructs and stores the Pauli-string objects resulting from
%   the product between each ansatz operator and a single fixed operator.
%   It is designed for terms that are independent across the ansatz index,
%   such as:
%
%       • ansatz_i × H
%       • ansatz_i × H^2
%       • ansatz_i × O
%
%   for a given Pauli operator O.
%
%   In the context of QITE, this structure is used to assemble the
%   right-hand-side vectors of the linear system arising from imaginary-time
%   evolution, where each ansatz operator is paired with Hamiltonian-derived
%   operators.
%
%   Only the real (Hermitian) part of the commutator is retained, since
%   imaginary components do not contribute to measurable expectation values
%   and therefore do not require quantum measurements.
%
% Key design choices:
%   • One-to-one mapping between ansatz operators and resulting Pauli strings
%   • Storage in a 1D cell array (cellPS{i})
%   • Pauli-string maps are precomputed for efficient reuse
%   • Measurement-free imaginary terms are explicitly discarded
%
% Main properties:
%   PDH         : Handle to pauliDictionary (shared Pauli algebra utilities)
%   size_ansatz : Number of operators in the ansatz
%   cellPS      : Cell array of resulting Pauli strings (length = size_ansatz)
%
% Typical usage:
%   ops = opsIndependent(PDH, ansatz, op);
%
%   where:
%     ansatz : row array of pauliString objects
%     op     : single pauliString (e.g. H or H^2 term)
%
% Notes:
%   This class complements opsCoefficient, which handles pairwise ansatz–
%   ansatz products. Together, they provide the algebraic backbone for
%   constructing QITE linear systems.
%
% =========================================================================

classdef opsIndependent
    properties
        % Handle to Pauli dictionary
        PDH    
        % General properties
        size_ansatz       % Series index length (new size linear system)
        cellPS            % Cell of resulting Pauli strings
    end
    methods
        % Constructor: ansatz is an array of pauliStrings, op is single pauliString
        % First argument is handle to Pauli dictionary
        function o = opsIndependent(PDH, ansatz, op)

            [m, o.size_ansatz] = size(ansatz);
            if(m ~= 1)
                error("Expected row array for relevant pauliStrings (operators)");
            end

            o.PDH = PDH;

            % Allocate array of pauliStrings
            o.cellPS = pauliString.empty;
            for i = 1 : o.size_ansatz
                % Pauli operator i
                tI = ansatz(i);
                % Product ansatz Operator(i) * given Operator
                o.cellPS{i} = realPS(commutator(op, tI)); % Ignore imaginary parts (no measurement needed)
                %fprintf("\nResulting PS (sig op): %s\n", o.cellPS{i}.string());
                % Calculate associated map
                if(o.cellPS{i}.len > 0)
                    o.cellPS{i} = o.cellPS{i}.calcMap();
                end
            end
        end
    end
end