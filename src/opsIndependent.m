% Custom class to store Pauli strings resulting fromHM ansatz(ind) * Operator
% where ind is the Pauli string in the ansatz arrsizeay
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