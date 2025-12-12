% pauliString : class defines two parallel arrays / cells to represent a Pauli string 
% such as 0.5 * XXYYZZ - 1 * XIIIIX + (1/3) * ZZZZZZ
% The first array stores the complex coefficients, the second (cell) the Pauli operators in string form.
% There are two optional arrays to represent the same operators in numerical form (1 = I, 2 = X, 3 = Y, 4 = Z)
% and in matrix form. An additional optional matrix stores the total sum(coef * op)
% This class overrides methods '+', '*' and '^2'
classdef pauliString
    properties
        % Handle to Pauli dictionary and status control
        PH             
        updateNeeded
        % Mandatory properties
        len            % Number of terms in Pauli string
        coef           % Complex coefficients (complex array)
        pStr           % Pauli operator in string form (cell)
    end % end properties
    methods
        % Class constructor. First argument is dictionary handle
        % In order to build a null object type pauliString([0],'IIII');
        function o = pauliString(h, coeffs, strings)
            szCoef = size(coeffs);
            szPStr = size(strings);
            if(szCoef(2) ~= 1 || szPStr(2) ~= 1 || szCoef(1) ~= szPStr(1))
                error("Inconsistent coefficient or string sizes to define pauliString");
            end
            o.PH = h;
            o.updateNeeded = true;
            o.len = szCoef(1);
            o.coef = coeffs;
            o.pStr = strings;
            % Filter string length to match number of qubits
            for i = 1 : o.len
                o.pStr{i, 1} = o.pStr{i, 1}(1 : o.PH.nqbs);
                % Check string consistency
                for j = 1 : o.PH.nqbs
                    if(o.pStr{i}(j) ~= 'I' && o.pStr{i}(j) ~= 'X' && o.pStr{i}(j) ~= 'Y' && o.pStr{i}(j) ~= 'Z')
                        error("Invalid character in Pauli string definition");
                    end
                end
            end
            % Filter duplicates
            o = o.intnSort();
        end % end class constructor
        % Internal sorting of the Pauli string
        % Function sorts charcter strings and sums coefficients if duplicates
        function o = intnSort(o)
            [tmpPStr, idx] = sort(o.pStr);
            tmpCoef = o.coef(idx);
            % Check if there are duplicates and sum their coefficients
            [newPStr, idx2] = unique(tmpPStr);

            % Sweep tmpCoef bottom-up and push coefficient's value one step
            % up if row is a duplicate (inexistent in idx2)
            for i = o.len : -1 : 2
                if(~ismember(i, idx2))
                    tmpCoef(i - 1) = tmpCoef(i - 1) + tmpCoef(i);
                end
            end
            newCoef = tmpCoef(idx2);
            % Filter entries with 0 coefficients
            newPStr = newPStr(newCoef ~= 0);
            newCoef = newCoef(newCoef ~= 0);
            [o.len, ~] = size(newCoef);
            o.pStr = newPStr;
            o.coef = newCoef;
        end % end intnSort
        % Internal multiplication of two character strings (with phases)
        function [ph, chs] = intnMult(o, chs1, chs2)
            % Initialize resulting string and phase
            chs = '';
            ph = 1;
            for qb = 1 : o.PH.nqbs
                if(chs1(qb) == chs2(qb))
                    chs = strcat(chs, 'I');
                    % ph unchanged
                elseif(chs1(qb) == 'I')
                    chs = strcat(chs, chs2(qb));
                    % ph unchanged
                elseif(chs2(qb) == 'I')
                    chs = strcat(chs, chs1(qb));  
                    % ph unchanged
                elseif(chs1(qb) == 'X' && chs2(qb) == 'Y')
                    chs = strcat(chs, 'Z');
                    ph = ph * 1j;
                elseif(chs1(qb) == 'Y' && chs2(qb) == 'Z')
                    chs = strcat(chs, 'X'); 
                    ph = ph * 1j;
                elseif(chs1(qb) == 'Z' && chs2(qb) == 'X')
                    chs = strcat(chs, 'Y');
                    ph = ph * 1j;
                elseif(chs1(qb) == 'Y' && chs2(qb) == 'X')
                    chs = strcat(chs, 'Z');
                    ph = - ph * 1j;
                elseif(chs1(qb) == 'Z' && chs2(qb) == 'Y')
                    chs = strcat(chs, 'X');    
                    ph = - ph * 1j;    
                elseif(chs1(qb) == 'X' && chs2(qb) == 'Z')
                    chs = strcat(chs, 'Y');
                    ph = - ph * 1j; 
                end
            end
        end % end intMult
        % Internal function to check if two character (Pauli) strings commute
        % They commute if non-coinciding characters different from I are even
        function iFlagComm = intnComm(o, chs1, chs2)
            counter = 0;
            for qb = 1 : o.PH.nqbs
                if(~(chs1(qb) == 'I' || chs2(qb) == 'I' || (chs1(qb) == chs2(qb))))
                    counter = counter + 1;
                end
            end
            if(mod(counter, 2) == 0)
                iFlagComm = 1;
            else
                iFlagComm = 0;
            end
        end % end intComm
        % function commStructure extracts commuting matrix of the component
        % strings of a Pauli string
        function [S, SMat] = commStructure(o)
            S = cell(0, 3);
            row = 0;
            SMat = eye(o.len);
            for i = 1 : o.len
                for j = i + 1 : o.len
                    row = row + 1;
                    S{row, 1} = strcat(sprintf('%d', i), " ", o.pStr{i});
                    S{row, 2} = strcat(sprintf('%d', j), " ", o.pStr{j});
                    S{row, 3} = o.intnComm(o.pStr{i}, o.pStr{j});
                    SMat(i, j) = S{row, 3};
                    SMat(j, i) = S{row, 3};
                end
            end
        end
        % Overriding operators
        % Sum of two Pauli strings
        function s = plus(obj1, obj2)
            totalPStr = [obj1.pStr; obj2.pStr];
            totalCoef = [obj1.coef; obj2.coef];
            if(~isempty(obj1))
                s = pauliString(obj1.PH, totalCoef, totalPStr);
            elseif(~isempty(obj2))
                s = pauliString(obj2.PH, totalCoef, totalPStr);
            else
                error("Unable to sum two empty pauliStrings");
            end
        end % end plus
        % Multiplication of two Pauli strings
        function m = mtimes(obj1, obj2)
            % Allow numeric * object and object * numeric (scalar scaling).
            if(isa(obj1,'pauliString') && isnumeric(obj2))
                if(obj1.len == 0)
                    m = obj1;
                else
                    m = pauliString(obj1.PH, obj1.coef * obj2, obj1.pStr);
                end
            elseif(isnumeric(obj1) && isa(obj2,'pauliString'))
                if(obj2.len == 0)
                    m = obj2;
                else
                    m = pauliString(obj2.PH, obj2.coef * obj1, obj2.pStr);
                end
            elseif(isa(obj1,'pauliString') && isa(obj2,'pauliString'))
                % Store temporarily the product of one term in obj1 times
                % complete obj2 in array of pauliString objects and then sum
                pss = [];
                % Loop first object
                for i = 1 : obj1.len
                    tmpStr = cell(0, 1);
                    tmpCoef = [];
                    % Loop second object
                    for j = 1 : obj2.len
                        % Obtain resulting Pauli string and its phase
                        [ph, chs] = obj2.intnMult(obj1.pStr{i}, obj2.pStr{j});
                        tmpStr{end + 1, 1} = chs;
                        tmpCoef = [tmpCoef; ph * obj1.coef(i) * obj2.coef(j)];
                    end
                    pss = [pss; pauliString(obj1.PH, tmpCoef, tmpStr)];
                end
                % Sum all temporary terms
                m = pss(1);
                for i = 2 : obj1.len
                    m = m + pss(i);
                end
            else
                error("INTERNAL ERROR: Unknown overload of * for pauñiString");
            end
        end % end mtimes
        % Real part of a Pauli String
        function obj2 = realPS(obj1)
            if(obj1.len == 0)
                obj2 = obj1;
            else
                obj2 = pauliString(obj1.PH, real(obj1.coef), obj1.pStr);
            end
        end
        % Imaginary part of a Pauli String
        function obj2 = imagPS(obj1)
            if(obj1.len == 0)
                obj2 = obj1;
            else            
                obj2 = pauliString(obj1.PH, imag(obj1.coef) * 1j, obj1.pStr);
            end
        end
        % Commutator of two Pauli strings
        function objC = commutator(obj1, obj2)
            objC = obj1 * obj2 - obj2 * obj1;
        end
        % Anticommutator of two Pauli strings
        function objA = anticommutator(obj1, obj2)
            objA = obj1 * obj2 + obj2 * obj1;
        end
        % Conjugate of a Pauli string
        function objD = dagger(obj1)
            objD = obj1;
            for i = 1 : objD.len
                objD.coef(i) = obj1.coef(i)';
            end
        end
        % Define square of a Pauli string
        function sq = mpower(obj1, n)
            if(n ~= 2)
                error("Only square powers are possible for pauliString objects");
            end
            % If performance issues this function can be optimized
            % checking the commuting or anticommuting Pauli strings
            sq = obj1 * obj1;
        end % end mpower
        % Multiply Pauli string by a scalar
        function t = times(a, obj1)
            t = pauliString(obj1.PH, a * obj1.coef, obj1.pStr);
        end % end times
        % Substract two Pauli strings
        function m = minus(obj1, obj2)
            m = obj1 + (-1) .* obj2;
        end % end minus
        % Display pauliString in console
        function disp(o)
            [size_o_n, size_o_m] = size(o);
            if(size_o_n > 1 || size_o_m > 1)
                fprintf("%d x %d array of pauliStrings\n", size_o_n, size_o_m);
                return;
            end
            if(isempty(o) || ~(o.len > 0))
                fprintf("[Empty pauliString]\n");
            else
                for i = 1 : o.len - 1;
                    if(real(o.coef(i)) == 0)
                        fprintf("(1j*%d) * %s + ", imag(o.coef(i)), o.pStr{i});
                    elseif(imag(o.coef(i)) == 0)
                        fprintf("(%d) * %s + ", real(o.coef(i)), o.pStr{i}); 
                    else%(real(o.coef(i)) == 0)
                        fprintf("(%d + 1j*%d) * %s + ", real(o.coef(i)), imag(o.coef(i)), o.pStr{i});
                    end
                end
                if(real(o.coef(o.len)) == 0)
                    fprintf("(1j*%d) * %s\n", imag(o.coef(o.len)), o.pStr{o.len});
                elseif(imag(o.coef(o.len)) == 0)
                    fprintf("(%d) * %s\n", real(o.coef(o.len)), o.pStr{o.len}); 
                else%(real(o.coef(o.len)) == 0)
                    fprintf("(%d + 1j*%d) * %s\n", real(o.coef(o.len)), imag(o.coef(o.len)), o.pStr{o.len});
                end
            end
        end
        function flag = isemptyorzero(obj1)
            if(isempty(obj1) || obj1.len == 0)
                flag = true;
            else
                flag = false;
            end
        end
        % Display pauliString in console
        function s = string(o)
            s = "";
            if(o.len == 0)
                s = s + sprintf("[Empty pauliString]\n");
            else
                for i = 1 : o.len - 1;
                    if(real(o.coef(i)) == 0)
                        s = s + sprintf("(1j*%d) * %s + ", imag(o.coef(i)), o.pStr{i});
                    elseif(imag(o.coef(i)) == 0)
                        s = s + sprintf("(%d) * %s + ", real(o.coef(i)), o.pStr{i}); 
                    else(real(o.coef(i)) == 0)
                        s = s + sprintf("(%d + 1j*%d) * %s + ", real(o.coef(i)), imag(o.coef(i)), o.pStr{i});
                    end
                end
                if(real(o.coef(o.len)) == 0)
                    s = s + sprintf("(1j*%d) * %s\n", imag(o.coef(o.len)), o.pStr{o.len});
                elseif(imag(o.coef(o.len)) == 0)
                    s = s + sprintf("(%d) * %s\n", real(o.coef(o.len)), o.pStr{o.len}); 
                else%(real(o.coef(o.len)) == 0)
                    s = s +sprintf("(%d + 1j*%d) * %s\n", real(o.coef(o.len)), imag(o.coef(o.len)), o.pStr{o.len});
                end
            end
        end        
        % Expected value calculation for a column state vector sv: sv' * pauliString * sv
        % or a weighted state vector 'matrix' sv: trace(sv' * pauliString * sv)
        % in which case each column's square norm in sv represents column's
        % probability
        % eValLS expectation value + measurement count increase for linear
        % system generation
        function e = eValLS(o, sv)
            [n, m] = size(sv);
            if(n ~= 2 ^ o.PH.nqbs)
                error("Invalid state vector in argument");
            end
            % If empty pauliString return 0 (no measurement needed)
            if(o.len == 0)
                e = 0;
                return;
            end
            if(o.updateNeeded)
                error("Internal error. Pauli string lacks map");
            end
            % Expected value is calculated for term Pauli operator as sv' * coef(term) mapT(term) * sv(mapIdx')
            % or reused form measurement dictionary if available (fresh)
            e = 0; % Initialize only once for the whole of the Pauli string
            % Loop on string terms
            for trm = 1 : o.len
                % Check if measurement is available
                stumea = o.PH.DICOM(o.pStr{trm});
                % If there is a fresh measurement available use it
                if(stumea.fresh)
                    e = e + o.coef(trm) * stumea.value;
                % else if measurement is not available, measure and update dictionary
                else % fresh = false
                    % Get map data structure from dictionary and measure
                    stu = o.PH.DIC(o.pStr{trm});
                    % Loop on weighted vectors
                    newmeasure = 0;
                    for j = 1 : m
                        svj = sv(:, j);
                        newmeasure = newmeasure + (svj' .* (double(stu.TR) + 1j * double(stu.TI)) * svj(stu.Idx'));
                        o.PH = o.PH.incMeasuresLS();
                    end     
                    e = e + o.coef(trm) * newmeasure;
                    stumea1.fresh = true;
                    stumea1.value = newmeasure;
                    o.PH.DICOM(o.pStr{trm}) = stumea1;
                    %DEBUG%
                    %fprintf("%s\n", o.pStr{trm});
                end % end if fresh
            end
        end % end eValLS
        % eValCO : idem as eValLS but increasers container measurement count
        function e = eValCO(o, sv)
            [n, m] = size(sv);
            if(n ~= 2 ^ o.PH.nqbs)
                error("Invalid state vector in argument");
            end
            % If empty pauliString return 0 (no measurement needed)
            if(o.len == 0)
                e = 0;
                return;
            end
            if(o.updateNeeded)
                error("Internal error. Pauli string lacks map");
            end
            % Expected value is calculated for term Pauli operator as sv' * coef(term) mapT(term) * sv(mapIdx')
            % or reused form measurement dictionary if available (fresh)
            e = 0; % Initialize only once for the whole of the Pauli string
            % Loop on string terms
            for trm = 1 : o.len
                % Check if measurement is available
                stumea = o.PH.DICOM(o.pStr{trm});
                % If there is a fresh measurement available use it
                if(stumea.fresh)
                    e = e + o.coef(trm) * stumea.value;
                % If measurement is not available, measure and update dictionary
                else % fresh = false
                    % Get map data structure from dictionary and measure
                    stu = o.PH.DIC(o.pStr{trm});
                    % Loop on weighted vectors
                    newmeasure = 0;
                    for j = 1 : m
                        svj = sv(:, j);
                        newmeasure = newmeasure + (svj' .* (double(stu.TR) + 1j * double(stu.TI)) * svj(stu.Idx'));
                        e = e + o.coef(trm) * newmeasure;
                        o.PH = o.PH.incMeasuresCO();
                    end     
                    stumea1.fresh = true;
                    stumea1.value = newmeasure;
                    o.PH.DICOM(o.pStr{trm}) = stumea1;
                end % end if fresh
            end
        end % end eValCO
        % eValAZ : idem as eValLS but increasers analyzeh measurement count
        % and doesn't reuse measurements as the statevector is different
        function e = eValAZ(o, sv)
            [n, m] = size(sv);
            if(n ~= 2 ^ o.PH.nqbs)
                error("Invalid state vector in argument");
            end
            % If empty pauliString return 0 (no measurement needed)
            if(o.len == 0)
                e = 0;
                return;
            end
            if(o.updateNeeded)
                error("Internal error. Pauli string lacks map");
            end
            % Expected value is calculated for term Pauli operator as sv' * coef(term) mapT(term) * sv(mapIdx')
            e = 0;
            % Loop on string terms
            for trm = 1 : o.len
                % Get map data structure from dictionary
                stu = o.PH.DIC(o.pStr{trm});
                % Loop on weighted vectors
                for j = 1 : m
                    svj = sv(:, j);
                    e = e + o.coef(trm) * (svj' .* (double(stu.TR) + 1j * double(stu.TI)) * svj(stu.Idx'));
                    o.PH = o.PH.incMeasuresAZ();
                end
            end
        end % end eValAZ        
        % Map calculation. WARNING the memory allocation below only admits
        % up to 15 qubits maximum
        function o = calcMap(o)
            % Allocate space for index and coefficients T
            % Allocate structure for Pauli operator defintion
            stu.Idx = zeros(1, o.PH.twoN, 'uint16');
            stu.TR = zeros(1, o.PH.twoN, 'int8'); % Pauli strings only have coefficients +1, -1, +1j or -1j
            stu.TI = zeros(1, o.PH.twoN, 'int8'); % Pauli strings only have coefficients +1, -1, +1j or -1j
            % Loop terms in Pauli string
            for trm = 1 : o.len
                if(~isKey(o.PH.DIC, o.pStr{trm}))
                    %fprintf("Calculating %s\n", o.pStr{trm});
                    % Loop rows from 1 to 2 ^ nqbs
                    for row = 1 : o.PH.twoN
                        % Calculate the "tensor" indices for this row
                        % Example, for 2 qubits a_ij * b_kl row 1 has tensor indices i = 1, k = 1
                        % Then deduce the only non zero element of the Pauli tensor product (j, l)
                        % And calculate the resulting index as (j - 1) * 2 + l
                        % Calculate sequence number in base 2
                        str = dec2base(row - 1, 2);
                        % Pad zeroes on the left up to no. of qubits
                        lenStr = length(str);
                        for c = (lenStr + 1) : o.PH.nqbs
                            str = strcat('0', str); 
                        end
                        % Retrieve indices in row vector
                        tensorIdxRow = zeros(1, o.PH.nqbs);
                        tensorIdxCol = zeros(1, o.PH.nqbs);
                        for q = 1 : o.PH.nqbs
                            tensorIdxRow(q) = 1 + str2num(str(q));
                        end
                        % Now deduce the non zero elements to infere column and accompanying coefficient
                        % Initialize coefficient
                        t = 1;
                        for q = 1 : o.PH.nqbs
                            % If operator is I or Z (diagonal) the column index is the same
                            if(o.pStr{trm}(q) == 'I')
                                tensorIdxCol(q) = tensorIdxRow(q);
                                t = t * o.PH.I(tensorIdxRow(q), tensorIdxCol(q));
                            elseif(o.pStr{trm}(q) == 'Z')
                                tensorIdxCol(q) = tensorIdxRow(q);
                                t = t * o.PH.Z(tensorIdxRow(q), tensorIdxCol(q));
                            % Otherwise take opposite index for X and Y
                            elseif(o.pStr{trm}(q) == 'X')
                                % Take opposite index (offset = 1)
                                tensorIdxCol(q) = not(tensorIdxRow(q) - 1) + 1;
                                t = t * o.PH.X(tensorIdxRow(q), tensorIdxCol(q));
                            elseif(o.pStr{trm}(q) == 'Y')
                                tensorIdxCol(q) = not(tensorIdxRow(q) - 1) + 1;
                                t = t * o.PH.Y(tensorIdxRow(q), tensorIdxCol(q));
                            end
                        end % end for q
                        % Now calculate column index 
                        col = tensorIdxCol(o.PH.nqbs);
                        for q = 1 : o.PH.nqbs - 1
                            col = col + (tensorIdxCol(q) - 1) * 2 ^ (o.PH.nqbs - q);
                        end
                        stu.Idx(row) = uint16(col);
                        % Coefficient storage in 1-byte arrays (real and imaginary parts separated)
                        if(imag(t) == 0)
                            stu.TR(row) = int8(real(t));
                        elseif(real(t) == 0)
                            stu.TI(row) = int8(imag(t));
                        else
                            error("Internal error. Complex number not allowed");
                        end
                    end % end for row
                    % Update dictionary through handle
                    o.PH.DIC(o.pStr{trm}) = stu;
                    % Create entry in measurement dictionary
                    stumea.fresh = false; % fresh will be true when first measurement is stored in value
                    stumea.value = 0;
                    o.PH.DICOM(o.pStr{trm}) = stumea;
                end % end if isKey
            end % end for trm
            o.updateNeeded = false;
        end % end calcMap
        % String to matrix translation
        function M = strToMatrix(o, inputStr)
            if(length(inputStr) ~= o.PH.nqbs)
                error("Internal error: string must have length nqbs");
            end
            % Allocate space for the Pauli operators in matrix form
            CMults = cell(1, o.PH.nqbs);
            for q = 1 : o.PH.nqbs
                if(inputStr(q) == 'I')
                    PauliIndex = 1;
                elseif(inputStr(q) == 'X')
                    PauliIndex = 2;
                elseif(inputStr(q) == 'Y')
                    PauliIndex = 3;
                elseif(inputStr(q) == 'Z')
                    PauliIndex = 4;      
                else
                    error("Internal error: string contains characters different from I, X, Y, Z");
                end
                CMults{q} = o.PH.PAULI{PauliIndex};
            end
            M = o.PH.kp(CMults);
        end % end strToMatrix        
        % Matrix representation of the total Pauli operator in object with coefficients
        function M = matrixRep(o)
            M = 0;
            for i = 1 : o.len
                M = M + o.coef(i) * o.strToMatrix(o.pStr{i});
            end
        end
    end % end methods
end