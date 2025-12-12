% Class definition hamiltonian : Hardcoded Hamiltonian models and parameters in
% master directory are processed according to user's choice. Class constructor
% permits to subgroup Hamiltonian terms according to row vectors with
% indices passed in the variable arguments (varargin). These indices
% correspond to Hamiltonian term definitions in file _HHF??.m or _HHS??.m
% An ansatz is provided in files _AZF??.m or _AZS??.m. For fermionic
% hamiltonians ansätze are typically a list of antihermitian single or
% double excitations and de-excitations. For spin hamiltonians, ansätze are typically 
% a list of antihermitian Pauli strings (e.g. -i * XXI, -i * XXZ etc.
% The class manages hamiltonian construction in the desired partition
% along with a specific ansatz per partition term. This is defined in
% domain files _DDF??.m or _DDS??.m
% For fermionic hamiltonians, the ansatz is specified by a list of ansatz
% operator numbers (as defined in ansatz file) per hamiltonian term (as
% defined in the hamiltonian file)
% For spin hamiltonians a string mask per hamiltonian term such as 'IIXXII' is passed in the
% domain file, and then strings corresponding to that mask in the ansatz
% file are selected. X stands for any Pauli string, whereas I forces
% strings matching exactly the same Is
% Domains are encoded as follows: ? means that the qubit is relevant, I that it is excluded
% Class manages the exact eigenvalue / eigenvector calculation for reference
classdef hamiltonian < handle
    properties
        EPSEIG = 1.0e-12; % Tolerance to consider that two numerically close eigenvalues are degenerate
        % Handle to Pauli dictionary
        PDH
        % Model type
        type         % "fermions" or "spins"
        % General properties
        model        % Hamiltonian name
        H            % Hamiltonian (array of Pauli strings, one per term in the partition, total M)   
        M            % Number of local hamiltonians composing H (array size) H(1),...,H(M)
        HT           % Total hamiltonian (grouped in one Pauli string)
        cellAZ       % Cell of ansätze (one array of pauliStrings per regrouped Hamiltonian term) AZ(1),...,AZ(M)
    end
    methods
        % Class constructor either in mode specified in config file: "fermions" or = "spins"
        % with varargin regroupments specified as vectors containing distinct idexes
        % Parameter p is model dependant. If not needed use value 1
        % First argument is a handle to the unique Pauli dictionary
        function o = hamiltonian(PDH, model, p, varargin)
            % Total number of terms in partition
            o.M = size(varargin, 2);
            % Handle to pauliDictionary
            o.PDH = PDH;
            % Retrieve configuration info and look up model
            currentPath = pwd;
            configPath = currentPath + "/params";
            masterPath = currentPath + "/master";
            addpath(currentPath);
            addpath(masterPath);
            % Load configuration structure
            stuConfig = jsondecode(fileread(configPath + "/cfg.json")); 

            allModels = string({stuConfig.model});
            % Find index where it matches
            idr = find(allModels == model, 1);

            if(isempty(idr))
                error("Inexistant hamiltonian model in configuration file");
            end
            if(stuConfig(idr).type ~= "fermions" && stuConfig(idr).type ~= "spins")
                error("Unknown model type. It must be 'fermions' or 'spins'");
            end
            if(isempty(stuConfig(idr).hamiltonian))
                error("Missing hamiltonian file definition");
            end
            if(isempty(stuConfig(idr).ansatz))
                error("Missing ansatz file definition");             
            end
            
            % Retrieve model and ansatz definitions from .m files
            run(masterPath + stuConfig(idr).hamiltonian); % Loads array of pauliStrings HH
            run(masterPath + stuConfig(idr).ansatz);  % Loads array of pauliStrings AZ
            % Read model and ansatz sizes as hardcoded in files
            MHH = size(HH, 2);
            MAZ = size(AZ, 2);

            % Check that hamiltonian terms are hermitian
            for t = 1 : MHH
                if(~isemptyorzero(HH(t) - dagger(HH(t))))
                    HH(t)
                    error("Non hermitian terms found in hamiltonian file definition");
                end
            end

            % Check that ansatz terms are anti-hermitian
            for t = 1 : MAZ
                if(~isemptyorzero(AZ(t) + dagger(AZ(t))))
                    error("Non anti-hermitian terms found in ansatz file definition");
                end
            end

            % Retrieve model parameter info if existant, else p = 1 everywhere
            if(~isempty(stuConfig(idr).parameter))
                run(masterPath + stuConfig(idr).parameter); % PP
                % Retrieve size
                MPP = size(PP, 2);
                if(MPP ~= MHH)
                    error("Parameter vector length inconsistent with hamiltonian model");
                end
                % Propagate parameter p in hamiltonian terms
                for i = 1 : MPP
                    if(PP(i) < 0 || PP(i) > 0)
                        HH(1, i) = (p * PP(i)) * HH(1, i);
                    else % If 0, do nothing
                        HH(1, i) = (1) * HH(1, i);
                    end
                end
            end

            % Define hamiltonian terms (one regrouped pauliString per term in partition)
            o.H = pauliString.empty;
            for m = 1 : o.M
                arrH = HH(varargin{m});
                size_arrH = size(arrH, 2);
                o.H(m) = arrH(1);
                for s = 2 : size_arrH
                    o.H(m) = o.H(m) + arrH(s);
                end
            end
            
            o.cellAZ = cell(0);
            % Define ansatz cell
            if(isempty(stuConfig(idr).domain))
                % If no doimains are specified, ansatz AZ is replicated for
                % each hamiltonian term, both for 'fermions' and 'spins'
                for m = 1 : o.M
                    o.cellAZ{end+1} = AZ;
                end
            else
                % Retrieve domain definitions from .m files
                % Beware: DD may be a cell (type == "fermions") or an array of pauliStrings (type == "spins")
                run(masterPath + stuConfig(idr).domain); % DD
                MDD = size(DD, 2);
                if(MDD ~= MHH)
                    error("Domain definitions inconsistent with model hamiltonian in size");
                end
                if(stuConfig(idr).type == "fermions")
                    % Merge domains for each partition term
                    for m = 1 : o.M
                        DDM = DD(varargin{m});
                        size_DDM = size(DDM, 2);
                        % Calculate effective ansatz indices
                        idd = DDM{1};
                        for i = 2 : size_DDM
                            idd = [idd, DDM{i}];
                        end
                        idd = sort(unique(idd));
                        o.cellAZ{end+1} = AZ(idd); % One ansatz per term m
                    end
                elseif(stuConfig(idr).type == "spins")
                    % For spins the ansatz is expected with one single pauliString per term
                    for z = 1 : MAZ
                        if(AZ(z).len ~= 1)
                            error("Incorrect ansatz file. For spins only one-term pauliStrings are expected");
                        end
                    end
                    % Merge domains for each partition term when their mask
                    % matches strings in ansatz
                    for m = 1 : o.M
                        DDM = DD(varargin{m});
                        size_DDM = size(DDM, 2);
                        includeAZterm = [];
                        % Compare every mask in partition to every
                        % pauliString in ansatz
                        for d = 1 : size_DDM
                            for z = 1 : MAZ
                                if(o.matchesI(DDM(d).pStr{1}, AZ(z).pStr{1}))
                                    includeAZterm = [includeAZterm, z];
                                end
                            end
                        end
                        includeAZterm = sort(unique(includeAZterm));
                        o.cellAZ{end+1} = AZ(includeAZterm);
                    end
                end
            end
            o.model = model;
            o.type = stuConfig(idr).type;
            % Total hamiltonian calculation
            o.HT = o.H(1);
            for m = 2 : o.M
                o.HT = o.HT + o.H(m);
            end
        end % constructor
        % Eigenvalues of the problem at hand
        function eigenValues = exact_eigenvalues(o)
            eigenValues = eig(o.HT.matrixRep());
        end
        % Minimum eigenvalue of the problem at hand
        function e = exact_eigenvalue(o)
            e = min(exact_eigenvalues(o));
        end
        % Eigenvector(s) of minimum eigenvalue of the problem at hand
        function eigenVectors = exact_eigenvectors(o)
            % We don't suppose that the smallest eigenvalue is not degenerate
            [Vs, Ds] = eig(o.HT.matrixRep());
            [Ds, Is] = sort(diag(Ds));
            Vs = Vs(:, Is);
            % Possible degenerate eigenstates
            idx = (min(Ds) * (1 - o.EPSEIG) >= Ds);
            eigenVectors = Vs(:, idx');
        end     
        function flag = matchesI(o, template, target)
            % Returns true iff every position that is 'I' in template
            % is also 'I' in target (target may have extra 'I's elsewhere).
            maskI = (template == 'I');
            flag = all(target(maskI) == 'I');
        end
    end % methods
end % class definition