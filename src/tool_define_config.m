current_path = pwd();
config_path = current_path + "/params";

i = 0;

% 1. Build an array of hamiltonian definitions
% 4 QUBIT MODELS
% H2 model, STO-3G basis
i = i + 1;
cfg(i).model = "H2";
cfg(i).type = "fermions";
cfg(i).hamiltonian = "/Q04_H2_HHF00.m";
cfg(i).ansatz = "/Q04_UCCSD0011_AZF00.m";
cfg(i).domain = "/Q04_H2_DDF00.m";
cfg(i).parameter = "/Q04_H2_PPF00.m";

% Ising model, open boundaries, domain 3
i = i + 1;
cfg(i).model = "ISIOB04";
cfg(i).type = "spins";
cfg(i).hamiltonian = "/Q04_ISIOB_HHS00.m";
cfg(i).ansatz = "/Q04_ISIOB_AZS00.m";
cfg(i).domain = "/Q04_ISIOB_DDS00.m";
cfg(i).parameter = "/Q04_ISIOB_PPS00.m";

% 6 QUBIT MODELS
% Ising model, open boundaries, domain 4
i = i + 1;
cfg(i).model = "ISIOB06";
cfg(i).type = "spins";
cfg(i).hamiltonian = "/Q06_ISIOB_HHS00.m";
cfg(i).ansatz = "/Q06_ISIOB_AZS00.m";
cfg(i).domain = "/Q06_ISIOB_DDS00.m";
cfg(i).parameter = "/Q06_ISIOB_PPS00.m";

% Ising model (symmetric), open boundaries, domain 4
i = i + 1;
cfg(i).model = "ISYOB06";
cfg(i).type = "spins";
cfg(i).hamiltonian = "/Q06_ISYOB_HHS00.m";
cfg(i).ansatz = "/Q06_ISYOB_AZS00.m";
cfg(i).domain = "/Q06_ISYOB_DDS00.m";
cfg(i).parameter = "/Q06_ISYOB_PPS00.m";

% Heisenberg model (symmetric), open boundaries, domain 4
i = i + 1;
cfg(i).model = "HEYOB06";
cfg(i).type = "spins";
cfg(i).hamiltonian = "/Q06_HEYOB_HHS00.m";
cfg(i).ansatz = "/Q06_HEYOB_AZS00.m";
cfg(i).domain = "/Q06_HEYOB_DDS00.m";
cfg(i).parameter = "/Q06_HEYOB_PPS00.m";

% Hubbard model, open boundaries, domain 4
i = i + 1;
cfg(i).model = "HUBOB06";
cfg(i).type = "spins";
cfg(i).hamiltonian = "/Q06_HUBOB_HHS00.m";
cfg(i).ansatz = "/Q06_HUBOB_AZS00.m";
cfg(i).domain = "/Q06_HUBOB_DDS00.m";
cfg(i).parameter = "/Q06_HUBOB_PPS00.m";

% 8 QUBIT MODELS
% Heisenberg model (symmetric), open boundaries, domain 5
i = i + 1;
cfg(i).model = "HEYOB08";
cfg(i).type = "spins";
cfg(i).hamiltonian = "/Q08_HEYOB_HHS00.m";
cfg(i).ansatz = "/Q08_HEYOB_AZS00.m";
cfg(i).domain = "/Q08_HEYOB_DDS00.m";
cfg(i).parameter = "/Q08_HEYOB_PPS00.m";

% Hydrogen chain, 4 sites, interatomic distance 0.60 Angstrom
i = i + 1;
cfg(i).model = "H4CHAIN060";
cfg(i).type = "fermions";
cfg(i).hamiltonian = "/Q08_H4CHAIN060_HHF00.m";
cfg(i).ansatz = "/Q08_UCCGSD_AZF00.m";
cfg(i).domain = string.empty; % Whole domain by default
cfg(i).parameter = string.empty; % By default 1 

% Hydrogen chain, 4 sites, interatomic distance 0.65 Angstrom
i = i + 1;
cfg(i).model = "H4CHAIN065";
cfg(i).type = "fermions";
cfg(i).hamiltonian = "/Q08_H4CHAIN065_HHF00.m";
cfg(i).ansatz = "/Q08_UCCGSD_AZF00.m";
cfg(i).domain = string.empty; % Whole domain by default
cfg(i).parameter = string.empty; % By default 1 

% Hydrogen chain, 4 sites, interatomic distance 0.70 Angstrom
i = i + 1;
cfg(i).model = "H4CHAIN070";
cfg(i).type = "fermions";
cfg(i).hamiltonian = "/Q08_H4CHAIN070_HHF00.m";
cfg(i).ansatz = "/Q08_UCCGSD_AZF00.m";
cfg(i).domain = string.empty; % Whole domain by default
cfg(i).parameter = string.empty; % By default 1 

% Hydrogen chain, 4 sites, interatomic distance 0.75 Angstrom
i = i + 1;
cfg(i).model = "H4CHAIN075";
cfg(i).type = "fermions";
cfg(i).hamiltonian = "/Q08_H4CHAIN075_HHF00.m";
cfg(i).ansatz = "/Q08_UCCGSD_AZF00.m";
cfg(i).domain  = string.empty; % Whole domain by default
cfg(i).parameter = string.empty; % By default 1 

% Hydrogen chain, 4 sites, interatomic distance 0.80 Angstrom
i = i + 1;
cfg(i).model = "H4CHAIN080";
cfg(i).type = "fermions";
cfg(i).hamiltonian = "/Q08_H4CHAIN080_HHF00.m";
cfg(i).ansatz = "/Q08_UCCGSD_AZF00.m";
cfg(i).domain = string.empty; % Whole domain by default
cfg(i).parameter = string.empty; % By default 1 

% Hydrogen chain, 4 sites, interatomic distance 0.85 Angstrom
i = i + 1;
cfg(i).model = "H4CHAIN085";
cfg(i).type = "fermions";
cfg(i).hamiltonian = "/Q08_H4CHAIN085_HHF00.m";
cfg(i).ansatz = "/Q08_UCCGSD_AZF00.m";
cfg(i).domain = string.empty; % Whole domain by default
cfg(i).parameter = string.empty; % By default 1 

% Hydrogen chain, 4 sites, interatomic distance 0.90 Angstrom
i = i + 1;
cfg(i).model = "H4CHAIN090";
cfg(i).type = "fermions";
cfg(i).hamiltonian = "/Q08_H4CHAIN090_HHF00.m";
cfg(i).ansatz = "/Q08_UCCGSD_AZF00.m";
cfg(i).domain = string.empty; % Whole domain by default
cfg(i).parameter = string.empty; % By default 1 

% Hydrogen chain, 4 sites, interatomic distance 0.95 Angstrom
i = i + 1;
cfg(i).model = "H4CHAIN095";
cfg(i).type = "fermions";
cfg(i).hamiltonian = "/Q08_H4CHAIN095_HHF00.m";
cfg(i).ansatz = "/Q08_UCCGSD_AZF00.m";
cfg(i).domain = string.empty; % Whole domain by default
cfg(i).parameter = string.empty; % By default 1 

% Hydrogen chain, 4 sites, interatomic distance 1.00 Angstrom
i = i + 1;
cfg(i).model = "H4CHAIN100";
cfg(i).type = "fermions";
cfg(i).hamiltonian = "/Q08_H4CHAIN100_HHF00.m";
cfg(i).ansatz = "/Q08_UCCGSD_AZF00.m";
cfg(i).domain = string.empty; % Whole domain by default
cfg(i).parameter = string.empty; % By default 1 

% Hydrogen chain, 4 sites, interatomic distance 1.05 Angstrom
i = i + 1;
cfg(i).model = "H4CHAIN105";
cfg(i).type = "fermions";
cfg(i).hamiltonian = "/Q08_H4CHAIN105_HHF00.m";
cfg(i).ansatz = "/Q08_UCCGSD_AZF00.m";
cfg(i).domain = string.empty; % Whole domain by default
cfg(i).parameter = string.empty; % By default 1 

% Hydrogen chain, 4 sites, interatomic distance 1.10 Angstrom
i = i + 1;
cfg(i).model = "H4CHAIN110";
cfg(i).type = "fermions";
cfg(i).hamiltonian = "/Q08_H4CHAIN110_HHF00.m";
cfg(i).ansatz = "/Q08_UCCGSD_AZF00.m";
cfg(i).domain = string.empty; % Whole domain by default
cfg(i).parameter = string.empty; % By default 1 

% Hydrogen chain, 4 sites, interatomic distance 1.15 Angstrom
i = i + 1;
cfg(i).model = "H4CHAIN115";
cfg(i).type = "fermions";
cfg(i).hamiltonian = "/Q08_H4CHAIN115_HHF00.m";
cfg(i).ansatz = "/Q08_UCCGSD_AZF00.m";
cfg(i).domain = string.empty; % Whole domain by default
cfg(i).parameter = string.empty; % By default 1 

% Hydrogen chain, 4 sites, interatomic distance 1.20 Angstrom
i = i + 1;
cfg(i).model = "H4CHAIN120";
cfg(i).type = "fermions";
cfg(i).hamiltonian = "/Q08_H4CHAIN120_HHF00.m";
cfg(i).ansatz = "/Q08_UCCGSD_AZF00.m";
cfg(i).domain = string.empty; % Whole domain by default
cfg(i).parameter = string.empty; % By default 1 

% 2. Encode as JSON (pretty form for readability is optional)
jsonText = jsonencode(cfg, "PrettyPrint", true);

% 3. Write to disk
fid = fopen(config_path + "/cfg.json",'w');
fwrite(fid, jsonText, 'char');
fclose(fid);
