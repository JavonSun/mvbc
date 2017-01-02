% evaluate on simulation data

% Javon, 2/3/2015
% Modified by Tingyang, 5/18/2015

clear;
clc;
%% Data
% Change different e here to run different simulation datasets. e can be 1,
% 0.8, 0.6, or 0.4
e=1;
switch e
    case 1
        % e=1 % NMI=0.6237
        folder='genotype2-es1';
        fprintf('Evaluating on simulation data with e=1...\n');
        % These are the parameters used for the first cluster. The values
        % of the parameters is determined by the PCA
        sz1 = 320;
        sv1 = [3; 10];
        ini_v1 = 1;
        % For the second cluster
        sz2 = 205;
        sv2 = [3; 10];
        ini_v2 = 4;
    case .8
        % e=0.8 %NMI=0.6226
        folder='genotype2-es.8';
        fprintf('Evaluating on simulation data with e=0.8...\n');
        sz1 = 328;
        sv1 = [3; 10];
        ini_v1 = 1;
        sz2 = 194;
        sv2 = [3; 10];
        ini_v2 = 1;
    case .6
        % e=0.6 %NMI=0.6125
        folder='genotype2-es.6';
        fprintf('Evaluating on simulation data with e=0.6...\n');
        sz1 = 330;
        sv1 = [3; 10];
        ini_v1 = 1;
        sz2 = 200;
        sv2 = [3; 10];
        ini_v2 = 4;
    case .4
        % e=0.4 %NMI=0.6099
        folder='genotype2-es.4';
        fprintf('Evaluating on simulation data with e=0.4...\n');
        sz1 = 328;
        sv1 = [3; 10];
        ini_v1 = 1;
        sz2 = 220;
        sv2 = [3; 10];
        ini_v2 = 4;
end
phe_data = dataset('XLSFile', ['..\sample-data\phe-gen\' folder '\phe_data.csv']);
gen_data = dataset('XLSFile', ['..\sample-data\phe-gen\' folder '\gen_data.csv']);
subtype_assign = dataset('XLSFile', ['..\sample-data\phe-gen\' folder '\subtype_assign.csv']);

n = size(phe_data, 1);
pd = size(phe_data, 2);
gd = size(gen_data, 2);
lbl = phe_data.lbl;

gen_lbl = zeros(size(subtype_assign, 1), 1);
gen_lbl(subtype_assign.gen_sub1 == 1) = 1;
gen_lbl(subtype_assign.gen_sub2 == 1) = 2;
nmi_views = nmi(lbl, gen_lbl);

M_phe_raw = double(phe_data(:, 1:end-1));
M_gen_raw = double(gen_data(:, 5:end))';
% normalize
M_phe_norm = normc(M_phe_raw);
M_gen_norm = normc(M_gen_raw);

M_phe = M_phe_norm;
M_gen = M_gen_norm;

M = cell(1, 1);
M{1} = M_phe;
M{2} = M_gen;

%% Proposed Method
[z1, U1, V1, obj1] = mvlrrl0(M, sz1, sv1, ini_v1);
fprintf('\nThe distribution of the true labels of the first identified cluster.\n');
results1=tbl(lbl(z1 ~= 0));
disp([[{'True Label: '};{'Counts: '}] num2cell(results1)]);

M2_dv = cell(2, 1);
M2_dv{1} = M_phe(z1 == 0, :);
M2_dv{2} = M_gen(z1 == 0, :);
lbl2 = lbl(z1 == 0);

[z2, U2, V2, obj2] = mvlrrl0(M2_dv, sz2, sv2, ini_v2);
fprintf('\nThe distribution of the true labels of the second identified cluster.\n');
results2=tbl(lbl(z2 ~= 0));
disp([[{'True Label: '};{'Counts: '}] num2cell(results2)]);

lbl_dv = double(z1 ~= 0);
lbl_dv(z1 == 0) = double(z2 ~= 0);
lbl_dv(z1 == 0 & lbl_dv == 1) = 2;
nmi_dv = nmi(lbl, lbl_dv); 

%% Save the results
rs = struct();
rs.sz1 = sz1;
rs.sv1 = sv1;
rs.z1 = z1;
rs.U1 = U1;
rs.V1 = V1;
rs.sz2 = sz2;
rs.sv2 = sv2;
rs.z2 = z2;
rs.U2 = U2;
rs.V2 = V2;
rs.lbl_dv = lbl_dv;
rs.nmi_dv = nmi_dv;
% save(sprintf('rs-es%g.mat',e), 'rs');

%% print the NMI
fprintf('NMI: %.4f\n',rs.nmi_dv);







