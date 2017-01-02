% evaluate on simulation data

% Javon, 12/30/2016

clear;
clc;
%% Data
% Change different e here to run different simulation datasets. e can be 1,
% 0.8, 0.6, or 0.4
e=0.4;
switch e
    case 1
        % e=1 % NMI=0.6237
        folder='genotype2-es1';
        fprintf('Evaluating on simulation data with e=1...\n');
    case .8
        % e=0.8 %NMI=0.6226
        folder='genotype2-es.8';
        fprintf('Evaluating on simulation data with e=0.8...\n');
    case .6
        % e=0.6 %NMI=0.6125
        folder='genotype2-es.6';
        fprintf('Evaluating on simulation data with e=0.6...\n');
    case .4
        % e=0.4 %NMI=0.6099
        folder='genotype2-es.4';
        fprintf('Evaluating on simulation data with e=0.4...\n');
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
% set parameters
lz1 = 0.01;
lv1 = [0.45; 1];
lz2 = 0.018;
lv2 = [0.5; 0.7];
% getting first cluster
[rowClus1, colClus1, z1, U1, V1, obj1] = mvssvd(M, lz1, lv1);
fprintf('\nThe distribution of the true labels of the first identified cluster.\n');
results1=tbl(lbl(rowClus1 ~= 0));
disp([[{'True Label: '};{'Counts: '}] num2cell(results1)]);

% getting second cluster
M2_dv = cell(2, 1);
M2_dv{1} = M_phe(rowClus1 == 0, :);
M2_dv{2} = M_gen(rowClus1 == 0, :);
lbl2 = lbl(rowClus1 == 0);

[rowClus2, colClus2, z2, U2, V2, obj2] = mvssvd(M2_dv, lz2, lv2);
fprintf('\nThe distribution of the true labels of the second identified cluster.\n');
results2=tbl(lbl2(rowClus2 ~= 0));
disp([[{'True Label: '};{'Counts: '}] num2cell(results2)]);

% generate final cluster assignment
lbl_dv = rowClus1;
lbl_dv(rowClus1 == 0) = rowClus2;
lbl_dv(rowClus1 == 0 & lbl_dv == 1) = 2;
nmi_dv = nmi(lbl, lbl_dv);

%% Save the results
rs = struct();
rs.lbl_dv = lbl_dv;
rs.nmi_dv = nmi_dv;
% save(sprintf('rs-es%g.mat',e), 'rs');

%% print the NMI
fprintf('NMI: %.4f\n',rs.nmi_dv);







