% simulation study

% Javon, 12/30/2016

clear;
clc;

%% load data
M = cell(2, 1);
M{1} = dlmread('../sample-data/gene-expr/view1-data.csv', ',');
M{2} = dlmread('../sample-data/gene-expr/view2-data.csv', ',');
trueClus = dlmread('../sample-data/gene-expr/true-clus.csv', ',');

%% multi-view bi-clustering
rowClus = zeros(size(M{1}, 1), 1); % row cluster shared accross datasets
colClus = cell(length(M), 1);
dimRatio = zeros(length(M), 1);
for i = 1:length(M)
    colClus{i} = [];
    dimRatio(i) = size(M{i}, 1) / size(M{i}, 2);
end

% first cluster
lz = 1.2;
lamb1_chose = 1.2;
lamb2_chose = 2;
iClus = 1;
[rowClus1, colClus1, z1, U, V, obj1] = mvslra(M, lz, lamb1_chose * ones(length(M), 1), ...
    dimRatio * lamb2_chose, size(M{1}, 1), 'weighted');
fprintf('none zero in u1: %d, none zero in u2: %d, shared none zero: %d\n', ...
    length(find(U(:, 1) ~= 0)), length(find(U(:, 2) ~= 0)), ...
    length(find(U(:, 1) ~= 0 & U(:, 2) ~= 0)));
% none zero in u1: 298, none zero in u2: 228, shared none zero: 228
% rowClusTmp = zeros(size(M{1}, 1), 1);
% rowClusTmp(U(:, 1) ~= 0 & U(:, 2) ~= 0) = iClus;
% sampleInClus = find(rowClusTmp == iClus);
% sampleInClus'
% for i = 1:length(M)
%     [colClus{i} V{i}]
% end

% assign class 
rowClus(U(:, 1) ~= 0 & U(:, 2) ~= 0) = iClus;
for i = 1:length(M)
    colClus{i} = [colClus{i} V{i}];
end

% second cluster
iClus = 2;
M_r = cell(length(M), 1);
for i = 1:length(M)
    M_r{i} = M{i}(rowClus == 0, :);
end

lz = 1;
lamb1_chose = 1.3;
lamb2_chose = 1.4;

[rowClus2, colClus2, z2, U, V, obj2] = mvslra(M_r, lz, lamb1_chose * ones(length(M_r{1}), 1), ...
    dimRatio * lamb2_chose, size(M_r{1}, 1), 'weighted');
fprintf('none zero in u1: %d, none zero in u2: %d, shared none zero: %d\n', ...
    length(find(U(:, 1) ~= 0)), length(find(U(:, 2) ~= 0)), ...
    length(find(U(:, 1) ~= 0 & U(:, 2) ~= 0)));
% none zero in u1: 202, none zero in u2: 322, shared none zero: 201
% rowClusTmp = rowClus;
tmp = zeros(size(M_r{1}, 1), 1);
tmp(U(:, 1) ~= 0 & U(:, 2) ~= 0) = iClus;
% rowClusTmp(rowClus == 0) = tmp;
% sampleInClus = find(rowClusTmp == iClus);
% sampleInClus'
% for i = 1:length(M_r)
%     [colClus{i} V{i}]
% end

% assign class 
rowClus(rowClus == 0) = tmp;
for i = 1:length(M_r)
    colClus{i} = [colClus{i} V{i}];
end

% third cluster
iClus = 3;
M_r = cell(length(M), 1);
for i = 1:length(M)
    M_r{i} = M{i}(rowClus == 0, :);
end

lz = 0.9;
lamb1_chose = 1.2;
lamb2_chose = 1.3;

[rowClus3, colClus3, z3, U, V, obj3] = mvslra(M_r, lz, lamb1_chose * ones(length(M_r{1}), 1), ...
    dimRatio * lamb2_chose, size(M_r{1}, 1), 'weighted');
fprintf('none zero in u1: %d, none zero in u2: %d, shared none zero: %d\n', ...
    length(find(U(:, 1) ~= 0)), length(find(U(:, 2) ~= 0)), ...
    length(find(U(:, 1) ~= 0 & U(:, 2) ~= 0)));
% none zero in u1: 219, none zero in u2: 233, shared none zero: 210
% rowClusTmp = rowClus;
tmp = zeros(size(M_r{1}, 1), 1);
tmp(U(:, 1) ~= 0 & U(:, 2) ~= 0) = iClus;
% rowClusTmp(rowClus == 0) = tmp;
% sampleInClus = find(rowClusTmp == iClus);
% sampleInClus'
% for i = 1:length(M_r)
%     [colClus{i} V{i}]
% end

% assign class 
rowClus(rowClus == 0) = tmp;
for i = 1:length(M_r)
    colClus{i} = [colClus{i} V{i}];
end

% calulate nmi
nmi_mvbi = nmi(trueClus, rowClus);

%% print the NMI
fprintf('NMI: %.4f\n',nmi_mvbi);




