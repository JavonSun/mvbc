function [z, U, V, obj] = mvslra_prox(M, sz, sv, iSeedV1) 

% Problem formulation:
% min (1/2) \sum ||M_i - (z \dot u_i)v_i^T||_F^2
% s.t. ||z||_0 <= sz, ||v_i||_0 <= sv(i)

% Reference
% Multi-view Sparse Co-clustering via Proximal Alternating Linearized Minimization
% Jiangwen Sun, Jin Lu, Tiangyang and Jinbo Bi
% In the Proceedings of The 32nd International Conference on Machine Learning (ICML), 2015

% Inputs:
%   M - a cell array of data matrixes from multiple views. Rows in each
%       matrix represent samples, columns represent features. It is assumed
%       samples are aligned well among all the matrix. In other words, rows with
%       same row index represent exactly the same sample but charaterized
%       from different views. All columns in matrixes are considered as
%       features in the analysis, so do not include sample id in the
%       matrixes.

%   sz - a scalar, hyperparameter, it controls the maximum non-zeros in vector z 

%   sv - a vector, hyperprameter, it constrols the maximum non-zeros in
%       each vector v_i, so it should have a length exact same as input M 

%   iSeedV1 - the index of the seed feature in view 1. This feature will be 
%       the one most likely remained active in view 1 (i.e., with a non-zero
%       component in v_1). This parameter helps the initialization and kick
%       off the iterative searching. Basically v1 is initialized with all zeros 
%       except this feature. 

% Outputs:
%   z - a vector, exactly the z in the formulation
%   U - a matrix, the i-th column is u_i in the formulation
%   V - a cell array, the i-th cell is v_i in the formulation
%   obj - the value of the objective function calculated with returned z,
%       U and V


% Javon, 2/3/2015

debug = false;

maxIter = 1e03;
threshold = 1e-05;

n = size(M{1}, 1);
nView = length(M);
d = zeros(nView, 1);
z = ones(n,1); % initialize z with all ones
U = ones(n, nView); % initialize u with all ones
V = cell(nView, 1);

gamma_z = 1.2;
gamma_u = 1.2;
gamma_v = 1.2;
for iView = 1:nView
    d(iView) = size(M{iView}, 2);
    V{iView} = ones(d(iView), 1); % initialize v with all ones
end

if debug
    fprintf('Initial objective: %2.2e\n', objective(M, z, U, V));
end

% first iteration
% starts from first view
if nargin > 3
    %initialize v1 with all zeros except 1 for the seed feature
    V{1}(1:end) = 0;
    V{1}(iSeedV1) = 1; 
else 
    V{1} = update_v(M{1}, z, U(:, 1), V{1}, gamma_v, sv(1)); % update v
end

U(:, 1) = update_u(M{1}, z, U(:, 1), V{1}, gamma_u); % update u
% U(:, 1) = map(U(:, 1), sz);
% update z
tmp = map(U(:, 1), sz);
z(tmp == 0) = 0;
% for the rest view
if nView > 1
    for iView = 2:nView
        V{iView} = update_v(M{iView}, z, U(:, iView), V{iView}, gamma_v, sv(iView)); % update v
        U(:, iView) = update_u(M{iView}, z, U(:, iView), V{iView}, gamma_u); % update u
    end
    % update z
    z = update_z(M, z, U, V, gamma_z, sz);
end

obj = objective(M, z, U, V);
if debug
    fprintf('Iter 1 -- objective: %2.2e\n', obj);
end

for iIter = 1:maxIter
    pre_obj = obj;
    
    pre_z = z;
%     pre_U = U;
%     pre_V = V;
    % given z, update U, V
    for iView = 1:nView
        V{iView} = update_v(M{iView}, z, U(:, iView), V{iView}, gamma_v, sv(iView)); % update v
        U(:, iView) = update_u(M{iView}, z, U(:, iView), V{iView}, gamma_u); % update u 
    end

    % given U, V, update z    
    z = update_z(M, z, U, V, gamma_z, sz);
    
    % calculate the objective
    obj = objective(M, z, U, V);
    
    if debug
        fprintf('Iter %d -- objective: %2.2e, change in obj: %2.2e, in z: %2.2e\n', ...
            iIter, obj, pre_obj - obj, norm(pre_z - z));
    end
    
    if norm(pre_z - z) < threshold
        break;
    end
end
end

function [obj] = objective(M, z, U, V)
% calculate objective value
obj = 0;
for i = 1:length(M)
    mat =  M{i} - (U(:, i) .* z) * V{i}';
    obj = obj + trace(mat' * mat) ;
end
end

function [u] = update_u(M, z, u, v, gamma) 
% given z, v, update u
grad = ((z .* u) * v' - M) * v .* z;
tmp = v' * v * z .* z;
lip = sqrt(tmp' * tmp);
if lip == 0
    u = zeros(length(u), 1);
    return;
end
u = u - 1 / (gamma * lip) * grad;
end

function [v] = update_v(M, z, u, v, gamma, s)
% given z, u, update v
grad = (z .* u * v' - M)' * (z .* u);
tmp = z .* u;
lip = tmp' * tmp;
if lip == 0
    v = zeros(length(v), 1);
    return;
end
v = v - 1 / (lip * gamma) * grad;
v = map(v, s);
end

function [z] = update_z(M, z, u, v, gamma, s) 
% given u, v, update z
grad = zeros(length(z), 1);
tmp = zeros(length(z), 1);
for i = 1:length(M)
    grad = grad + (z .* u(:, i) * v{i}' - M{i}) * v{i} .* u(:, i);
    tmp = tmp + v{i}' * v{i} * u(:, i) .* u(:, i);
end
lip = sqrt(tmp' * tmp);
if lip == 0 
    z = zeros(length(z), 1);
    return;
end
z = z - 1 / (lip * gamma) * grad;
z = map(z, s);
end

function u = map(u, s)
% map u to its closest that has <= s none zeros 
u_a = abs(u);
[v, I] = sort(u_a, 1, 'descend');
u(I((s + 1):end)) = 0;
end

