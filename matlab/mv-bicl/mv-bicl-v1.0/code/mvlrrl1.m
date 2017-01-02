function [rowClus, colClus, z, U, V, obj] = mvlrrl1(M, lz, lu, lv, maxRowClus, view_weight, optimal_z) 
% Optimization problem:
% \sum_i(||M_i - (u_i \dot z) * v_i^T||_F^2) + Lz * ||z||_1 + \sum_i(Lu_i * ||u_i||_1)
% + \sum_i(Lv_i * ||v_i||_1)
% s.t. ||z||_0 <= k (maxRowClus)

% Reference:
% A cross-species bi-clustering approach to identifying conserved co-regulated genes
% Jiangwen Sun, Zongliang Jiang, Xiuchun Tian and Jinbo Bi
% Bioinformatics, 32 (12), i137-i146, 2016

% Arguments:
%   M - a cell array containing the list of data matrixes from multiple
%       views
%   Lz - lambda in front of z
%   Lu - lambda in front of u
%   Lv - lambda in front of v
%   maxRowClus - upper bound on the size of row clusters, i.e., maximum non
%       zeros in z
%   view_weight: indicating whether all the views are equal weighted. Two
%       options: 'equal', all views are equal weighted; 'weighted', views
%       are not equal weighted, in this case, the first view is assumed to
%       have the highest weight and will be used to initialize the
%       searching process
%   optimal_z: if known (for example in simulation), used
%       to calculate the optimal objective value. This is for evaluating
%       this algorithm or debuging

% Javon, 10/1/2014

debug = false;

maxIter = 1e03;
threshold = 1e-04;

n = size(M{1}, 1);
nView = length(M);
d = zeros(nView, 1); % number of fetures in each view
U = zeros(n, nView);
V = cell(nView, 1);
for i = 1:nView
    d(i) = size(M{i}, 2);
    V{i} = zeros(d(i), 1);
end
% stack all M
cM = M{1};
for i = 2:nView
    cM = [cM M{i}];
end

if nargin > 6
    opt_z = optimal_z;
    opt_U = zeros(n, nView);
    opt_V = cell(nView, 1);
    for iView = 1:nView 
        [opt_U(:, iView)] = solve_uv(M{iView}, opt_z, lu(iView), lv(iView), opt_z ~= 0);
    end
    if debug
        [obj, err, pnl] = objective(M, opt_sigma, opt_z, opt_U, opt_V, lz, lv);
        fprintf('MVBC5_2: optimal obj: %2.2e, err: %2.2e, pen: %2.2e\n', obj, err, pnl);
    end
    % U = opt_U;
    % V = opt_V;
    % sigma = opt_sigma;
end

% initialization
if strcmp(view_weight, 'equal') 
    for i = 1:nView
        [U(:, i), V{i}] = sparseDecomp(M{i}, lu(i), lv(i));
    end
elseif strcmp(view_weight, 'weighted')
    fnorm2 = zeros(nView, 1);
    for i = 1:nView
        fnorm2(i) = trace(M{i} * M{i}');
    end
    [U(:, 1) V{1}] = sparseDecomp(M{1}, lz * fnorm2(1) / sum(fnorm2), lv(1));
    z = U(:, 1) ~= 0;
    for i = 2:nView
        [U(:, i) V{i}] = solve_uv(M{i}, z, lu(i), lv(i), U(:, 1));
    end
else 
    error('unknown value of parameter view_weight');
end

% given U, V, solve z
[z] = solve_z(cM, U, V, lz, maxRowClus); % given U, V, solve z
if debug
    [obj, err, pnl] = objective(M, z, U, V, lz, lu, lv);
    fprintf('MVBC5_2: Initial obj: %2.2e, err: %2.2e, pen: %2.2e, n0 in z: %d\n', obj, err, pnl, length(find(z ~= 0)));
end

% alternatively find the optimal
for iIter = 1:maxIter
    pre_z = z;
    pre_U = U;
    pre_V = V;
    % given z, solve U, V
    for iView = 1:nView % can run in parallel
        % fix v, solve u
        U(:, iView) = solve_u(M{iView}, z, V{iView}, lu(iView));
        % fix u, solve v
        V{iView} = solve_v(M{iView}, z, U(:, iView), lv(iView));
    end
    
    % given U, V, solve z
    [z] = solve_z(cM, U, V, lz, maxRowClus); % given U, V, solve z
    
    if debug
        pre_obj = obj;
        [obj, err, pnl] = objective(M, z, U, V, lz, lu, lv);
        fprintf('MVBC5_2: Iter %d -- change in u1: %2.2e, change in obj: %2.2e, n0 in z: %d, obj: %2.2e, err: %2.2e, pe: %2.2e\n', ...
            iIter, norm(pre_U(:,1) - U(:, 1)), pre_obj - obj, length(find(z ~= 0)), obj, err, pnl);
    end
    if norm(pre_U(:,1) - U(:, 1)) < threshold
        break;
    end
end
if iIter == maxIter
    warning('MVBC5_2 does not converge, try to increase the maximum iteration limit');
end

% generate clusters
rowClus = zeros(n, 1); % row clusters
rowClus(z ~= 0) = 1;
colClus = cell(nView, 1); % column clusters
for i = 1:nView
    colClus{i} = zeros(d(i), 1);
    colClus{i}(V{i} ~= 0) = 1;
end
% calculate the objective
obj = objective(M, z, U, V, lz, lu, lv);
end

function [obj err pnl] = objective(M, z, U, V, lz, lu, lv)
    obj = 0;
    err = 0;
    pnl = 0;
    for i = 1:length(M)
        mat =  M{i} - (U(:, i) .* z) * V{i}';
        err = err + trace(mat' * mat);
        pnl = pnl + norm(U(:, i), 1);
        pnl = pnl + norm(V{i}, 1);
        obj = obj + trace(mat' * mat) + lu(i) * norm(U(:, i), 1) + lv(i) * norm(V{i}, 1);
    end
    obj = obj + lz * norm(z, 1);
    pnl = pnl + lz * norm(z, 1);
end

% given z, interatively solve u and v
function [u, v] = solve_uv(M, z, lu, lv, u0)
debug = false;

maxIter = 1e04;
threshold = 1e-05;

u = u0;
v = solve_v(M, z, u, lv);
% alternatively find the optimal
for iIter = 1:maxIter
    pre_v = v;
    % given v, solve U
    [u] = solve_u(M, z, v, lu);
    % given u, solve v
    [v] = solve_v(M, z, u, lv);
    
    if debug
        fprintf('    solve_uv: Iter %d -- change in v: %2.2e\n', iIter, norm(pre_v - v));
    end
    if norm(pre_v - v) < threshold
        break;
    end
end

end

function [z] = solve_z(cM, U, V, lz, maxRowClus) 
% given U, V, solve z
n = size(cM, 1);
nView = length(V);

E = [];
for i = 1:nView
    E = [E (U(:, i) * V{i}')];
end
row_2sum_E = zeros(size(E, 1), 1);
for i = 1:size(E, 1)
    row_2sum_E(i) = E(i,:) * E(i,:)';
end;
nonZero = row_2sum_E > 0; % if ||E_i||_2^2 = 0, z_i = 0

% calculate alpha
alpha = diag(cM * E');
alpha(~nonZero) = 0;
alpha(nonZero) = alpha(nonZero) ./ row_2sum_E(nonZero);
alpha_abs = abs(alpha);

% calculate beta
beta = lz * ones(n, 1);
beta(~nonZero) = 0;
beta(nonZero) = beta(nonZero) ./ (2 * row_2sum_E(nonZero));

% compute z
z = zeros(n, 1);
nonZero = alpha_abs > beta;
z(nonZero) = sign(alpha(nonZero)) .* (alpha_abs(nonZero) - beta(nonZero));

% enforce the constraint of maximum non-zero in z 
if (length(find(nonZero)) > maxRowClus)
    % keep only the largest maxRowClus number of none zero in z
    [Y, I] = sort(abs(nonZero), 'descend');
    z(I((maxRowClus + 1) : end)) = 0;
end
end

function [u] = solve_u(M, z, v, lu) 
% given z, v, solve u
v_2norm_2 = v' * v;

if v_2norm_2 == 0
    u = zeros(size(M, 1), 1);
    return;
end

% calculate alpha, beta
alpha = zeros(length(z), 1);
beta = zeros(length(z), 1);
nonZero = z ~= 0;
mv = M * v;
alpha(nonZero) = mv(nonZero) ./ (z(nonZero) * v_2norm_2);
beta(nonZero) = lu ./ (2 * (z(nonZero).^2) * v_2norm_2);
alpha_abs = abs(alpha);

% calculate u
u = zeros(length(z), 1);
nonZero = alpha_abs > beta;
u(nonZero) = sign(alpha(nonZero)) .* (alpha_abs(nonZero) - beta(nonZero));
end

function [v] = solve_v(M, z, u, lv)
zu = z .* u;
zu_2sum = zu' * zu;

alpha = M' * zu / zu_2sum;
alpha_abs = abs(alpha);
beta = lv * ones(size(M, 2), 1) / (2 * zu_2sum);

% compute v
v = zeros(size(M, 2), 1);
nonZero = alpha_abs > beta;
v(nonZero) = sign(alpha(nonZero)) .* (alpha_abs(nonZero) - beta(nonZero));
end

function [u, v] = sparseDecomp(M, lu, lv)
[UI SI VI] = svd(M,'econ');
u = UI(:, 1) * sqrt(SI(1, 1));
v = VI(:, 1) * sqrt(SI(1, 1));
% c = sqrt(sqrt(norm(M, 'fro')));
% u = c * ones(size(M, 1), 1);
% v = c * ones(size(M, 2), 1);

debug = false;

maxIter = 1e05;
threshold = 1e-05;

if debug
    [obj, los, pen] = sv_obj(M, u, v, lu, lv);
    fprintf('Initial, obj = %.2f, los = %.2f, pen = %.2f\n', obj, los, pen);
end

% alternatively find the optimal
for iIter = 1:maxIter
    pre_u = u;
    % given u, solve v
    [v] = solve_v(M, ones(length(u), 1), u, lv);
    if debug
        [obj, los, pen] = sv_obj(M, u, v, lu, lv);
        fprintf('After optimizing v, obj = %.2f, los = %.2f, pen = %.2f\n', obj, los, pen);
    end
    % given v, solve U
    [u] = solve_v(M', ones(length(v), 1), v, lu);
    if debug
        [obj, los, pen] = sv_obj(M, u, v, lu, lv);
        fprintf('After optimizing u, obj = %.2f, los = %.2f, pen = %.2f\n', obj, los, pen);
    end
    
    if debug
        fprintf('Iter %d -- change in u: %2.2e, none zero in u: %d, none zero in v: %d\n', ...
            iIter, norm(pre_u - u), length(find(u ~= 0)), length(find(v ~= 0)));
    end
    if norm(pre_u - u) < threshold
        break;
    end
end

end

function [obj, los, pen] = sv_obj(M, u, v, lu, lv)
dist = M - u * v';
los = trace(dist' * dist);
pen = lu * norm(u, 1) + lv * norm(v, 1);
obj = los + pen;
end