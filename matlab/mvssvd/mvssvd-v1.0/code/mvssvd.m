function [rowClus, colClus, z, U, V, obj] = mvssvd(M, lz, lv, exp_z) 

% Optimization problem:
% \sum_i(||M_i - \sigma_i * (u_i \dot z) * v_i^T||_F^2) + \lambda_z * ||\sigma_1 * z||_1
% + \sum_i(\lambda_v_i * ||\sigma_i * v_i||_1)
% s.t. ||v_i||_2 = 1, ||z \dot u_i||_2 = 1

% Reference
% Multi-view Sparse Co-clustering via Proximal Alternating Linearized Minimization
% Jiangwen Sun, Jin Lu, Tiangyang and Jinbo Bi
% In the Proceedings of The 32nd International Conference on Machine Learning (ICML), 2015

% Arguments:
%   M - the list of data
%   lz - lambda on z
%   lv - lambda on v
%   exp_z - expected (true) z, if known (for example in simulation), used
%           to calculate the optimal objective value

debug = false;

maxIter = 1e05;
threshold = 1e-05;

n = size(M{1}, 1);
nView = length(M);
d = zeros(nView, 1);
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

if nargin > 3
    opt_z = exp_z;
    opt_U = zeros(n, nView);
    opt_V = cell(nView, 1);
    opt_sigma = zeros(nView, 1);
    for iView = 1:nView 
        [opt_U(:, iView), opt_V{iView}, opt_sigma(iView)] = solve_uv(M{iView}, opt_z, lv(iView), opt_z ~= 0);
    end
    if debug
        [obj, err, pnl] = objective(M, opt_sigma, opt_z, opt_U, opt_V, lz, lv);
        fprintf('MVBC3: optimal objective: %2.2e, err: %2.2e, penalty: %2.2e\n', obj, err, pnl);
    end
    % U = opt_U;
    % V = opt_V;
    % sigma = opt_sigma;
end

% initialization
fnorm2 = zeros(nView, 1);
for i = 1:nView
    fnorm2(i) = trace(M{i} * M{i}');
end
sigma = zeros(nView, 1);
[U(:, 1) V{1} sigma(1)] = ssvd(M{1}, lz * fnorm2(1) / sum(fnorm2), lv(1));
% [z V{1} sigma(1)] = ssvd(M{1}, 0.12, 0.9);
z = U(:, 1) ~= 0;
for i = 2:nView
    [U(:, i) V{i} sigma(i)] = solve_uv(M{i}, z, lv(i), U(:, 1));
end

% given U, V, solve z
[z U sigma] = solve_z(cM, U, V, sigma, lz); % given U, V, solve z
if debug
    [obj, err, pnl] = objective(M, sigma, z, U, V, lz, lv);
    fprintf('MVBC3: Initial objective: %2.2e, err: %2.2e, penalty: %2.2e, nonzero in z: %d\n', obj, err, pnl, length(find(z ~= 0)));
end

% alternatively find the optimal
for iIter = 1:maxIter
    pre_z = z;
    pre_U = U;
    pre_V = V;
    % given z, V, solve U
    for iView = 1:nView % can run in parallel
        [U(:, iView), V{iView}, sigma(iView)] = solve_uv(M{iView}, z, lv(iView), U(:, iView));        
    end
    
    % given U, V, solve z
    [z U sigma] = solve_z(cM, U, V, sigma, lz); % given U, V, solve z
    
    if debug
        [obj, err, pnl] = objective(M, sigma, z, U, V, lz, lv);
        fprintf('MVBC3: Iter %d -- change in u1: %2.2e, nonzero in z: %d, objective: %2.2e, err: %2.2e, penalty: %2.2e\n', iIter, norm(pre_U(:,1) - U(:, 1)), length(find(z ~= 0)), obj, err, pnl);
    end
    if norm(pre_U(:,1) - U(:, 1)) < threshold
        break;
    end
end
if iIter == maxIter
    warning('MVBC3 does not converge, try to increase the maximum iteration limit');
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
obj = objective(M, sigma, z, U, V, lz, lv);
end

function [obj err pnl] = objective(M, sigma, z, U, V, lz, lv)
    obj = 0;
    err = 0;
    pnl = 0;
    for i = 1:length(M)
        mat =  M{i} - sigma(i) * (U(:, i) .* z) * V{i}';
        err = err + trace(mat' * mat);
        pnl = pnl + norm(sigma(i) * V{i}, 1);
        obj = obj + trace(mat' * mat) + lv(i) * norm(sigma(i) * V{i}, 1);
    end
    obj = obj + lz * norm(sigma(1) * z, 1);
    pnl = pnl + lz * norm(sigma(1) * z, 1);
end

% ssvd, singular decompose M with both u and v being regularized
function [u, v, sigma] = ssvd(M, lu, lv)
[UI SI VI] = svd(M,'econ');
u = UI(:, 1);
v = VI(:, 1);

debug = false;

maxIter = 1e05;
threshold = 1e-05;

% alternatively find the optimal
for iIter = 1:maxIter
    pre_u = u;
    % given u, solve v
    [v] = solve_v(M, ones(length(u), 1), u, lv);
    % given v, solve U
    [u, sigma] = solve_v(M', ones(length(v), 1), v, lu);
    
    if debug
        fprintf('    SSVD: Iter %d -- change in u: %2.2e\n', iIter, norm(pre_u - u));
    end
    if norm(pre_u - u) < threshold
        break;
    end
end

end

% given z, interatively solve u and v
function [u, v, sigma] = solve_uv(M, z, lv, u0)
debug = false;

maxIter = 1e04;
threshold = 1e-05;

u = u0;
v = solve_v(M, z, u, lv);
% alternatively find the optimal
for iIter = 1:maxIter
    pre_v = v;
    % given v, solve U
    [u] = solve_u(M, z, v);
    % given u, solve v
    [v sigma] = solve_v(M, z, u, lv);
    
    if debug
        fprintf('    solve_uv: Iter %d -- change in v: %2.2e\n', iIter, norm(pre_v - v));
    end
    if norm(pre_v - v) < threshold
        break;
    end
end

end

function [z U sigma] = solve_z(cM, U, V, sigma, lz) 
% given U, V, solve z
n = size(cM, 1);
nView = length(V);

if sigma(1) == 0
    z = zeros(n, 1);
    return;
end
coef = sigma / sigma(1);
U = U * diag(coef); % for minimize ||sigma_1 * z||_1
E = [];
for i = 1:nView
    E = [E (U(:, i) * V{i}')];
end
row_2sum_E = diag(E * E'); 
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

% convert z to non-zero indicator vector
U = diag(z) * U;
z = z ~= 0;
% normalize to get sigmas
% normalize z.* U_i to unit 1
sigma = sqrt(diag(U' * U));
U = U * diag(1 ./ sigma);

end

function [u sigma] = solve_u(M, z, v) 
% given z, v, solve u
if v' * v == 0
    u = zeros(size(M, 1), 1);
    return;
end

nonZero = z ~= 0;
u = M * v;
u(~nonZero) = 0;
u(nonZero) = u(nonZero) ./ (v' * v * z(nonZero));

% normalize u to 1
sigma = norm(u);
if sigma ~= 0
    u = u ./sigma;
end
end

function [v, sigma] = solve_v(M, z, u, lv)
zu = z .* u;
zu_2sum = zu' * zu;

alpha = M' * zu / zu_2sum;
alpha_abs = abs(alpha);
beta = lv * ones(size(M, 2), 1) / (2 * zu_2sum);

% compute v
v = zeros(size(M, 2), 1);
nonZero = alpha_abs > beta;
v(nonZero) = sign(alpha(nonZero)) .* (alpha_abs(nonZero) - beta(nonZero));

% normalize v to get sigma
sigma = norm(v);
if sigma ~= 0
    v = v ./ sigma;
end
end
