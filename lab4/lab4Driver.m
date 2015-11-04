function lab4Driver(min_pts, eps)

if ~exist('min_pts', 'var') || isempty(min_pts), min_pts = 1; end
if ~exist('eps', 'var') || isempty(eps), eps = 0.05; end

num_exec = length(min_pts) * length(eps);
param_set = zeros(num_exec, 2);
idx = 1;
for i=1:length(min_pts)
    for j=1:length(eps)
        param_set(idx,1) = min_pts(i);
        param_set(idx,2) = eps(j);e
        idx = idx + 1;
    end
end

for i=1:num_exec
    disp(['Minimun number neighbor: ', num2str(param_set(i,1)) ...
                , ', Epsilon: ', num2str(param_set(i,2))]);
    [Euc_DB_sil_score, Euc_DB_entropy, euc_DB_var, ...
    Man_DB_sil_score, Man_DB_entropy, man_DB_var, ...
    Euc_K_sil_score, Euc_K_entropy, euc_K_var, ...
    Man_K_sil_score, Man_K_entropy, man_K_var] = clusterReuter ...
                                                (param_set(i,1), param_set(i,2))
end

end