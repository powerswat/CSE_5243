function [est_sims] = kMinHash(tot_mat, tot_vec_lbl, k)

%% Generate an input matrix (# of shingles * # of docs)
input_mat = tot_mat';
N = size(input_mat,1);

%% Linear Permutation Hashing to generate a sigmature matrix M
num_perm = k;
[num_sh, num_doc] = size(input_mat);
M = zeros(num_perm, num_doc);
for i=1:num_perm
    perm = randperm(N);
    for j=1:num_doc
        M(i,j) = min(perm(find(input_mat(:,j))));
    end
end

%% Calculate the similarities bewteen each two documents
est_sims = 1-pdist(M', 'jaccard')';

end