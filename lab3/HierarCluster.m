function [H_cluster] = HierarCluster(tot_mat)

% Generate a similarity matrix
sim_mat = zeros(size(tot_mat,1));

% Fill the Jaccard similarity information for each cell in the matrix
nz_idcs = cell(size(tot_mat,1),1);
for i = 1:length(sim_mat)
    nz_idcs{i} = find(tot_mat(i,:));
end

tic;
for i = 1:length(sim_mat)
    for j = i:length(sim_mat)
        if i==j
            sim_mat(i,j) = 1;
            continue;
        end
        nz_idx_i = nz_idcs{i};
        nz_idx_j = nz_idcs{j};
        intsct_ij = intersect(nz_idx_i, nz_idx_j);
        sim_mat(i,j) = length(intsct_ij) / ...
                (length(nz_idx_i) + length(nz_idx_j) - length(intsct_ij));
    end
%     if i==100
%         toc
        disp(['Similarity Matrix Iteration: ', num2str(i)]);
%         a = 1;
%     end
    
end


H_cluster = 0;

end