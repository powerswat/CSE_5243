function [H_cluster] = HierarCluster(tot_mat)

% Generate a similarity matrix
sim_mat = zeros(size(tot_mat,1));

% Fill the Jaccard similarity information for each cell in the matrix
% nz_maps = zeros(size(tot_mat,1),size(tot_mat,2));
for i = 1:length(sim_mat)
    nz_maps(i,:) = tot_mat(i,:)>0;
end

tic;
for i = 1:length(sim_mat)
    for j = i:length(sim_mat)
        if i==j
            sim_mat(i,j) = 1;
            continue;
        end
%         nz_idx_i = nz_idcs{i};
%         nz_idx_j = nz_idcs{j};
        intsct_ij = length(find(bsxfun(@and, nz_maps(i,:), nz_maps(j,:))));
        union_ij = length(find(bsxfun(@or, nz_maps(i,:), nz_maps(j,:))));
        sim_mat(i,j) = intsct_ij / union_ij;
    end
    if mod(i,100)==0
        toc
        disp(['Similarity Matrix Iteration: ', num2str(i)]);
        a = 1;
    end 
end

H_cluster = 0;

end
