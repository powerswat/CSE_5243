function [centroids, dist_mat] = FindInitCentroids(tot_mat, num_centroids);

init_idx = randperm(size(tot_mat,1), 1);
dist_mat = zeros(size(tot_mat,1), size(tot_mat,1));

for i=1:num_centroids
    
end

end