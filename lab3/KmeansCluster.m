function [Euc_K_cluster, Man_K_cluster] = KmeansCluster(tot_mat, euc_num_centroids, ...
                                                man_num_centroids)

a = 1;

% Determine initial centroids
[euc_centroids, euc_dist_mat] = FindInitCentroids(tot_mat, euc_num_centroids);
[man_centroids, man_dist_mat] = FindInitCentroids(tot_mat, man_num_centroids);

end