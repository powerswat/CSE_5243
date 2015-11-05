function [Euc_DB_sil_score, Euc_DB_entropy, euc_DB_var, euc_db_time, euc_db_num_clust, ...
    Man_DB_sil_score, Man_DB_entropy, man_DB_var, man_db_time, man_db_num_clust, ...
    Euc_K_sil_score, Euc_K_entropy, euc_K_var, euc_k_time, euc_DB_num_cluster, ...
    Man_K_sil_score, Man_K_entropy, man_K_var, man_k_time, man_DB_num_cluster] ...
        = clusterReuter(min_pts, eps, tot_mat, rest_tot_mat)

%% Run DBScan Clustering
tic;

% DBScan clustering using Euclidean distance
disp('Execute DBScan Clustering (Euclidean Distance)');
Euc_DB_cluster = DBScanning(tot_mat, min_pts, eps, 2);

% Assigning the remaining points into the appropriate cluster
[Euc_DB_cluster] = AssignRest(Euc_DB_cluster, tot_mat, rest_tot_mat, 2);

euc_db_num_clust = max(Euc_DB_cluster);

toc
euc_db_time = toc;

% DBScan clustering using Manhattan distance
tic;
disp('Execute DBScan Clustering (Manhattan Distance)');
Man_DB_cluster = DBScanning(tot_mat, min_pts, eps, 1);

% Assigning the remaining points into the appropriate cluster
[Man_DB_cluster] = AssignRest(Man_DB_cluster, tot_mat, rest_tot_mat, 2);

man_db_num_clust = max(Man_DB_cluster);

toc
man_db_time = toc;

% Generate some statistics for the clustering results
[euc_DB_hist_map, euc_DB_num_cluster, euc_DB_var] ...
                = checkClusterDistr(Euc_DB_cluster, 1);
[man_DB_hist_map, man_DB_num_cluster, man_DB_var] ...
                = checkClusterDistr(Man_DB_cluster, 4);


%% Run Kmeans Clustering
tic;

% Kmeans clustering using Euclidean distance
disp('Execute Kmeans Clustering (Euclidean Distance)');
[Euc_K_cluster] = KmeansCluster(tot_mat, euc_DB_num_cluster, 2);

% Assigning the remaining points into the appropriate cluster
[Euc_K_cluster] = AssignRest(Euc_K_cluster, tot_mat, rest_tot_mat, 2);
toc
euc_k_time = toc;

% Kmeans clustering using Manhattan distance
tic;
disp('Execute Kmeans Clustering (Manhattan Distance)');
[Man_K_cluster] = KmeansCluster(tot_mat, man_DB_num_cluster, 1);

% Assigning the remaining points into the appropriate cluster
[Man_K_cluster] = AssignRest(Man_K_cluster, tot_mat, rest_tot_mat, 2);
toc
man_k_time = toc;

% Generate some statistics for the clustering results
[euc_K_hist_map, euc_K_num_cluster, euc_K_var] ...
                = checkClusterDistr(Euc_K_cluster, 1);
[man_K_hist_map, man_K_num_cluster, man_K_var] ...
                = checkClusterDistr(Man_K_cluster, 4);


%% Analyze results
tic;
disp('Evaluate Clustering results');
tot_mat = [tot_mat; rest_tot_mat];
[Euc_DB_sil_score, Euc_DB_entropy] = analyzeResult(tot_mat, Euc_DB_cluster);
euc_DB_var;
[Man_DB_sil_score, Man_DB_entropy] = analyzeResult(tot_mat, Man_DB_cluster);
man_DB_var;
[Euc_K_sil_score, Euc_K_entropy] = analyzeResult(tot_mat, Euc_K_cluster);
euc_K_var;
[Man_K_sil_score, Man_K_entropy] = analyzeResult(tot_mat, Man_K_cluster);
man_K_var;
toc
 
end


%% Generate some statistics for the clustering results
function [hist_map, num_cluster, variance] =  ...
                checkClusterDistr(cluster, min_appearance)

hist_map = zeros(max(cluster), 1);
for i=1:length(cluster)
    if cluster(i) > 0
        hist_map(cluster(i)) = hist_map(cluster(i)) + 1;
    end
end
num_cluster = length(find(hist_map > min_appearance));
variance = var(hist_map);

end