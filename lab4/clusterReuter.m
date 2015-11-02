function clusterReuter()

%% Read all the necessary input data
baseDir = 'C:\Temp\CSE_5243\';
disp('Read all the necessary input data');
tic;
[tot_mat, tot_vec_lbl, bdyVectLabel] = readInputMat(baseDir);
toc


%% Run DBScan Clustering
tic;
disp('Execute DBScan Clustering (Euclidean Distance)');
euc_min_pts = 1;  
euc_eps = 0.05;
man_min_pts = 1;
man_eps = 0.05;

Euc_DB_cluster = DBScanning(tot_mat, euc_min_pts, euc_eps, 2);
toc
tic;
disp('Execute DBScan Clustering (Manhattan Distance)');
Man_DB_cluster = DBScanning(tot_mat, man_min_pts, man_eps, 1);
toc

% Generate some statistics for the clustering results
[euc_DB_hist_map, euc_DB_num_cluster, euc_DB_var] ...
                = checkClusterDistr(Euc_DB_cluster, 1);
[man_DB_hist_map, man_DB_num_cluster, man_DB_var] ...
                = checkClusterDistr(Man_DB_cluster, 4);


%% Run Kmeans Clustering
tic;
disp('Execute Kmeans Clustering (Euclidean Distance)');
[Euc_K_cluster] = KmeansCluster(tot_mat, euc_DB_num_cluster, 2);
toc
tic;
disp('Execute Kmeans Clustering (Manhattan Distance)');
[Man_K_cluster] = KmeansCluster(tot_mat, man_DB_num_cluster, 1);
toc

% Generate some statistics for the clustering results
[euc_K_hist_map, euc_K_num_cluster, euc_K_var] ...
                = checkClusterDistr(Euc_K_cluster, 1);
[man_K_hist_map, man_K_num_cluster, man_K_var] ...
                = checkClusterDistr(Man_K_cluster, 1);


%% Analyze results
tic;
disp('Evaluate Clustering results');
[Euc_DB_sil_score, Euc_DB_entropy] = analyzeResult(tot_mat, Euc_DB_cluster)
euc_DB_var
[Man_DB_sil_score, Man_DB_entropy] = analyzeResult(tot_mat, Man_DB_cluster)
man_DB_var
[Euc_K_sil_score, Euc_K_entropy] = analyzeResult(tot_mat, Euc_K_cluster)
euc_K_var
[Man_K_sil_score, Man_K_entropy] = analyzeResult(tot_mat, Man_K_cluster)
man_K_var
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