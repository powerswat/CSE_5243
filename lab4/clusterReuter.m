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

euc_hist_map = zeros(max(Euc_DB_cluster),1);
man_hist_map = zeros(max(Man_DB_cluster),1);

for i=1:length(Euc_DB_cluster)
    if Euc_DB_cluster(i) > 0
        euc_hist_map(Euc_DB_cluster(i)) = euc_hist_map(Euc_DB_cluster(i)) + 1;
    end
end
euc_num_cluster = length(find(euc_hist_map>1));

for i=1:length(Man_DB_cluster)
    if Man_DB_cluster(i) > 0
        man_hist_map(Man_DB_cluster(i)) = man_hist_map(Man_DB_cluster(i)) + 1;
    end
end
man_num_cluster = length(find(man_hist_map>4));


%% Run Kmeans Clustering
tic;
disp('Execute Kmeans Clustering (Euclidean Distance)');
[Euc_K_cluster] = KmeansCluster(tot_mat, euc_num_cluster, 2);
toc
tic;
disp('Execute Kmeans Clustering (Manhattan Distance)');
[Man_K_cluster] = KmeansCluster(tot_mat, man_num_cluster, 1);
toc

%% Analyze results
tic;
disp('Evaluate Clustering results');
[Euc_K_sil_score, Euc_K_entropy] = analyzeResult(tot_mat, Euc_K_cluster)
[Man_K_sil_score, Man_K_entropy] = analyzeResult(tot_mat, Man_K_cluster)
[Euc_K_sil_score, Euc_K_entropy] = analyzeResult(tot_mat, Euc_K_cluster)
[Man_K_sil_score, Man_K_entropy] = analyzeResult(tot_mat, Man_K_cluster)
toc
 
end