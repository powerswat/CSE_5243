function clusterReuter()

%% Read all the necessary input data
baseDir = 'C:\Temp\CSE_5243\';
disp('Read all the necessary input data');
tic;
[tot_mat, tot_vec_lbl, bdyVectLabel] = readInputMat(baseDir);
toc

%% Run DBScan Clustering
% tic;
% disp('Execute DBScan Clustering');
% min_pts = 3;
% eps = 2;
% 
% Euc_DB_cluster = DBScanning(tot_mat, min_pts, eps, '');
% Man_DB_cluster = DBScanning(tot_mat, min_pts, eps, 'manhattan');
% 
% toc
% 
% euc_hist_map = zeros(max(Euc_DB_cluster),1);
% man_hist_map = zeros(max(Man_DB_cluster),1);
% 
% for i=1:length(Euc_DB_cluster)
%     euc_hist_map(Euc_DB_cluster(i)) = euc_hist_map(Euc_DB_cluster(i)) + 1;
% end
% euc_num_cluster = length(find(euc_hist_map>1));
% 
% for i=1:length(Man_DB_cluster)
%     man_hist_map(Man_DB_cluster(i)) = man_hist_map(Man_DB_cluster(i)) + 1;
% end
% man_num_cluster = length(find(man_hist_map>1));

tic;
disp('Execute DBScan Clustering');
min_pts = 3;
eps = 1;

Euc_DB_cluster = dbscan(tot_mat, min_pts, eps);
toc
Man_DB_cluster = dbscan(tot_mat, min_pts, eps);

toc


%% Run Kmeans Clustering
tic;
disp('Execute Kmeans Clustering');

% Temp code --> Erase!!
euc_num_centorids = 71;
man_num_centorids = 51;

% euc_num_centorids = euc_num_cluster;
% man_num_centorids = man_num_cluster;

[Euc_K_cluster, Man_K_cluster] = KmeansCluster(tot_mat, euc_num_centorids, ...
                                                man_num_centorids);
toc

%% Analyze results

 
end