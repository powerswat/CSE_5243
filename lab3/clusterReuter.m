function clusterReuter()

%% Read all the necessary input data
baseDir = 'C:\Temp\CSE_5243\';
disp('Read all the necessary input data');
tic;
[tot_mat, tot_vec_lbl] = readInputMat(baseDir);
toc

%% Min-Max normalization
tic;
disp('Execute Min-Max Normalization');
[tot_mat] = minMaxNormalize(tot_mat);
toc

%% Run DBScan



%% Run Hierarchical Clustering
tic;
disp('Execute Hierarchical Clustering');
[H_cluster] = HierarCluster(tot_mat);
toc

% Analyze results
 
end