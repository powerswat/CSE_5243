function lab5Driver()

%% Read all the necessary input data
baseDir = 'C:\Temp\CSE_5243\';
disp('Read all the necessary input data');
tic;
[tot_mat, tot_vec_lbl] = readInput(baseDir);
toc

%% Create a baseline by calculating Jaccard similarity for every pair of documents
tic;
if ~exist([baseDir, 'jac_sims.mat'], 'file')
    jac_sims = 1-pdist(tot_mat, 'jaccard')';
    save([baseDir, 'jac_sims.mat'], 'jac_sims');
else
    load([baseDir, 'jac_sims.mat']);
end
toc

%% Call the K-min hash oparation
tic;
kMinHash(tot_mat, tot_vec_lbl);
toc

end