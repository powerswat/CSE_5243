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
k = [8 16 32 64 128 256]
len_k = length(k);
mses = zeros(len_k,1);
for i=1:len_k
    tic;
    [est_sims] = kMinHash(tot_mat, tot_vec_lbl, k(i));
    mses(i) = sum(bsxfun(@minus, jac_sims, est_sims).^2);
    toc
end

rmses = sqrt(mses./length(est_sims))
save([baseDir, 'rmse_results.mat'], 'rmses');

plot(k,rmses);
xlabel('No. of Minhash (K)')
ylabel('RMSE')

end