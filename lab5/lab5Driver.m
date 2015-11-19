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
orig_time = toc;
disp('Calculating Jaccard similarity on the raw feature vector');
toc

%% Call the K-min hash oparation
k = [8 16 32 64 128 256]
len_k = length(k);
mses = zeros(len_k,1);
minhash_time = zeros(len_k,1);
for i=1:len_k
    tic;
    disp(['Running ', num2str(k(i)), '-minwise hash']);
    [est_sims] = kMinHash(tot_mat, tot_vec_lbl, k(i));
    mses(i) = sum(bsxfun(@minus, jac_sims, est_sims).^2);
    minhash_time(i) = toc;
    toc
end

disp('Elapsed time of Jaccard similarity calculation on K-min hash method');
for i=1:len_k
    disp(['K = ', num2str(k(i)), ' -> ', num2str(minhash_time(i)), ' seconds']);
end

rmses = sqrt(mses./length(est_sims));
disp('RMSE of each K-min hash method');
for i=1:len_k
    disp(['K = ', num2str(k(i)), ' -> ', num2str(rmses(i))]);
end

save([baseDir, 'rmse_results.mat'], 'rmses');

figure;
plot(k,minhash_time);
xlabel('No. of Minhash (K)')
ylabel('Sec.')
print([baseDir, 'elap_time'], '-dpng');

figure;
plot(k,rmses);
xlabel('No. of Minhash (K)')
ylabel('RMSE')
print([baseDir, 'rmse'], '-dpng');

end