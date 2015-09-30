function [off_eff, on_eff, accu] = textClassify(trn_ratio, TFnewIdx, ...
                bdyVectLabel, tpcVectLabel, bdy_is, bdy_js, bdy_vals, ...
                tpc_is, tpc_js, TFIDF, tpcFeatMat, bdyFeatMat)

if ~exist('trn_ratio', 'var') || isempty(trn_ratio), trn_ratio = 0.8; end

% isFullNghbr = 0;
% if exist([baseDir, 'nghbr_full.mat'], 'file')
%     disp('Load existing full neighborhood');
%     load([baseDir, 'nghbr_full.mat']);
%     isFullNghbr = 1;
% end
% isHalfNghbr = 0;
% if exist([baseDir, 'nghbr_half.mat'], 'file')
%     disp('Load existing half neighborhood');
%     load([baseDir, 'nghbr_half.mat']);
%     isHalfNghbr = 1;
% end


%% Run the first classifier KNN that is from built-in libary in Matlab
% http://www.mathworks.com/help/stats/classification-using-nearest-neighbors.html

% if isFullNghbr == 0 || isHalfNghbr == 0
    % Split the dataset into a training and a testing dataset
    [trnX,valX,cmplDocIds,~,trnT,valT,trnDocID,valDocID] = genLearningDataset(tpcFeatMat, bdyFeatMat, ...
                            tpcVectLabel, bdyVectLabel, 0, 0, 0, trn_ratio);
    [trnSmX,valSmX,~,smBdyVectLabel] = genLearningDataset(tpcFeatMat, bdyFeatMat, ...
                            tpcVectLabel, bdyVectLabel, 1, cmplDocIds, TFIDF, trn_ratio);

    % Compares the document at trnDocID(nghbrIDs(:,1)) and that at
    % trnDocID(valDocID) based on the content of each body text
    tic;
    [nghbrIDs, dists] = knnsearch(trnX, valX, 'K', 5);
    disp('[Offline cost of the full dataset by KNN] ');
    toc
    KNN_off_full_effcn = toc;
%     save([baseDir, 'nghbr_full.mat'], 'nghbrIDs', 'dists', 'trnT', 'valT', ...
%         'trnDocID', 'valDocID', 'trnX', 'valX', 'cmplDocIds');

    tic;
    [nghbrSmIDs, smDists] = knnsearch(trnSmX, valSmX, 'K', 5);
    disp('[Offline cost of the half dataset by KNN] ');
    toc
    KNN_off_half_effcn = toc;
%     save([baseDir, 'nghbr_half.mat'], 'nghbrSmIDs', 'smDists', 'trnSmX', 'valSmX');
% end

% Evaluate the performance (accuracy) of the KNN model on the both feature
% vectors (FV-1 and FV-2)
tic;
[crrct_KNN_vec_full, accu_KNN_full] = evalKNN(nghbrIDs, trnT, tpcVectLabel, valT);
disp('[Online cost of the full dataset by KNN] ');
toc
KNN_on_full_effcn = toc;
disp(['KNN accuracy on the full dataset: ' num2str(accu_KNN_full)]);

tic;
[crrct_KNN_vec_sm, accu_KNN_sm] = evalKNN(nghbrSmIDs, trnT, tpcVectLabel, valT);
disp('[Online cost of the half dataset by KNN] ');
toc
KNN_on_half_effcn = toc;
disp(['KNN accuracy on the half dataset: ' num2str(accu_KNN_sm)]);

%% Run the second classifier Naive Bayes

% Compute all the necessary conditional probabilities (Training)
tic;
slctd_full_tpc = nvBayes(trnX, trnT, valX, valT, tpcVectLabel);
disp('[Offline cost of the full dataset by Naive Bayes] ');
toc
NB_off_full_effcn = toc;

tic;
slctd_sm_tpc = nvBayes(trnSmX, trnT, valSmX, valT, tpcVectLabel);
disp('[Offline cost of the full dataset by Naive Bayes] ');
toc
NB_off_half_effcn = toc;

% Evaluate the performance (accuracy) of the Naive Bayes model on the both feature
% vectors (FV-1 and FV-2)
tic;
[crrct_NB_vec_full, accu_NB_full] = evalNvBayes(slctd_full_tpc, valT, tpcVectLabel);
disp('[Online cost of the full dataset by Naive Bayes] ');
toc
NB_on_full_effcn = toc;
disp(['Naive Bayes accuracy on the full dataset: ' num2str(accu_NB_full)]);

tic;
[crrct_NB_vec_sm, accu_NB_sm] = evalNvBayes(slctd_sm_tpc, valT, tpcVectLabel);
disp('[Online cost of the half dataset by Naive Bayes] ');
toc
NB_on_half_effcn = toc;
disp(['Naive Bayes accuracy on the half dataset: ' num2str(accu_NB_sm)]);

off_eff = [KNN_off_full_effcn KNN_off_half_effcn NB_on_full_effcn NB_on_half_effcn];
on_eff = [KNN_on_full_effcn KNN_on_half_effcn NB_on_full_effcn NB_on_half_effcn];
accu = [accu_KNN_full accu_KNN_sm accu_NB_full accu_NB_sm];

a = 1;