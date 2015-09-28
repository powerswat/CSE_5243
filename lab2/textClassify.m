function textClassify

% Read the existing feature vectors
baseDir = 'C:\Temp\CSE_5243\';
if exist([baseDir, 'bdyFeat_fin.mat'], 'file')
    disp('Load existing body features');
    load([baseDir, 'bdyFeat_fin.mat']);
end
if exist([baseDir, 'tpcFeat_fin.mat'], 'file')
    disp('Load existing topic features');
    load([baseDir, 'tpcFeat_fin.mat']);
end
if exist([baseDir, 'TFIDF_fin.mat'], 'file')
    disp('Load existing TFIDF');
    load([baseDir, 'TFIDF_fin.mat']);
end
isFullNghbr = 0;
if exist([baseDir, 'nghbr_full.mat'], 'file')
    disp('Load existing full neighborhood');
    load([baseDir, 'nghbr_full.mat']);
    isFullNghbr = 1;
end
isHalfNghbr = 0;
if exist([baseDir, 'nghbr_half.mat'], 'file')
    disp('Load existing half neighborhood');
    load([baseDir, 'nghbr_half.mat']);
    isHalfNghbr = 1;
end

%%  Transform the read datasets so it is able to be analyzed
[bdyFeatMat, tpcFeatMat, bdyVectLabel] = transformDataSet(TFnewIdx, bdyVectLabel, tpcVectLabel, ...
                                bdy_is, bdy_js, bdy_vals, tpc_is, tpc_js);  


% %% Run the first classifier KNN that is from built-in libary in Matlab
% % http://www.mathworks.com/help/stats/classification-using-nearest-neighbors.html
% 
% if isFullNghbr == 0 || isHalfNghbr == 0
%     % Split the dataset into a training and a testing dataset
%     [trnX,valX,cmplDocIds,~,trnT,valT,trnDocID,valDocID] = genLearningDataset(tpcFeatMat, bdyFeatMat, ...
%                             tpcVectLabel, bdyVectLabel, 0, 0);
%     [trnSmX,valSmX,~,smBdyVectLabel] = genLearningDataset(tpcFeatMat, bdyFeatMat, ...
%                             tpcVectLabel, bdyVectLabel, 1, cmplDocIds, TFIDF);
% 
%     % Compares the document at trnDocID(nghbrIDs(:,1)) and that at
%     % trnDocID(valDocID) based on the content of each body text
%     tic;
%     [nghbrIDs, dists] = knnsearch(trnX, valX, 'K', 5);
%     disp('[Full dataset] ');
%     toc
%     save([baseDir, 'nghbr_full.mat'], 'nghbrIDs', 'dists', 'trnT', 'valT', ...
%         'trnDocID', 'valDocID', 'trnX', 'valX', 'cmplDocIds');
% 
%     tic;
%     [nghbrSmIDs, smDists] = knnsearch(trnSmX, valSmX, 'K', 5);
%     disp('[Half dataset] ');
%     toc
%     save([baseDir, 'nghbr_half.mat'], 'nghbrSmIDs', 'smDists', 'trnSmX', 'valSmX');
% end
% 
% % Evaluate the performance (accuracy) of the KNN model on the both feature
% % vectors (FV-1 and FV-2)
% [crrct_vec_full, accu_full] = evalKNN(nghbrIDs, trnT, tpcVectLabel, valT);
% [crrct_vec_sm, accu_sm] = evalKNN(nghbrSmIDs, trnT, tpcVectLabel, valT);

%% Run the second classifier Naive Bayes

% Compute all the necessary conditional probabilities (Training)
P_H = (sum(trnT)./length(trnT))';
H_idcs = cell(length(P_H),1);
for i=1:length(P_H)
    H_idcs{i} = num2cell(find(trnT(:,i)));
end

[trnX_is,trnX_js,trnX_vals] = find(trnX);
for i=1:length(trnX_is)
    trnX(trnX_is(i),trnX_js(i)) = 1;
end
[valX_is,valX_js,valX_vals] = find(valX);
for i=1:length(valX_is)
    valX(valX_is(i),valX_js(i)) = 1;
end

PX_H = ones(size(trnX,2),size(trnT,2));
for i=1:size(PX_H,2)
    h_idx = cell2mat(H_idcs{i});
    if length(h_idx)==0
        continue;
    end
    
    x_h = trnX(h_idx,:);
    cnt_x_h = sum(x_h,1);
    PX_H(:,i) = (cnt_x_h./length(h_idx))';
    zero_idcs = find(PX_H(:,i)==0);
    PX_H(zero_idcs,i) = 0.0001;
end

prob_vec = zeros(length(P_H),1);
slctd_tpc = zeros(length(valX),2);
for i=1:size(valX,1)
    feat_idcs = find(valX(i,:)==1)';
    for j=1:length(P_H)
        prob_x = PX_H(feat_idcs, j);
        non_zero_idx = find(prob_x>0);
        prov_vec(j) = prod(prob_x(non_zero_idx))*P_H(j);
    end
    [slctd_tpc(i,1), slctd_tpc(i,2)] = max(prov_vec);
end
a = 1;