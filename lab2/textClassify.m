function textClassify

rng(22);

% Read the existing feature vectors
baseDir = 'C:\Temp\CSE_5243\';
isBodyFeat = 0;
if exist([baseDir, 'bdyFeat_fin.mat'], 'file')
    disp('Load existing body features');
    load([baseDir, 'bdyFeat_fin.mat']);
    isBodyFeat = 1;
end
isTpcFeat = 0;
if exist([baseDir, 'tpcFeat_fin.mat'], 'file')
    disp('Load existing topic features');
    load([baseDir, 'tpcFeat_fin.mat']);
    isTpcFeat = 1;
end
isTFIDF = 0;
if exist([baseDir, 'TFIDF_fin.mat'], 'file')
    disp('Load existing TFIDF');
    load([baseDir, 'TFIDF_fin.mat']);
    isTFIDF = 1;
end
isPlc = 0;
if exist([baseDir, 'plcFeat_fin.mat'], 'file')
    disp('Load existing place feature');
    load([baseDir, 'plcFeat_fin.mat']);
    isPlc = 1;
end

% Remove redundant features and sort it based on the decreasing order of
% TFIDF value
tpcLen = length(tpcVectLabel);
for i=1:tpcLen
    dupIdx = find(strcmp(bdyVectLabel(:,1), tpcVectLabel{i}));
    if ~isempty(dupIdx)
        bdyVectLabel(dupIdx,:) = [];
        bdyFeatMat(:,dupIdx) = [];
    end
end
[~,rnd_idx] = sort(cell2mat(bdyVectLabel(:,2)), 'descend');
bdyVectLabel = bdyVectLabel(rnd_idx,:);
bdyFeatMat = bdyFeatMat(:,rnd_idx);

%% Run the first classifier KNN like
% http://www.mathworks.com/help/stats/classification-using-nearest-neighbors.html

% Split the dataset into a training and a testing dataset
[tpcRow, ~] = find(tpcFeatMat);
tpcRow = unique(tpcRow);
tpcIdcs = [1:length(tpcFeatMat)]';
[bdyRow, ~] = find(bdyFeatMat);
bdyRow = unique(bdyRow);
bdyIdcs = [1:length(bdyFeatMat)]';
trnDocIds = intersect(tpcRow, bdyRow);
X = bdyFeatMat(trnDocIds,:);
T = tpcFeatMat(trnDocIds,:);
docIDLoc = zeros(length(X),1);
rnd_idx = randperm(length(X))';
docIDLoc = trnDocIds(rnd_idx);
X = X(rnd_idx,:);
T = T(rnd_idx,:);
trnEdIdx = floor(length(X)*0.7);
trnX = X(1:trnEdIdx,:);
trnDocID = docIDLoc(1:trnEdIdx);
valEdIdx = trnEdIdx + floor(length(X)*0.1);
valX = X(trnEdIdx+1:valEdIdx,:);
valDocID = docIDLoc(trnEdIdx+1:valEdIdx);
tstX = X(valEdIdx+1:length(X),:);
tstDocID = docIDLoc(valEdIdx+1:length(X));

% Compares the document at trnDocID(nghbrIDs(:,1)) and that at
% trnDocID(valDocID) based on the content of each body text
[nghbrIDs, dists] = knnsearch(trnX, valX, 'K', 5);


%% Run the second classifier Naive bayes like
% http://www.mathworks.com/help/stats/classification-naive-bayes.html



% Select topic features


a = 1;