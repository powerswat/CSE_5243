function [trnX,tstX] = genLearningDataset(tpcFeatMat, bdyFeatMat)

[tpcRow, ~] = find(tpcFeatMat);
tpcRow = unique(tpcRow);
tpcIdcs = [1:length(tpcFeatMat)]';
[bdyRow, ~] = find(bdyFeatMat);
bdyRow = unique(bdyRow);
bdyIdcs = [1:length(bdyFeatMat)]';
cmplDocIds = intersect(tpcRow, bdyRow);
X = bdyFeatMat(cmplDocIds,:);
T = tpcFeatMat(cmplDocIds,:);
docIDLoc = zeros(length(X),1);
rnd_idx = randperm(length(X))';
docIDLoc = cmplDocIds(rnd_idx);
X = X(rnd_idx,:);
T = T(rnd_idx,:);
trnEdIdx = floor(length(X)*0.8);
trnX = X(1:trnEdIdx,:);
trnDocID = docIDLoc(1:trnEdIdx);
% valEdIdx = trnEdIdx + floor(length(X)*0.1);
% valX = X(trnEdIdx+1:valEdIdx,:);
% valDocID = docIDLoc(trnEdIdx+1:valEdIdx);
tstX = X(trnEdIdx+1:length(X),:);
tstDocID = docIDLoc(trnEdIdx+1:length(X));
