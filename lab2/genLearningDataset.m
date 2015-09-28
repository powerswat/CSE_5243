function [trnX, valX, cmplDocIds, bdyVectLabel, trnT, valT, trnDocID, valDocID] = ...
                genLearningDataset(tpcFeatMat, bdyFeatMat, ...
                    tpcVectLabel, bdyVectLabel, isSmallSet, cmplDocIds, TFIDF)

if ~exist('TFIDF', 'var') || isempty(TFIDF), TFIDF = 0; end

rng(22);

[tpcRow, ~] = find(tpcFeatMat);
tpcRow = unique(tpcRow);
[bdyRow, ~] = find(bdyFeatMat);
bdyRow = unique(bdyRow);

% Pick 512 features if it is asked to make the smaller set
if isSmallSet == 1
    isFilledInFullMat = zeros(length(bdyFeatMat),1);
    isFilledInFullMat(cmplDocIds) = 1;
    isFilledInSubMat = zeros(length(bdyFeatMat),1);
    candBdyFeatMat = zeros(size(bdyFeatMat,1), size(bdyFeatMat,2));
    i = 1;
    while length(find(isFilledInSubMat(cmplDocIds)==1)) < length(cmplDocIds) ...
            && i <= 6900
        candBdyFeatMat(:,i) = bdyFeatMat(:,i);
        isFilledInSubMat(find(bdyFeatMat(:,i))) = 1;
%          plot(isFilledInSubMat);
%          drawnow;
        i = i + 1;
    end
    
    for j=1:size(tpcFeatMat,2)
        curTpcDocIdx = find(tpcFeatMat(:,j)==1);
        subBdyFeatMat = candBdyFeatMat(curTpcDocIdx,:);
        if sum(sum(subBdyFeatMat)) == 0
            tmpBdyFeatMat = bdyFeatMat(curTpcDocIdx,:);
            [~,cols] = find(tmpBdyFeatMat);
            if isempty(cols)
                continue;
            end
            i = i + 1;
            candBdyFeatMat(:,i) = bdyFeatMat(:,cols(1));
            bdyVectLabel(i) = bdyVectLabel(cols(1));
            isFilledInSubMat(find(bdyFeatMat(:,cols(1)))) = 1;
            if (length(find(isFilledInSubMat(cmplDocIds)==1)) >= length(cmplDocIds))
                break;
            end
        end
    end
    
    % cmplDocIds is already provided using a parameter in this case
    X = candBdyFeatMat(cmplDocIds,1:i);
else
    cmplDocIds = intersect(tpcRow, bdyRow);
    X = bdyFeatMat(cmplDocIds,:);
end

T = tpcFeatMat(cmplDocIds,:);
docIDLoc = zeros(length(X),1);
rnd_idx = randperm(size(X,1))';
docIDLoc = cmplDocIds(rnd_idx);
X = X(rnd_idx,:);
T = T(rnd_idx,:);

trnEdIdx = floor(size(X,1)*0.8);
trnX = X(1:trnEdIdx,:);
trnT = T(1:trnEdIdx,:);
trnDocID = docIDLoc(1:trnEdIdx);

valX = X(trnEdIdx+1:size(X,1),:);
valT = T(trnEdIdx+1:size(X,1),:);
valDocID = docIDLoc(trnEdIdx+1:size(X,1));

a = 1;
