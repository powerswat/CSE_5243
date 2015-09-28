function [bdyFeatMat, tpcFeatMat, bdyVectLabel] = transformDataSet(TFnewIdx, bdyVectLabel, tpcVectLabel, ...
                                        bdy_is, bdy_js, bdy_vals, tpc_is, tpc_js)

% Convert the topic and body lists to matrices
bdyFeatMat = zeros(length(TFnewIdx),length(bdyVectLabel));
for i=1:length(bdy_js)
    if exist('bdy_vals','var')
        bdyFeatMat(bdy_is(i),bdy_js(i)) = bdy_vals(i);
    else
        bdyFeatMat(bdy_is(i),bdy_js(i)) = 1;
    end
end
tpcFeatMat = zeros(length(TFnewIdx),length(tpcVectLabel));
for i=1:length(tpc_js)
    tpcFeatMat(tpc_is(i),tpc_js(i)) = 1;
end

% Remove redundant features and sort it based on the decreasing order of
% TFIDF value
for i=1:length(tpcVectLabel);
    dupIdx = find(strcmp(bdyVectLabel(:,1), tpcVectLabel{i}));
    if ~isempty(dupIdx)
        bdyVectLabel(dupIdx,:) = [];
        bdyFeatMat(:,dupIdx) = [];
    end
end
[~,rnd_idx] = sort(cell2mat(bdyVectLabel(:,2)), 'descend');
bdyVectLabel = bdyVectLabel(rnd_idx,:);
bdyFeatMat = bdyFeatMat(:,rnd_idx);
% bdyVectLabel(:,13701:end) = [];
% bdyFeatMat(:,13701:end) = [];