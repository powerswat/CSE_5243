function [bdyFeatMat, tpcFeatMat, bdyVectLabel] = transformTpcBdyFeat(bdyVectLabel, tpcVectLabel, ...
                                        bdy_is, bdy_js, bdy_vals, tpc_is, tpc_js)

% Convert the topic and body lists to matrices
bdyFeatMat = zeros(max(bdy_is),length(bdyVectLabel));
for i=1:length(bdy_js)
    if exist('bdy_vals','var')
        bdyFeatMat(bdy_is(i),bdy_js(i)) = bdy_vals(i);
    else
        bdyFeatMat(bdy_is(i),bdy_js(i)) = 1;
    end
end
tpcFeatMat = zeros(max(bdy_is),length(tpcVectLabel));
for i=1:length(tpc_js)
    tpcFeatMat(tpc_is(i),tpc_js(i)) = 1;
end

% Remove redundant features and sort it based on the decreasing order of
% TFIDF value
dup_idx_map = zeros(length(bdyVectLabel(:,1)),1);
for i=1:length(tpcVectLabel);
    tmp_idx = find(strcmp(bdyVectLabel(:,1), tpcVectLabel{i}));
    dup_idx_map(tmp_idx,1) = 1;
end

dup_idcs = find(dup_idx_map);
bdyFeatMat(:,dup_idcs) = [];
bdyVectLabel(dup_idcs,:) = [];

[~,s_idx] = sort(cell2mat(bdyVectLabel(:,2)), 'descend');
bdyVectLabel = bdyVectLabel(s_idx,:);
bdyFeatMat = bdyFeatMat(:,s_idx);

end