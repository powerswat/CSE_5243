function [crrct_vec, accu] = evalKNN(nghbrIDs, trnT, tpcVectLabel, valT)

% Evaluate the performance (accuracy) of the KNN model on the both feature
% vectors (FV-1 and FV-2)

tpcFullMap = cell(size(nghbrIDs,1), size(nghbrIDs,2));
maxNumElem = 0;
for i=1:size(nghbrIDs,1)
    elemNum = 0;
    for j=1:size(nghbrIDs,2)
        tpcIdcs = find(trnT(nghbrIDs(i,j),:));
        tpcFullMap{i,j} = tpcVectLabel(tpcIdcs);
        elemNum = elemNum + length(tpcIdcs);
    end
    maxNumElem = max(maxNumElem, elemNum);
end

tpcVec = cell(maxNumElem, 1);
crrct_vec = zeros(size(tpcFullMap,1), 1);
for i=1:size(tpcFullMap,1)
    col = 1;
    for j=1:size(tpcFullMap,2)
        if length(tpcFullMap{i,j}) > 1
            for k=1:length(tpcFullMap{i,j})
                tpcVec(i,col) = tpcFullMap{i,j}(k);
                col = col + 1;
            end
        else
            tpcVec(i,col) = tpcFullMap{i,j};
            col = col + 1;
        end
    end
    tpcDict = tpcVec(i,:)';
    keep_idx = find(~cellfun(@isempty, tpcDict));
    tpcDict = tpcDict(keep_idx);
    [tpcDict,~,ind] = unique(tpcDict(keep_idx));
    freq_data = histc(ind,1:numel(tpcDict));
    [~,m_idx] = max(freq_data);
    pred_tpc = char(tpcDict(m_idx));
%     pred_tpc = char(tpcVec(i,1));
    
    % Check whether the predicted topic is included in the set of ground
    % truth topics of the validation article
    val_vec = valT(i,:);
    tpc_idx = find(val_vec==1);
    true_tpcs = tpcVectLabel(tpc_idx);
    num_match = sum(strcmp(true_tpcs, pred_tpc));
    if num_match > 0
        crrct_vec(i) = 1;
    end
end

accu = sum(crrct_vec==1)/length(crrct_vec);

end