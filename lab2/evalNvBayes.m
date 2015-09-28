function [crrct_vec, accu] = evalNvBayes(slctd_tpc, valT, tpcVectLabel)

% % Evaluate the performance (accuracy) of the Naive Bayes model on the both feature
% vectors (FV-1 and FV-2)

crrct_vec = zeros(size(valT,1), 1);
for i=1:length(valT) 
    val_vec = valT(i,:);
    tpc_idx = find(val_vec==1);
    true_tpcs = tpcVectLabel(tpc_idx);
    num_match = sum(strcmp(true_tpcs, slctd_tpc{i,2}));
    if num_match > 0
        crrct_vec(i) = 1;
    end
end

accu = sum(crrct_vec==1)/length(crrct_vec);

end