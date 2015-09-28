function [slctd_tpc] = nvBayes(trnX, trnT, valX, valT, tpcVectLabel)

% Naive Bayes classifier
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
slctd_tpc = cell(size(valX,1),2);
for i=1:size(valX,1)
    feat_idcs = find(valX(i,:)==1)';
    for j=1:length(P_H)
        prob_x = PX_H(feat_idcs, j);
        non_zero_idx = find(prob_x>0);
        prov_vec(j) = prod(prob_x(non_zero_idx))*P_H(j);
    end
    [slctd_tpc{i,1}, slctd_tpc{i,2}] = max(prov_vec);
    slctd_tpc{i,2} = tpcVectLabel(slctd_tpc{i,2});
end
