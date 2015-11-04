function [tot_mat, tot_vec_lbl, rest_tot_mat, rest_vec_lbl] = genHalfDataset(tot_mat, ...
                tot_vec_lbl, bdyVectLabel)
            
rng(28);
% rand_idx = randperm(size(tot_mat,1), 3000);
rand_idx = randperm(size(tot_mat,1));
rand_idx = rand_idx(1:300);
rest_tot_mat = tot_mat;
rest_tot_mat(rand_idx,:) = [];
rest_vec_lbl = tot_vec_lbl;
tot_mat = tot_mat(rand_idx,:);

filled_idcs = zeros(size(tot_mat,1),1);
cut_feat_num = floor(size(tot_mat,2));
tmp_tot_mat = zeros(size(tot_mat,1), size(tot_mat,2));

tmp_tot_mat(:,1:cut_feat_num) = tot_mat(:,1:cut_feat_num);
nz_idcs = find(sum(tmp_tot_mat, 2));
filled_idcs(nz_idcs) = 1;
add_filled_idx = zeros(size(tot_mat,2),1);
for i=1:size(tot_mat,1)
    
    if filled_idcs(i) == 1
        continue;
    end
    
    nz_slct_idx = min(find(tot_mat(i,:)));
    if(add_filled_idx(nz_slct_idx)==1)
        continue;
    else
        add_filled_idx(i) = 1;
    end
    tmp_tot_mat(:,cut_feat_num+i) = tot_mat(:,nz_slct_idx);
    idx_to_fill = find(tot_mat(:,nz_slct_idx));
    filled_idcs(idx_to_fill) = 1;
    
    if sum(filled_idcs) == size(tot_mat,1)
        break;
    end
end

% Muitlply TF-IDF values to each column
tpc_weight_mat = repmat(max(cell2mat(bdyVectLabel(:,2))), ...
            size(tmp_tot_mat,1), size(tot_vec_lbl,1) - size(bdyVectLabel,1));
arr_for_mult = [tpc_weight_mat, cell2mat(repmat(bdyVectLabel(:,2)', size(tmp_tot_mat,1),1))];
tmp_tot_mat = bsxfun(@times, tmp_tot_mat, arr_for_mult);

% Remove zero column vectors in tmp_tot_max and re-assign to the orignal
% tot_mat matrix
rmv_idcs = find(sum(tmp_tot_mat)<100);
rest_tot_mat(:,rmv_idcs) = [];
tmp_tot_mat(:,rmv_idcs) = [];
rest_vec_lbl(rmv_idcs) = [];
tot_vec_lbl(rmv_idcs) = [];
tot_mat = tmp_tot_mat;

end