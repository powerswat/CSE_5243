function [tot_mat, tot_vec_lbl] = genSmallerDataset(tot_mat, ...
                tot_vec_lbl, bdyVectLabel)
            
% Reduce the size of the feature vector by removing less important features
filled_idcs = zeros(size(tot_mat,1),1);
cut_feat_num = length(find(cell2mat(bdyVectLabel(:,2)) > 3.5));
tmp_tot_mat = zeros(size(tot_mat,1), size(tot_mat,2));

tmp_tot_mat(:,1:cut_feat_num) = tot_mat(:,1:cut_feat_num);
nz_idcs = find(sum(tmp_tot_mat, 2));
filled_idcs(nz_idcs) = 1;
idcs_to_fill = find(filled_idcs==0);
len_empty_idcs = length(idcs_to_fill);

for i=1:len_empty_idcs
    row_to_fill = idcs_to_fill(i);
    
    if(filled_idcs(row_to_fill)==1)
        continue;
    else
        filled_idcs(row_to_fill) = 1;
    end
    
    tmp_nz_idcs = find(tot_mat(row_to_fill,:));
    picked_idx = ceil(length(tmp_nz_idcs) * 0.5);
    nz_slct_col_idx = tmp_nz_idcs(picked_idx);
    tmp_tot_mat(:,cut_feat_num+i) = tot_mat(:,nz_slct_col_idx);
    idx_to_fill = find(tot_mat(:,nz_slct_col_idx));
    filled_idcs(idx_to_fill) = 1;
    
    if sum(filled_idcs) == size(tot_mat,1)
        break;
    end
end

zero_freq_feats = find(sum(tmp_tot_mat)==0);
tmp_tot_mat(:,zero_freq_feats) = [];
tot_mat = tmp_tot_mat;
tot_vec_lbl(zero_freq_feats) = [];

end