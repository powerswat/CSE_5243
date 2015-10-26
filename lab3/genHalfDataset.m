function [tot_mat] = genHalfDataset(tot_mat)

filled_idcs = zeros(size(tot_mat,1),1);
cut_feat_num = floor(size(tot_mat,2)/80);
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
    
    if sum(filled_idcs) == (size(tot_mat,1)*0.95)
        break;
    end
end

% Remove zero column vectors in tmp_tot_max and re-assign to the orignal
% tot_mat matrix
rmv_idcs = find(sum(tmp_tot_mat)==0);
tmp_tot_mat(:,rmv_idcs) = [];
tot_mat = tmp_tot_mat;

end