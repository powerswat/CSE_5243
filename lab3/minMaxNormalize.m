function [tot_mat] = minMaxNormalize(tot_mat)

% Calculate the maximum/minimum value for each column
max_col = max(tot_mat)';
min_col = min(tot_mat)';

min_mat = repmat(min_col, [1, size(tot_mat,1)])';
max_mat = repmat(max_col, [1, size(tot_mat,1)])';

tot_mat = bsxfun(@minus, tot_mat, min_mat) ./ bsxfun(@minus, max_mat, min_mat);

end 