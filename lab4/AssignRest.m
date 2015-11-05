function [tot_cluster] = AssignRest(cluster, tot_mat, rest_tot_mat, dist_order)

[num_rest_rows, num_feat] = size(rest_tot_mat);
rest_cluster = zeros(num_rest_rows, 1);

for i=1:num_rest_rows
   dist_mat = calcDist(rest_tot_mat(i,:), tot_mat, dist_order)';
   [~, min_idx] = min(dist_mat);
   rest_cluster(i) = cluster(min_idx);
   if mod(i,1000)==0
       disp(['Assigned rest docs: ', num2str(i), ' / ', num2str(num_rest_rows)]);
   end
end

tot_cluster = [cluster; rest_cluster];

end