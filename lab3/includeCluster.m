function [expd_cluster, cluster] = includeCluster(tot_mat, cluster, i, ...
                                                cluster_id, eps, min_pts, order)

dist_mat = calcDist(+tot_mat(i,:), +tot_mat, order);

neighbors = find(dist_mat < eps);

if size(neighbors, 2) < min_pts
    cluster(i) = -1;
    expd_cluster = 0;
else
    while ~isempty(neighbors)
        curr_pt = neighbors(1);
        dist_mat = calcDist(+tot_mat(curr_pt,:), +tot_mat, order);
        res = find(dist_mat <= eps);
        
        if length(res) >= min_pts
            cluster(curr_pt) = cluster_id;
            res_unclsfd = res(find(cluster(res)==0));
            res_noise = res(find(cluster(res)==-1));
            cluster([res_unclsfd, res_noise]) = -2;
            neighbors = union(neighbors, res_unclsfd);
        end
        
        neighbors = neighbors(2:size(neighbors,2));
    end
    
    expd_cluster = 1;
end

end

