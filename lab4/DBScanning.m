function [cluster] = DBScanning(tot_mat, min_pts, eps, order)

doc_num = size(tot_mat,1);
cluster = zeros(doc_num,1);
cluster_id = 1;
for i=1:doc_num
    if cluster(i) == 0
        [expd_cluster, cluster] = includeCluster(tot_mat, cluster, i, ...
                                                cluster_id, eps, min_pts, order);
        if expd_cluster
            cluster_id = cluster_id + 1;
        end
    end
end

score = cluster;
core_idx = find(score > 0);
border_pts = find(cluster == -2);
for i=1:length(border_pts)
    curr_b = border_pts(i);
    neighbors = calcDist(+tot_mat(curr_b), +tot_mat(core_idx), order);
    [tmp nearest_core] = min(neighbors);
    nearest_core_idx = core_idx(nearest_core);
    cluster(curr_b) = cluster(nearest_core_idx);
end

end
