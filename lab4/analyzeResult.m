function [sil_score, entropy] = analyzeResult(tot_mat, cluster)

%% Silhouette Coefficient
sil_scores = zeros(max(cluster), 1);
[num_cluster, ~] = size(sil_scores);

for i=1:num_cluster
    
    % Calculate the intra-distance of each cluster
    curr_cluster_idcs = find(cluster == i);
    if isempty(curr_cluster_idcs)
        sil_scores(i) = 0;
        continue;
    end
    dist_calc_template = tot_mat(curr_cluster_idcs,:);
    
    num_intra_docs = length(curr_cluster_idcs);
%     picked_idx = randperm(num_intra_docs, 1);
    picked_idx = randperm(num_intra_docs);
    picked_idx = picked_idx(1);
    
    if num_intra_docs == 1
        sil_scores(i) = 0;
        continue;
    else
        picked_others_dists = calcDist(dist_calc_template(picked_idx,:), dist_calc_template, 2);
        intra_dist = mean(picked_others_dists);
    end
    
    % Calculate average inter-distance between the selected point and the
    % points in the other clusters
    intra_pt = dist_calc_template(picked_idx, :);
    other_cluster_dists = zeros(num_cluster, 1);
    for j=1:num_cluster
        if i == j
            other_cluster_dists(i) = 0;
            continue;
        end
        curr_cluster_idcs = find(cluster == j);
        if isempty(curr_cluster_idcs)
            other_cluster_dists(j) = 0;
        else
            dist_calc_template = tot_mat(curr_cluster_idcs,:);
            other_cluster_dists(j) = mean(calcDist(intra_pt, dist_calc_template, 2));        
        end
    end
    inter_dist = mean(other_cluster_dists);
    
    sil_scores(i) = 1 - intra_dist / inter_dist;
end

% Calculate the final (= average) Silhouette score
sil_score = mean(sil_scores);

%% Information Entropy

entropies = zeros(num_cluster, 1);
num_rows = size(tot_mat, 1);
for i=1:num_cluster
    num_clst_dps = size(find(cluster == i),1);
    if num_clst_dps == 0
        entropies(i) = 0;
    else
        entropies(i) =  (num_clst_dps/num_rows) * (-(num_clst_dps / num_rows)*log2(num_clst_dps / num_rows) ...
                    -((num_rows - num_clst_dps) / num_rows)*log2((num_rows - num_clst_dps) / num_rows));
    end
end

entropy = sum(entropies);

end